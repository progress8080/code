"""Welding deformation of the assembled hull and the resulting block gaps.

A hull built from prefabricated blocks is not straight after welding. Two
mechanisms are modelled, both as an upward keel lift z(x):

* angular distortion φ at an erection joint x_j (a kink: the slope jumps by φ),
      z_j(x) = φ · |x − x_j| / 2
* a uniform curvature κ over a hull block [xa, xb] from differential shrinkage
  of deck and bottom (z'' = κ inside the block, straight outside),
      z(x) = κ · ψ(x),  ψ = 0 (x < xa), (x − xa)²/2 (xa ≤ x ≤ xb),
                            (xb − xa)(x − (xa + xb)/2) (x > xb)

Positive φ and κ make the keel concave up, i.e. they lift the hull ends.
Alternatively a surveyed keel line is interpolated, or used to identify the
joint angles by regularised linear least squares.

Only the deviation from a straight line matters (a straight line is a rigid
tilt absorbed by the support solvers), so the gap at block i is the lift
relative to the least-squares line through the block positions, shifted so
that the lowest block point has zero gap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from .model import MeasuredProfile, WeldingDeformation, WeldJoint


def kink_shape(x: np.ndarray, xj: float) -> np.ndarray:
    """Lift per radian of a joint kink at xj (symmetric form)."""
    return 0.5 * np.abs(np.asarray(x, dtype=float) - xj)


def zone_shape(x: np.ndarray, xa: float, xb: float) -> np.ndarray:
    """Lift per unit curvature of a zone [xa, xb]: z'' = 1 inside, 0 outside."""
    x = np.asarray(x, dtype=float)
    inside = (x - xa) ** 2 / 2.0
    after = (xb - xa) * (x - 0.5 * (xa + xb))
    return np.where(x < xa, 0.0, np.where(x <= xb, inside, after))


def interpolate_profile(profile: MeasuredProfile, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if profile.interpolation == "linear":
        # linear extrapolation with the end slopes
        y = np.interp(x, profile.x, profile.lift)
        lo, hi = x < profile.x[0], x > profile.x[-1]
        s0 = (profile.lift[1] - profile.lift[0]) / (profile.x[1] - profile.x[0])
        s1 = (profile.lift[-1] - profile.lift[-2]) / (profile.x[-1] - profile.x[-2])
        y[lo] = profile.lift[0] + s0 * (x[lo] - profile.x[0])
        y[hi] = profile.lift[-1] + s1 * (x[hi] - profile.x[-1])
        return y
    if profile.interpolation == "spline" and profile.x.size >= 3:
        return CubicSpline(profile.x, profile.lift, bc_type="natural", extrapolate=True)(x)
    return PchipInterpolator(profile.x, profile.lift, extrapolate=True)(x)


@dataclass(frozen=True)
class Identification:
    """Result of fitting joint angles to a surveyed keel line."""

    joints: tuple[WeldJoint, ...]
    rms_residual: float          # m
    max_residual: float          # m


def identify_joint_angles(profile: MeasuredProfile, joint_x, regularization: float = 0.0
                          ) -> Identification:
    """Least-squares fit  z ≈ c0 + c1·x + Σ φ_j |x − x_j|/2  with Tikhonov weight on φ.

    The straight-line terms c0, c1 absorb the unknown survey datum and tilt.
    """
    joint_x = np.asarray(joint_x, dtype=float)
    x, z = profile.x, profile.lift
    scale = max(np.ptp(x), 1.0)
    A = np.column_stack([np.ones_like(x), x / scale] +
                        [kink_shape(x, xj) for xj in joint_x])
    rhs = z.copy()
    if regularization > 0:
        reg = np.zeros((joint_x.size, A.shape[1]))
        reg[:, 2:] = np.sqrt(regularization) * np.eye(joint_x.size) * scale
        A = np.vstack([A, reg])
        rhs = np.concatenate([rhs, np.zeros(joint_x.size)])
    coef, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    fitted = A[:x.size] @ coef
    residual = z - fitted
    joints = tuple(WeldJoint(float(xj), float(phi)) for xj, phi in zip(joint_x, coef[2:]))
    return Identification(joints, float(np.sqrt(np.mean(residual**2))),
                          float(np.max(np.abs(residual))))


def effective_model(deformation: WeldingDeformation) -> tuple[WeldingDeformation, Identification | None]:
    """The deformation model actually used for the gaps.

    * survey with identify_joints: the joint angles are fitted to the survey after
      the known zone-shrinkage part has been subtracted (the survey contains it);
    * survey without identification: the survey is the total keel shape, so any
      joints and zones given as well are not added on top of it;
    * no survey: joints and zones as given.
    """
    m = deformation.measured
    if m is None:
        return deformation, None
    if not m.identify_joints:
        return WeldingDeformation(measured=m), None
    zone_part = sum((z.curvature * zone_shape(m.x, z.xa, z.xb) for z in deformation.zones),
                    np.zeros_like(m.x))
    residual = MeasuredProfile(m.x, m.lift - zone_part, m.interpolation)
    ident = identify_joint_angles(residual, [j.x for j in deformation.joints], m.regularization)
    return WeldingDeformation(joints=ident.joints, zones=deformation.zones), ident


def keel_lift(deformation: WeldingDeformation, x) -> np.ndarray:
    """Upward keel lift z(x) of the welded hull (defined up to a straight line)."""
    x = np.asarray(x, dtype=float)
    z = np.zeros_like(x)
    for j in deformation.joints:
        z += j.angle * kink_shape(x, j.x)
    for zone in deformation.zones:
        z += zone.curvature * zone_shape(x, zone.xa, zone.xb)
    if deformation.measured is not None:
        z += interpolate_profile(deformation.measured, x)
    return z


def gaps_from_lift(lift_at_blocks: np.ndarray, block_x: np.ndarray) -> np.ndarray:
    """Block gaps (m, ≥ 0): lift relative to the least-squares line through the
    block positions, shifted so the lowest point has zero gap."""
    lift = np.asarray(lift_at_blocks, dtype=float)
    x = np.asarray(block_x, dtype=float)
    if np.ptp(x) > 0:
        slope, intercept = np.polyfit(x, lift, 1)
        lift = lift - (intercept + slope * x)
    return lift - lift.min()


def lift_components(deformation: WeldingDeformation, x_curve, block_x) -> dict[str, np.ndarray]:
    """Keel lift split into its parts (joint kinks, zone shrinkage, survey), each with
    the least-squares line through the block positions removed. Their sum plus a
    constant equals detrended_lift."""
    block_x = np.asarray(block_x, dtype=float)
    x_curve = np.asarray(x_curve, dtype=float)
    parts = {}
    if deformation.joints:
        parts["joint kinks"] = WeldingDeformation(joints=deformation.joints)
    if deformation.zones:
        parts["zone shrinkage"] = WeldingDeformation(zones=deformation.zones)
    if deformation.measured is not None:
        parts["survey (interpolated)"] = WeldingDeformation(measured=deformation.measured)
    out = {}
    for name, part in parts.items():
        zb = keel_lift(part, block_x)
        slope, intercept = np.polyfit(block_x, zb, 1) if np.ptp(block_x) > 0 else (0.0, 0.0)
        out[name] = keel_lift(part, x_curve) - (intercept + slope * x_curve)
    return out


def detrended_lift(deformation: WeldingDeformation, x_curve, block_x) -> np.ndarray:
    """The same normalisation as gaps_from_lift, applied to a curve for plotting."""
    block_x = np.asarray(block_x, dtype=float)
    zb = keel_lift(deformation, block_x)
    slope, intercept = np.polyfit(block_x, zb, 1) if np.ptp(block_x) > 0 else (0.0, 0.0)
    shift = (zb - (intercept + slope * block_x)).min()
    x_curve = np.asarray(x_curve, dtype=float)
    return keel_lift(deformation, x_curve) - (intercept + slope * x_curve) - shift
