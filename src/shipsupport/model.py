"""Input data model.

All quantities are SI (N, m, Pa). The longitudinal coordinate x is measured
from the aft perpendicular (AP), positive forward. Vertical displacement is
positive downward, into the blocks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


class InputError(ValueError):
    """Invalid or inconsistent input data."""


def _require_positive(value, what):
    if value is None or not math.isfinite(value) or value <= 0:
        raise InputError(f"{what} must be a positive finite number (got {value!r})")


def _optional_positive(value, what):
    if value is not None:
        _require_positive(value, what)


def _as_array(values, what):
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        raise InputError(f"{what} must be a non-empty list of finite numbers")
    return arr


# --------------------------------------------------------------------------- #
# Weights and hull girder
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class WeightItem:
    """A weight spread as a trapezoid over [xa, xb], or a point load if xa == xb.

    The trapezoid preserves both the weight and its LCG, which is only possible
    with non-negative end intensities when the LCG lies in the middle third.
    """

    name: str
    weight: float   # N (negative allowed for deductions)
    lcg: float      # m from AP
    xa: float       # m
    xb: float       # m

    def __post_init__(self):
        if not math.isfinite(self.weight) or self.weight == 0:
            raise InputError(f"Weight item '{self.name}': weight must be non-zero and finite")
        if self.xb < self.xa:
            raise InputError(f"Weight item '{self.name}': xb ({self.xb}) < xa ({self.xa})")
        if self.is_point:
            if abs(self.lcg - self.xa) > 1e-9:
                raise InputError(f"Weight item '{self.name}': a point load needs lcg == xa == xb")
            return
        third = self.length / 3.0
        xbar = self.lcg - self.xa
        tol = 1e-9 * self.length
        if not (third - tol <= xbar <= 2.0 * third + tol):
            raise InputError(
                f"Weight item '{self.name}': LCG {self.lcg:.3f} m must lie in the middle third of "
                f"its extent [{self.xa:.3f}, {self.xb:.3f}] m for a non-negative trapezoid; "
                f"widen the extent or split the item")

    @property
    def length(self) -> float:
        return self.xb - self.xa

    @property
    def is_point(self) -> bool:
        return self.length <= 1e-9


@dataclass(frozen=True, eq=False)
class WeightCurve:
    """Piecewise-linear weight distribution q(x) ≥ 0 (N/m) through the points (x_k, q_k).

    The load is zero outside the first and last point. x must not decrease; the same x
    twice is a jump of q. Every segment is an exact trapezoid, so the curve is converted
    into one WeightItem per loaded segment (see items()) and the load model is unchanged.
    """

    name: str
    x: np.ndarray        # m from AP
    q: np.ndarray        # N/m

    def __post_init__(self):
        object.__setattr__(self, "x", np.asarray(self.x, dtype=float).ravel())
        object.__setattr__(self, "q", np.asarray(self.q, dtype=float).ravel())
        if self.x.size != self.q.size or self.x.size < 2:
            raise InputError(f"Weight curve '{self.name}': give at least two points (x, q)")
        for k, (x, q) in enumerate(zip(self.x, self.q, strict=True)):
            where = f"Weight curve '{self.name}', point {k + 1} (x = {x:g} m)"
            if not (math.isfinite(x) and math.isfinite(q)):
                raise InputError(f"{where}: x and q must be finite numbers")
            if q < 0:
                raise InputError(f"{where}: q = {q:g} N/m must not be negative")
            if k and x < self.x[k - 1]:
                raise InputError(f"{where}: x must not be smaller than at point {k} "
                                 f"(x = {self.x[k - 1]:g} m); give the points from aft to forward")
        if not self.items():
            raise InputError(f"Weight curve '{self.name}' carries no weight")

    def items(self) -> list[WeightItem]:
        """One exact trapezoidal item per segment with load: weight (q_k + q_k+1) L / 2 and
        LCG x_k + L (q_k + 2 q_k+1) / (3 (q_k + q_k+1)). Jumps and unloaded segments are skipped."""
        out = []
        n = self.x.size - 1
        for k in range(n):
            xa, xb, qa, qb = self.x[k], self.x[k + 1], self.q[k], self.q[k + 1]
            length = xb - xa
            if length <= 1e-9 or qa + qb <= 0:
                continue
            lcg = xa + length * (qa + 2.0 * qb) / (3.0 * (qa + qb))
            out.append(WeightItem(f"{self.name} {k + 1}/{n}", 0.5 * (qa + qb) * length,
                                  float(lcg), float(xa), float(xb)))
        return out


@dataclass(frozen=True)
class Section:
    """Hull girder properties, constant over [xa, xb]."""

    xa: float
    xb: float
    EI: float                          # N·m²
    kGA: float | None = None           # N, effective shear stiffness G·As (None = no shear deformation)
    z_deck: float | None = None        # m³
    z_keel: float | None = None        # m³
    shear_area: float | None = None    # m², effective shear area for τ = V / As

    def __post_init__(self):
        if self.xb <= self.xa:
            raise InputError(f"Section [{self.xa}, {self.xb}]: xb must be greater than xa")
        _require_positive(self.EI, f"Section [{self.xa}, {self.xb}] EI")
        _optional_positive(self.kGA, f"Section [{self.xa}, {self.xb}] kGA")
        _optional_positive(self.z_deck, f"Section [{self.xa}, {self.xb}] z_deck")
        _optional_positive(self.z_keel, f"Section [{self.xa}, {self.xb}] z_keel")
        _optional_positive(self.shear_area, f"Section [{self.xa}, {self.xb}] shear_area")


# --------------------------------------------------------------------------- #
# Blocks
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Layer:
    thickness: float   # m
    E: float           # Pa, modulus in the loading direction (perpendicular to grain for timber)

    def __post_init__(self):
        _require_positive(self.thickness, "Layer thickness")
        _require_positive(self.E, "Layer E")


@dataclass(frozen=True)
class BlockType:
    """A block build-up. Layers act as springs in series over the contact area.

    If crush_stress is given the block is bilinear: stiffness k up to the crush
    force, then post_crush_ratio * k.
    """

    name: str
    contact_width: float                     # m, transverse
    contact_length: float                    # m, longitudinal
    layers: tuple[Layer, ...] = ()
    stiffness_override: float | None = None  # N/m, replaces the layer calculation
    capacity: float | None = None            # N, rated block load
    allowable_pressure: float | None = None  # Pa, allowable bearing pressure
    crush_stress: float | None = None        # Pa
    post_crush_ratio: float = 0.05

    def __post_init__(self):
        _require_positive(self.contact_width, f"Block type '{self.name}' contact_width")
        _require_positive(self.contact_length, f"Block type '{self.name}' contact_length")
        if not self.layers and self.stiffness_override is None:
            raise InputError(f"Block type '{self.name}': give layers or a stiffness override")
        _optional_positive(self.stiffness_override, f"Block type '{self.name}' stiffness")
        _optional_positive(self.capacity, f"Block type '{self.name}' capacity")
        _optional_positive(self.allowable_pressure, f"Block type '{self.name}' allowable_pressure")
        _optional_positive(self.crush_stress, f"Block type '{self.name}' crush_stress")
        if not 0.0 <= self.post_crush_ratio <= 1.0:
            raise InputError(f"Block type '{self.name}': post_crush_ratio must be in [0, 1]")

    @property
    def area(self) -> float:
        return self.contact_width * self.contact_length

    @property
    def stiffness(self) -> float:
        """Initial vertical stiffness, N/m."""
        if self.stiffness_override is not None:
            return self.stiffness_override
        compliance = sum(layer.thickness / layer.E for layer in self.layers)
        return self.area / compliance

    @property
    def crush_force(self) -> float | None:
        return None if self.crush_stress is None else self.crush_stress * self.area


@dataclass(frozen=True)
class Block:
    x: float                 # m from AP
    type: BlockType
    gap: float = 0.0         # m, positive = block top below the keel reference line
    label: str = ""
    strength: float | None = None   # N, force this block can carry; None = its type's capacity

    def __post_init__(self):
        _optional_positive(self.strength, f"Block at x = {self.x} strength")

    @property
    def effective_strength(self) -> float | None:
        """Force the block can carry (N): its own strength, else its type's capacity."""
        return self.strength if self.strength is not None else self.type.capacity


# --------------------------------------------------------------------------- #
# Permissible hull girder limits
# --------------------------------------------------------------------------- #

@dataclass(frozen=True, eq=False)
class BendingEnvelope:
    """Permissible still-water bending moment; sag and hog are positive magnitudes (N·m)."""

    x: np.ndarray
    sag: np.ndarray
    hog: np.ndarray

    def __post_init__(self):
        for name in ("x", "sag", "hog"):
            object.__setattr__(self, name, _as_array(getattr(self, name), f"Bending envelope {name}"))
        if not (self.x.size == self.sag.size == self.hog.size):
            raise InputError("Bending envelope: x, sag and hog must have equal length")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Bending envelope: x must be strictly increasing")
        if np.any(self.sag <= 0) or np.any(self.hog <= 0):
            raise InputError("Bending envelope: limits must be positive magnitudes")

    def at(self, x):
        return np.interp(x, self.x, self.sag), np.interp(x, self.x, self.hog)


@dataclass(frozen=True, eq=False)
class ShearEnvelope:
    """Permissible still-water shear force magnitude (N)."""

    x: np.ndarray
    limit: np.ndarray

    def __post_init__(self):
        for name in ("x", "limit"):
            object.__setattr__(self, name, _as_array(getattr(self, name), f"Shear envelope {name}"))
        if self.x.size != self.limit.size:
            raise InputError("Shear envelope: x and limit must have equal length")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Shear envelope: x must be strictly increasing")
        if np.any(self.limit <= 0):
            raise InputError("Shear envelope: limits must be positive")

    def at(self, x):
        return np.interp(x, self.x, self.limit)


# --------------------------------------------------------------------------- #
# Welding deformation of the assembled hull
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class WeldJoint:
    """Angular distortion (kink) at an erection joint between two hull blocks."""

    x: float        # m from AP
    angle: float    # rad; positive = concave-up kink, which lifts both hull ends

    def __post_init__(self):
        if not (math.isfinite(self.x) and math.isfinite(self.angle)):
            raise InputError("Weld joint x and angle must be finite")


@dataclass(frozen=True)
class ShrinkageZone:
    """Uniform welding-shrinkage curvature over a hull block [xa, xb]."""

    xa: float
    xb: float
    curvature: float   # 1/m; positive = concave up (deck shortens more than bottom), lifts the ends

    def __post_init__(self):
        if self.xb <= self.xa:
            raise InputError(f"Shrinkage zone [{self.xa}, {self.xb}]: xb must be greater than xa")
        if not math.isfinite(self.curvature):
            raise InputError("Shrinkage zone curvature must be finite")


PROFILE_INTERPOLATIONS = ("pchip", "linear", "spline")


@dataclass(frozen=True, eq=False)
class MeasuredProfile:
    """Surveyed keel line after welding: upward lift (m) at stations x (m)."""

    x: np.ndarray
    lift: np.ndarray
    interpolation: str = "pchip"
    identify_joints: bool = False   # fit joint angles to the survey instead of interpolating it
    regularization: float = 0.0     # Tikhonov weight for the joint-angle fit

    def __post_init__(self):
        object.__setattr__(self, "x", _as_array(self.x, "Measured profile x"))
        object.__setattr__(self, "lift", _as_array(self.lift, "Measured profile lift"))
        if self.x.size != self.lift.size or self.x.size < 2:
            raise InputError("Measured profile: x and lift need the same length (at least 2 points)")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Measured profile: x must be strictly increasing")
        if self.interpolation not in PROFILE_INTERPOLATIONS:
            raise InputError(f"Measured profile: interpolation must be one of {PROFILE_INTERPOLATIONS}")
        if self.regularization < 0:
            raise InputError("Measured profile: regularization must be >= 0")


@dataclass(eq=False)
class WeldingDeformation:
    """Initial (stress-free) keel shape of the welded hull.

    The shape is the sum of joint kinks, zone curvatures and an interpolated survey.
    Only its deviation from a straight line matters: a straight-line part is a
    rigid-body tilt that the support solvers absorb.
    """

    joints: tuple[WeldJoint, ...] = ()
    zones: tuple[ShrinkageZone, ...] = ()
    measured: MeasuredProfile | None = None

    def __post_init__(self):
        self.joints = tuple(self.joints)
        self.zones = tuple(self.zones)
        if not self.joints and not self.zones and self.measured is None:
            raise InputError("Welding deformation needs joints, zones or a measured profile")
        if self.measured is not None and self.measured.identify_joints and not self.joints:
            raise InputError("Identifying joint angles from the survey needs the joint positions")


# --------------------------------------------------------------------------- #
# Ship, settings, case
# --------------------------------------------------------------------------- #

@dataclass(eq=False)
class Ship:
    name: str
    lpp: float
    weights: list[WeightItem]
    sections: list[Section]
    x_aft: float = 0.0                     # m, aft end of the hull
    x_fwd: float | None = None             # m, forward end of the hull (default Lpp)
    declared_weight: float | None = None   # N, from the docking plan, for cross-checking
    declared_lcg: float | None = None      # m
    structure_positions: np.ndarray | None = None   # m, floors / bulkheads / web frames
    permissible_bm: BendingEnvelope | None = None
    permissible_sf: ShearEnvelope | None = None
    # the weight curve of the input, if any; its segments are already among the weights
    weight_curve: WeightCurve | None = None

    def __post_init__(self):
        _require_positive(self.lpp, "Ship Lpp")
        if self.x_fwd is None:
            self.x_fwd = self.lpp
        if self.x_fwd <= self.x_aft:
            raise InputError("Ship x_fwd must be greater than x_aft")
        if not self.weights:
            raise InputError("Ship needs at least one weight item")
        if not self.sections:
            raise InputError("Ship needs at least one hull section")
        if self.total_weight <= 0:
            raise InputError("Ship total weight must be positive")
        if self.structure_positions is not None:
            self.structure_positions = np.sort(_as_array(self.structure_positions, "Structure positions"))

    @property
    def total_weight(self) -> float:
        return float(sum(w.weight for w in self.weights))

    @property
    def lcg(self) -> float:
        return float(sum(w.weight * w.lcg for w in self.weights)) / self.total_weight


@dataclass
class Settings:
    methods: tuple[str, ...] = ("rigid", "fem")
    max_element_length: float = 0.5               # m
    shear_deformation: bool = True                 # Timoshenko elements where kGA is given
    allowable_bending_stress: float | None = None  # Pa
    allowable_shear_stress: float | None = None    # Pa
    structure_tolerance: float = 0.15              # m, allowed block offset from structure
    weight_tolerance: float = 0.005                # relative, declared vs computed weight
    max_iterations: int = 200

    def __post_init__(self):
        _require_positive(self.max_element_length, "max_element_length")
        _optional_positive(self.allowable_bending_stress, "allowable_bending_stress")
        _optional_positive(self.allowable_shear_stress, "allowable_shear_stress")
        if self.structure_tolerance < 0:
            raise InputError("structure_tolerance must be >= 0")
        if self.weight_tolerance < 0:
            raise InputError("weight_tolerance must be >= 0")
        if self.max_iterations < 1:
            raise InputError("max_iterations must be >= 1")


@dataclass
class SpacingDesign:
    """Spacing design: the positions of a given number of blocks with given strengths
    (see spacing.py). The block count and the order of the blocks are kept."""

    method: str = "fem"                      # support model: fem | rigid
    target_utilisation: float = 1.0          # required max of reaction / strength
    min_spacing: float = 1.0                 # m
    max_spacing: float | None = None         # m, None = 2 × mean spacing of the given line
    x_min: float | None = None               # m, aft limit of the block line (None = first block)
    x_max: float | None = None               # m, forward limit (None = last block)
    fix_end_blocks: bool = False             # keep the first and last block where they are
    snap_to_structure: bool = True           # final positions at structure positions
    grid: float | None = None                # m, rounding step when not snapping (None = none)
    max_iterations: int = 20                 # sequential-LP iterations

    def __post_init__(self):
        if self.method not in ("rigid", "fem"):
            raise InputError("Spacing design method must be 'rigid' or 'fem'")
        _require_positive(self.target_utilisation, "Spacing design target_utilisation")
        _require_positive(self.min_spacing, "Spacing design min_spacing")
        if self.max_spacing is not None and self.max_spacing < self.min_spacing:
            raise InputError("Spacing design max_spacing must not be smaller than min_spacing")
        if self.x_min is not None and self.x_max is not None and self.x_max <= self.x_min:
            raise InputError("Spacing design x_max must be greater than x_min")
        _optional_positive(self.grid, "Spacing design grid")
        if self.max_iterations < 0:
            raise InputError("Spacing design max_iterations must be >= 0")


@dataclass(eq=False)
class Case:
    ship: Ship
    blocks: list[Block]
    settings: Settings = field(default_factory=Settings)
    deformation: WeldingDeformation | None = None
    spacing: SpacingDesign | None = None

    def __post_init__(self):
        if len(self.blocks) < 2:
            raise InputError("At least two blocks are required")
