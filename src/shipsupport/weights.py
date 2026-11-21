"""Longitudinal weight distribution.

Weight items are spread as trapezoids onto a mesh whose nodes include every
item boundary, so each element carries an exactly linear load and all
integrals are exact.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from .model import InputError, WeightItem

NODE_TOL = 1e-6  # m, points closer than this are the same node


def trapezoid_ends(item: WeightItem) -> tuple[float, float]:
    """End intensities (N/m) of the trapezoid that preserves the item's weight and LCG.

    From W = L(wa + wb)/2 and W·x̄ = L²(wa + 2wb)/6.
    """
    L = item.length
    xbar = item.lcg - item.xa
    wa = 2.0 * item.weight / L * (2.0 - 3.0 * xbar / L)
    wb = 2.0 * item.weight / L * (3.0 * xbar / L - 1.0)
    return wa, wb


def build_nodes(key_points: Iterable[float], x_aft: float, x_fwd: float,
                max_element_length: float) -> np.ndarray:
    """Sorted mesh nodes containing every key point, no element longer than max_element_length."""
    if max_element_length <= 0:
        raise InputError("max_element_length must be positive")
    pts = np.asarray(list(key_points) + [x_aft, x_fwd], dtype=float)
    outside = (pts < x_aft - NODE_TOL) | (pts > x_fwd + NODE_TOL)
    if outside.any():
        raise InputError(f"x = {pts[outside][0]:.3f} m lies outside the hull extent "
                         f"[{x_aft:.3f}, {x_fwd:.3f}] m")
    pts = np.sort(np.clip(pts, x_aft, x_fwd))
    pts = pts[np.concatenate([[True], np.diff(pts) > NODE_TOL])]
    pts[0], pts[-1] = x_aft, x_fwd

    nodes = [pts[:1]]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(1, math.ceil((b - a) / max_element_length - 1e-9))
        nodes.append(a + (b - a) * np.arange(1, n + 1) / n)
    return np.concatenate(nodes)


def node_index(nodes: np.ndarray, x) -> np.ndarray:
    """Index of the node at each x; raises if x is not a mesh node."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    idx = np.clip(np.searchsorted(nodes, x), 1, len(nodes) - 1)
    idx = np.where(np.abs(x - nodes[idx - 1]) <= np.abs(nodes[idx] - x), idx - 1, idx)
    off = np.abs(nodes[idx] - x) > NODE_TOL
    if off.any():
        raise InputError(f"x = {x[off][0]:.4f} m is not a mesh node")
    return idx


@dataclass(frozen=True, eq=False)
class LoadModel:
    """Piecewise-linear distributed load plus nodal point loads (downward positive)."""

    nodes: np.ndarray
    q1: np.ndarray      # N/m at the aft end of each element
    q2: np.ndarray      # N/m at the forward end of each element
    point: np.ndarray   # N at each node

    @property
    def lengths(self) -> np.ndarray:
        return np.diff(self.nodes)

    @property
    def total(self) -> float:
        return float(np.sum(self.lengths * (self.q1 + self.q2) / 2.0) + self.point.sum())

    @property
    def moment(self) -> float:
        """First moment of the load about x = 0, N·m."""
        x1, x2, L = self.nodes[:-1], self.nodes[1:], self.lengths
        distributed = np.sum(L * (self.q1 * (2.0 * x1 + x2) + self.q2 * (x1 + 2.0 * x2)) / 6.0)
        return float(distributed + np.dot(self.point, self.nodes))

    @property
    def lcg(self) -> float:
        return self.moment / self.total

    def intensity_curve(self) -> tuple[np.ndarray, np.ndarray]:
        """(x, w) pairs per element, suitable for plotting jumps in w(x)."""
        x = np.column_stack([self.nodes[:-1], self.nodes[1:]]).ravel()
        w = np.column_stack([self.q1, self.q2]).ravel()
        return x, w


def build_load_model(items: Sequence[WeightItem], nodes: np.ndarray) -> LoadModel:
    x1, x2 = nodes[:-1], nodes[1:]
    xm = 0.5 * (x1 + x2)
    q1 = np.zeros(xm.size)
    q2 = np.zeros(xm.size)
    point = np.zeros(nodes.size)

    for item in items:
        if item.is_point:
            point[node_index(nodes, item.lcg)[0]] += item.weight
            continue
        node_index(nodes, [item.xa, item.xb])  # elements must not straddle item ends
        wa, wb = trapezoid_ends(item)
        slope = (wb - wa) / item.length
        inside = (xm > item.xa) & (xm < item.xb)
        q1[inside] += wa + slope * (x1[inside] - item.xa)
        q2[inside] += wa + slope * (x2[inside] - item.xa)

    return LoadModel(nodes, q1, q2, point)
