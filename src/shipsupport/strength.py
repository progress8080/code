"""Hull girder shear force, bending moment, stresses and block checks.

Sign conventions:
    V(x) = Σ upward forces − Σ downward loads to the left (aft) of x
    M(x) = ∫ V dx, so sagging is positive and hogging negative
    stresses are positive in tension (sagging puts the deck in compression)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .model import Block, Section, Settings, Ship
from .solvers.contact import STATE_NAMES
from .solvers.result import SupportResult
from .weights import LoadModel


@dataclass(frozen=True)
class Check:
    name: str
    unit: str                  # "N", "N·m" or "Pa"
    demand: float              # governing value (signed)
    capacity: float | None     # limit at the governing location
    utilisation: float         # |demand| / capacity, nan if no limit is defined
    x: float                   # m, governing location

    @property
    def passed(self) -> bool | None:
        return None if np.isnan(self.utilisation) else bool(self.utilisation <= 1.0)


@dataclass(eq=False)
class BlockTable:
    x: np.ndarray
    label: list[str]
    reaction: np.ndarray           # N
    pressure: np.ndarray           # Pa
    pressure_util: np.ndarray      # nan if no allowable
    capacity_util: np.ndarray      # nan if no capacity
    state: list[str]
    structure_offset: np.ndarray   # m to nearest structure, nan if none defined
    off_structure: np.ndarray      # bool


@dataclass(eq=False)
class StrengthResult:
    method: str
    x: np.ndarray                  # sample points; nodes appear twice to show jumps in V
    V: np.ndarray                  # N
    M: np.ndarray                  # N·m
    sigma_deck: np.ndarray         # Pa
    sigma_keel: np.ndarray         # Pa
    tau: np.ndarray                # Pa
    bm_util: np.ndarray
    sf_util: np.ndarray
    closure_shear: float           # V at the forward end / max |V|
    closure_moment: float          # M at the forward end / max |M|
    blocks: BlockTable
    checks: list[Check]

    @property
    def passed(self) -> bool | None:
        results = [c.passed for c in self.checks if c.passed is not None]
        return None if not results else all(results)


def shear_moment_curves(load: LoadModel, block_nodes: np.ndarray, reactions: np.ndarray,
                        samples_per_element: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact V(x) and M(x) for linear element loads, nodal point loads and block reactions."""
    nodes = load.nodes
    nodal = -load.point.astype(float)
    np.add.at(nodal, block_nodes, reactions)
    frac = np.arange(1, samples_per_element) / samples_per_element

    xs, Vs, Ms = [], [], []
    V = M = 0.0
    for j in range(nodes.size):
        xs.append([nodes[j], nodes[j]])
        Vs.append([V, V + nodal[j]])
        Ms.append([M, M])
        V += nodal[j]
        if j == nodes.size - 1:
            break
        L = nodes[j + 1] - nodes[j]
        q1 = load.q1[j]
        dq = load.q2[j] - q1
        s = L * frac
        xs.append(nodes[j] + s)
        Vs.append(V - (q1 * s + dq * s**2 / (2 * L)))
        Ms.append(M + V * s - (q1 * s**2 / 2 + dq * s**3 / (6 * L)))
        M += V * L - (q1 * L**2 / 2 + dq * L**2 / 6)
        V -= q1 * L + dq * L / 2
    return np.concatenate(xs), np.concatenate(Vs), np.concatenate(Ms)


def section_property(sections: Sequence[Section], x: np.ndarray, attr: str) -> np.ndarray:
    out = np.full(x.shape, np.nan)
    for s in sections:
        value = getattr(s, attr)
        if value is None:
            continue
        m = np.isnan(out) & (x >= s.xa - 1e-9) & (x <= s.xb + 1e-9)
        out[m] = value
    return out


def governing(name: str, unit: str, x, demand, capacity) -> Check | None:
    """The location with the highest utilisation, or the largest |demand| if no limits exist."""
    x = np.asarray(x, dtype=float)
    demand = np.asarray(demand, dtype=float)
    capacity = np.broadcast_to(np.asarray(capacity, dtype=float), demand.shape)
    finite = np.isfinite(demand)
    if not finite.any():
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        util = np.abs(demand) / capacity
    valid = finite & np.isfinite(util)
    if valid.any():
        i = np.flatnonzero(valid)[np.argmax(util[valid])]
        return Check(name, unit, float(demand[i]), float(capacity[i]), float(util[i]), float(x[i]))
    i = np.flatnonzero(finite)[np.argmax(np.abs(demand[finite]))]
    return Check(name, unit, float(demand[i]), None, float("nan"), float(x[i]))


def block_table(support: SupportResult, blocks: Sequence[Block], ship: Ship,
                settings: Settings) -> BlockTable:
    R = support.reactions
    area = np.array([b.type.area for b in blocks])
    allowable = np.array([np.nan if b.type.allowable_pressure is None else b.type.allowable_pressure
                          for b in blocks])
    capacity = np.array([np.nan if b.effective_strength is None else b.effective_strength for b in blocks])
    pressure = R / area
    if ship.structure_positions is not None:
        offset = np.min(np.abs(support.block_x[:, None] - ship.structure_positions[None, :]), axis=1)
        off = offset > settings.structure_tolerance + 1e-9
    else:
        offset = np.full(R.size, np.nan)
        off = np.zeros(R.size, dtype=bool)
    return BlockTable(
        x=support.block_x,
        label=[b.label for b in blocks],
        reaction=R,
        pressure=pressure,
        pressure_util=pressure / allowable,
        capacity_util=R / capacity,
        state=[STATE_NAMES[int(s)] for s in support.states],
        structure_offset=offset,
        off_structure=off,
    )


def compute_strength(load: LoadModel, support: SupportResult, blocks: Sequence[Block], ship: Ship,
                     settings: Settings, samples_per_element: int = 8) -> StrengthResult:
    x, V, M = shear_moment_curves(load, support.block_nodes, support.reactions, samples_per_element)

    z_deck = section_property(ship.sections, x, "z_deck")
    z_keel = section_property(ship.sections, x, "z_keel")
    shear_area = section_property(ship.sections, x, "shear_area")
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_deck = -M / z_deck
        sigma_keel = M / z_keel
        tau = V / shear_area

    if ship.permissible_bm is not None:
        sag, hog = ship.permissible_bm.at(x)
        bm_cap = np.where(M >= 0.0, sag, hog)
    else:
        bm_cap = np.full(x.shape, np.nan)
    sf_cap = ship.permissible_sf.at(x) if ship.permissible_sf is not None else np.full(x.shape, np.nan)
    bm_util = np.abs(M) / bm_cap
    sf_util = np.abs(V) / sf_cap

    table = block_table(support, blocks, ship, settings)
    allow_b = np.nan if settings.allowable_bending_stress is None else settings.allowable_bending_stress
    allow_s = np.nan if settings.allowable_shear_stress is None else settings.allowable_shear_stress
    allowable_p = np.array([np.nan if b.type.allowable_pressure is None else b.type.allowable_pressure
                            for b in blocks])
    capacity = np.array([np.nan if b.effective_strength is None else b.effective_strength for b in blocks])

    checks = [
        governing("Hull bending moment", "N·m", x, M, bm_cap),
        governing("Hull shear force", "N", x, V, sf_cap),
        governing("Deck bending stress", "Pa", x, sigma_deck, allow_b),
        governing("Keel bending stress", "Pa", x, sigma_keel, allow_b),
        governing("Hull shear stress", "Pa", x, tau, allow_s),
        governing("Block bearing pressure", "Pa", table.x, table.pressure, allowable_p),
        governing("Block load vs capacity", "N", table.x, table.reaction, capacity),
    ]

    v_max = max(np.max(np.abs(V)), 1e-12)
    m_max = max(np.max(np.abs(M)), 1e-12)
    return StrengthResult(
        method=support.method,
        x=x, V=V, M=M,
        sigma_deck=sigma_deck, sigma_keel=sigma_keel, tau=tau,
        bm_util=bm_util, sf_util=sf_util,
        closure_shear=float(V[-1] / v_max),
        closure_moment=float(M[-1] / m_max),
        blocks=table,
        checks=[c for c in checks if c is not None],
    )
