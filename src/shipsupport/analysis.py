"""Run a complete docking analysis for a Case."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from .deformation import Identification, effective_model, gaps_from_lift, keel_lift
from .model import Block, Case, InputError, Ship, WeldingDeformation
from .solvers import OPEN, BlockSprings, SupportResult, solve_fem, solve_rigid
from .spacing import SpacingResult, design_spacing
from .strength import StrengthResult, compute_strength
from .units import MM, to_tonnes
from .weights import LoadModel, build_load_model, build_nodes

METHODS = ("rigid", "fem")


@dataclass(eq=False)
class AnalysisResult:
    case: Case
    load: LoadModel
    springs: BlockSprings
    supports: dict[str, SupportResult]
    strength: dict[str, StrengthResult]
    warnings: list[str] = field(default_factory=list)
    weld_gap: np.ndarray | None = None              # m per block, from welding deformation
    deformation: WeldingDeformation | None = None   # model actually used (after identification)
    identification: Identification | None = None
    spacing: SpacingResult | None = None

    @property
    def design_method(self) -> str:
        return "fem" if "fem" in self.supports else next(iter(self.supports))

    @property
    def block_supports(self) -> dict[str, SupportResult]:
        """Support solutions on the given block plan (the designed spacing is reported
        separately)."""
        return {m: r for m, r in self.supports.items() if not m.endswith("+spacing")}


@dataclass(frozen=True, eq=False)
class SpacingEvaluation:
    """One analysed candidate block layout (what the spacing design's evaluate() returns)."""

    load: LoadModel
    support: SupportResult
    springs: BlockSprings          # with the welding gaps at the candidate block positions

    @property
    def reactions(self) -> np.ndarray:
        return self.support.reactions


def _load_model(ship: Ship, blocks: list[Block], max_element_length: float) -> LoadModel:
    """Mesh with a node at every weight-item end, section end and block, and the load on it."""
    keys = ([w.xa for w in ship.weights] + [w.xb for w in ship.weights]
            + [s.xa for s in ship.sections] + [s.xb for s in ship.sections]
            + [b.x for b in blocks])
    return build_load_model(ship.weights,
                            build_nodes(keys, ship.x_aft, ship.x_fwd, max_element_length))


def _springs(blocks: list[Block], deformation: WeldingDeformation | None) -> BlockSprings:
    """Block springs; the welded hull's keel lift at the block positions adds to the gaps."""
    springs = BlockSprings.from_blocks(blocks)
    if deformation is not None:
        springs = dataclasses.replace(
            springs, gap=springs.gap + gaps_from_lift(keel_lift(deformation, springs.x), springs.x))
    return springs


def spacing_evaluator(case: Case, deformation: WeldingDeformation | None, method: str):
    """evaluate(blocks) for the spacing design: mesh, welding gaps at the candidate block
    positions and support solution of a candidate layout, exactly as run_case does for the
    given layout."""
    ship, settings = case.ship, case.settings

    def evaluate(blocks: list[Block]) -> SpacingEvaluation:
        load = _load_model(ship, blocks, settings.max_element_length)
        springs = _springs(blocks, deformation)
        if method == "rigid":
            support = solve_rigid(load, springs, settings.max_iterations)
        else:
            support = solve_fem(load, ship.sections, springs, settings.shear_deformation,
                                settings.max_iterations)
        return SpacingEvaluation(load, support, springs)

    return evaluate


def run_case(case: Case, methods: Iterable[str] | None = None) -> AnalysisResult:
    methods = tuple(methods or case.settings.methods)
    unknown = set(methods) - set(METHODS)
    if not methods or unknown:
        raise InputError(f"Unknown method(s) {sorted(unknown)}; choose from {METHODS}")

    ship, settings = case.ship, case.settings
    # sorted by x; blocks built in code may have no label, so give them B01, B02, … too
    blocks = [b if b.label else dataclasses.replace(b, label=f"B{i + 1:02d}")
              for i, b in enumerate(sorted(case.blocks, key=lambda b: b.x))]
    load = _load_model(ship, blocks, settings.max_element_length)
    warnings = _check_weights(case, load)

    x = np.array([b.x for b in blocks])
    if not x.min() <= load.lcg <= x.max():
        raise InputError(f"LCG {load.lcg:.3f} m lies outside the block line "
                         f"[{x.min():.3f}, {x.max():.3f}] m; the ship would tip")

    # Welding deformation: the welded hull's keel is not straight, which adds to the gaps.
    weld_gap = deformation = identification = None
    if case.deformation is not None:
        given = case.deformation
        if given.measured is not None and not given.measured.identify_joints and \
                (given.joints or given.zones):
            warnings.append("Welding deformation: the surveyed keel line is used as the total "
                            "shape; joints and zones given as well are ignored (set "
                            "identify_joint_angles to fit the joint angles to the survey)")
        deformation, identification = effective_model(given)
        weld_gap = gaps_from_lift(keel_lift(deformation, x), x)
    springs = _springs(blocks, deformation)

    def solve(method: str) -> SupportResult:
        if method == "rigid":
            return solve_rigid(load, springs, settings.max_iterations)
        return solve_fem(load, ship.sections, springs, settings.shear_deformation,
                         settings.max_iterations)

    supports: dict[str, SupportResult] = {m: solve(m) for m in methods}
    strength = {m: compute_strength(load, r, blocks, ship, settings) for m, r in supports.items()}

    for method, r in supports.items():
        if r.n_open:
            lifted = ", ".join(blocks[i].label for i in np.flatnonzero(r.states == OPEN))
            cause = " (welding end lift / open gap)" if weld_gap is not None else \
                " (lift-off or open gap)"
            warnings.append(f"{method}: {r.n_open} block(s) carry no load{cause}: {lifted}")
        if r.n_crushed:
            warnings.append(f"{method}: {r.n_crushed} block(s) loaded beyond crush force")
    if weld_gap is not None:
        warnings.append(f"Welding deformation: keel lift relative to the block line up to "
                        f"{weld_gap.max() / MM:.1f} mm (aft end {weld_gap[0] / MM:.1f} mm, "
                        f"forward end {weld_gap[-1] / MM:.1f} mm)")
    design = "fem" if "fem" in supports else methods[0]
    table = strength[design].blocks
    for i in [i for i, off in enumerate(table.off_structure) if off]:
        warnings.append(f"Block {table.label[i]} at x = {table.x[i]:.3f} m is "
                        f"{table.structure_offset[i]:.3f} m from the nearest structure")
    for method, s in strength.items():
        if abs(s.closure_shear) > 1e-6 or abs(s.closure_moment) > 1e-6:
            warnings.append(f"{method}: shear/moment closure error {s.closure_shear:.2e} / "
                            f"{s.closure_moment:.2e}")

    # Spacing design: positions of the given blocks (count and strengths kept); every
    # candidate is analysed with the welding gaps at its own block positions.
    spacing = None
    if case.spacing is not None:
        sd = case.spacing
        x_lo = max(ship.x_aft, sd.x_min if sd.x_min is not None else blocks[0].x)
        x_hi = min(ship.x_fwd, sd.x_max if sd.x_max is not None else blocks[-1].x)
        spacing = design_spacing(blocks, spacing_evaluator(case, deformation, sd.method),
                                 dataclasses.replace(sd, x_min=x_lo, x_max=x_hi), (x_lo, x_hi),
                                 ship.structure_positions)
        key = f"{sd.method}+spacing"
        supports[key] = dataclasses.replace(spacing.final_evaluation.support, method=key)
        strength[key] = compute_strength(spacing.final_evaluation.load, supports[key],
                                         spacing.final_blocks, ship, settings)
        if not spacing.fits:
            warnings.append(f"Spacing design: {spacing.message}")

    return AnalysisResult(case, load, springs, supports, strength, warnings,
                          weld_gap, deformation, identification, spacing)


def _check_weights(case: Case, load: LoadModel) -> list[str]:
    ship, settings = case.ship, case.settings
    warnings = []
    if abs(load.total - ship.total_weight) > 1e-9 * abs(ship.total_weight) or \
            abs(load.lcg - ship.lcg) > 1e-9 * ship.lpp:
        raise RuntimeError("Internal error: distributed load does not reproduce the weight items")
    if ship.declared_weight is not None:
        diff = (ship.total_weight - ship.declared_weight) / ship.declared_weight
        if abs(diff) > settings.weight_tolerance:
            warnings.append(f"Sum of weight items {to_tonnes(ship.total_weight):.1f} t differs from "
                            f"declared {to_tonnes(ship.declared_weight):.1f} t ({100 * diff:+.2f} %)")
    if ship.declared_lcg is not None:
        diff = ship.lcg - ship.declared_lcg
        if abs(diff) > settings.weight_tolerance * ship.lpp:
            warnings.append(f"LCG of weight items {ship.lcg:.3f} m differs from declared "
                            f"{ship.declared_lcg:.3f} m ({diff:+.3f} m)")
    return warnings
