"""Spacing design: the positions of a given number of blocks with given strengths.

Problem
    N blocks in a fixed order, block i with strength S_i. Find the positions
    x_1 < … < x_N within the allowed extent [x_lo, x_hi], with
    s_min ≤ x_{i+1} − x_i ≤ s_max, that minimise the largest utilisation

        U(x) = max_i R_i(x) / S_i,

    where R_i(x) are the block reactions of the support model. The layout fits when
    U ≤ u_target. The block count and the order of the blocks (and so which block has
    which strength) are not changed; only the spacing is designed.

Method (three phases, every candidate is analysed with the full support model)
    1. Strength-weighted equidistribution, a fixed-point iteration: from the reactions
       of the current layout, a load density along the block line is estimated and the
       blocks are placed so that the load between two neighbours is proportional to
       their mean strength. Fast, global, gives a good start.
    2. Sequential linear programming (SLP) on the minimax problem: the sensitivities
       ∂U_i/∂x_k are computed by finite differences and the linearised problem
           min t  s.t.  U_i + Σ_k J_ik Δ_k ≤ t,  spacing and extent limits,  |Δ_k| ≤ δ
       is solved as a linear programme (HiGHS). A step is accepted only if the analysed
       U decreases; the trust radius δ grows after success and halves after failure.
    3. Rounding to admissible positions (structure positions or a grid) and a 1-opt
       local search that moves the governing block or a neighbour by one grid step.

The analysis is injected as ``evaluate(blocks)`` (dependency injection; in the program it
is analysis.spacing_evaluator, which also recomputes the welding gaps at the candidate
block positions), so the design can be tested with a substitute analysis.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import linprog

from .model import Block, InputError, SpacingDesign

_TOL = 1e-9   # m, positions closer than this are the same position


# --------------------------------------------------------------------------- helpers

def utilisation(reactions: np.ndarray, blocks: Sequence[Block]) -> np.ndarray:
    """R_i / strength_i for every block; NaN where a block has no strength."""
    strength = np.array([np.nan if b.effective_strength is None else b.effective_strength
                         for b in blocks], dtype=float)
    with np.errstate(invalid="ignore"):
        return np.asarray(reactions, dtype=float) / strength


def relabel(blocks: Sequence[Block]) -> list[Block]:
    """Sort by position and label B01, B02, … from aft to forward."""
    ordered = sorted(blocks, key=lambda b: b.x)
    return [dataclasses.replace(b, label=f"B{i + 1:02d}") for i, b in enumerate(ordered)]


def _admissible(y: np.ndarray, s_min: float, s_max: float,
                grid: np.ndarray | None = None) -> np.ndarray | None:
    """Enforce the minimum and maximum spacing (and grid / structure positions) with fixed
    end blocks: one sweep from aft to forward, each block kept in the window that its
    predecessor and the room needed by the remaining blocks leave for it. None if no
    admissible layout is found."""
    n = y.size
    lo, hi = y[0], y[-1]
    if hi - lo < (n - 1) * s_min - _TOL or hi - lo > (n - 1) * s_max + _TOL:
        return None
    y = y.copy()
    for i in range(1, n - 1):
        rest = n - 1 - i                              # gaps still to come
        lower = max(y[i - 1] + s_min, hi - rest * s_max)
        upper = min(y[i - 1] + s_max, hi - rest * s_min)
        if upper < lower - _TOL:
            return None
        if grid is None:
            y[i] = min(max(y[i], lower), upper)
            continue
        ok = grid[(grid >= lower - _TOL) & (grid <= upper + _TOL)]
        if ok.size == 0:
            return None
        y[i] = ok[np.argmin(np.abs(ok - y[i]))]
    gap = hi - y[-2] if n > 2 else hi - lo
    if gap < s_min - _TOL or gap > s_max + _TOL:
        return None
    return y


@dataclass(frozen=True)
class SpacingStep:
    phase: str                 # start | equidistribution | slp | grid
    iteration: int
    max_utilisation: float
    min_utilisation: float
    evaluations: int           # analyses performed so far


@dataclass(eq=False)
class SpacingResult:
    status: str                # fits | does_not_fit
    message: str
    initial_blocks: list[Block]
    final_blocks: list[Block]
    initial_evaluation: object
    final_evaluation: object
    history: list[SpacingStep] = field(default_factory=list)
    evaluations: int = 0
    target: float = 1.0
    suggested_count: int | None = None     # blocks probably needed when it does not fit
    total_strength: float = 0.0            # N, Σ S_i
    total_load: float = 0.0                # N, Σ R_i

    @property
    def fits(self) -> bool:
        return self.status == "fits"

    @property
    def initial_utilisation(self) -> np.ndarray:
        return utilisation(self.initial_evaluation.reactions, self.initial_blocks)

    @property
    def final_utilisation(self) -> np.ndarray:
        return utilisation(self.final_evaluation.reactions, self.final_blocks)

    @property
    def spacing(self) -> np.ndarray:
        """Designed distance from each block to the next (N − 1 values)."""
        return np.diff([b.x for b in self.final_blocks])

    @property
    def initial_spacing(self) -> np.ndarray:
        return np.diff([b.x for b in self.initial_blocks])


Evaluate = Callable[[Sequence[Block]], object]   # returns an object with .reactions


class _Counter:
    """Wraps evaluate(): counts the analyses and returns (evaluation, U_i)."""

    def __init__(self, evaluate: Evaluate, blocks: Sequence[Block], strength: np.ndarray):
        self.evaluate, self.blocks, self.strength = evaluate, list(blocks), strength
        self.count = 0

    def blocks_at(self, x: np.ndarray) -> list[Block]:
        return [dataclasses.replace(b, x=round(float(v), 6)) for b, v in zip(self.blocks, x, strict=True)]

    def __call__(self, x: np.ndarray):
        self.count += 1
        evaluation = self.evaluate(self.blocks_at(x))
        return evaluation, np.maximum(np.asarray(evaluation.reactions, dtype=float), 0.0) / self.strength


def _peak(u: np.ndarray) -> float:
    return float(np.max(u))


def _score(u: np.ndarray) -> tuple[float, float]:
    """Largest utilisation first, then a smooth measure of the high ones (tie-break)."""
    return round(float(np.max(u)), 9), float(np.sum(u ** 8))


# --------------------------------------------------------------------------- phase 1

def equidistribution_step(x: np.ndarray, reactions: np.ndarray, strength: np.ndarray) -> np.ndarray:
    """New positions with the ends fixed: the load between two neighbouring blocks is made
    proportional to their mean strength.

    The reactions spread uniformly over the tributary cells give a load density; its
    cumulative Q(x) is inverted at the cumulative interval weights (S_i + S_{i+1}) / 2.
    """
    edges = np.concatenate([[x[0]], 0.5 * (x[1:] + x[:-1]), [x[-1]]])
    load = np.maximum(reactions, 0.0)
    load = np.maximum(load, 0.05 * max(float(np.mean(load)), 1e-12))     # lifted-off blocks
    Q = np.concatenate([[0.0], np.cumsum(load)])
    w = 0.5 * (strength[:-1] + strength[1:])
    F = np.concatenate([[0.0], np.cumsum(w)]) / np.sum(w)
    y = np.interp(F * Q[-1], Q, edges)
    y[0], y[-1] = x[0], x[-1]
    return y


# --------------------------------------------------------------------------- phase 2

def _sensitivities(run: _Counter, x: np.ndarray, u: np.ndarray, free: np.ndarray,
                   lo: float, hi: float, s_min: float, s_max: float, h: float) -> np.ndarray:
    """J[i, k] = ∂U_i / ∂x_k by one-sided finite differences (admissible direction)."""
    n = x.size
    J = np.zeros((n, n))
    for k in np.flatnonzero(free):
        for step in (h, -h):
            y = x.copy()
            y[k] += step
            gaps = np.diff(y)
            if y[0] < lo - _TOL or y[-1] > hi + _TOL or np.any(gaps < s_min - _TOL) or \
                    np.any(gaps > s_max + _TOL):
                continue
            _, uk = run(y)
            J[:, k] = (uk - u) / step
            break
    return J


def _slp_step(x: np.ndarray, u: np.ndarray, J: np.ndarray, free: np.ndarray, lo: float,
              hi: float, s_min: float, s_max: float, radius: float, smoothing: float):
    """Solve the linearised minimax problem

        min  t + μ Σ_i |x'_{i+2} − 2 x'_{i+1} + x'_i|          (x' = x + Δ)
        s.t. U_i + J_i Δ ≤ t,   s_min ≤ x'_{i+1} − x'_i ≤ s_max,   x'_1 ≥ lo,  x'_N ≤ hi,
             |Δ_k| ≤ radius

    as a linear programme. The small second term (μ = ``smoothing``) prefers a spacing
    that changes gradually among layouts with practically the same largest utilisation.
    Returns (Δx, predicted max U) or None.
    """
    n = x.size
    idx = np.flatnonzero(free)
    m = idx.size
    if m == 0:
        return None
    k2 = max(n - 2, 0)
    # variables: Δ+ (m), Δ− (m), a (n − 2, |second difference|), t;   Δ = Δ+ − Δ−
    nv = 2 * m + k2 + 1
    c = np.zeros(nv)
    c[:2 * m] = 1e-6                                       # no needless moves
    c[2 * m:2 * m + k2] = smoothing
    c[-1] = 1.0
    col = {k: j for j, k in enumerate(idx)}

    def delta(k, sign=1.0):
        r = np.zeros(nv)
        if k in col:
            r[col[k]] += sign
            r[m + col[k]] -= sign
        return r

    rows, rhs = [], []
    Jf = J[:, idx]
    for i in range(n):                                     # U_i + J_i Δ − t ≤ 0
        r = np.zeros(nv)
        r[:m], r[m:2 * m], r[-1] = Jf[i], -Jf[i], -1.0
        rows.append(r)
        rhs.append(-u[i])
    for i in range(n - 1):
        gap = x[i + 1] - x[i]
        d = delta(i + 1) - delta(i)                         # change of gap i
        if d.any():
            rows.append(-d)
            rhs.append(gap - s_min)
            if np.isfinite(s_max):
                rows.append(d)
                rhs.append(s_max - gap)
    for i in range(k2):                                    # ± second difference ≤ a_i
        d2 = x[i + 2] - 2 * x[i + 1] + x[i]
        r = delta(i + 2) - 2 * delta(i + 1) + delta(i)
        a = np.zeros(nv)
        a[2 * m + i] = 1.0
        rows.append(r - a)
        rhs.append(-d2)
        rows.append(-r - a)
        rhs.append(d2)
    ends = (0, n - 1)                                      # only the end blocks touch the extent
    up = [max(min(radius, hi - x[k]) if k in ends else radius, 0.0) for k in idx]
    down = [max(min(radius, x[k] - lo) if k in ends else radius, 0.0) for k in idx]
    bounds = [(0.0, b) for b in up] + [(0.0, b) for b in down] + [(0.0, None)] * k2 + [(None, None)]
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), bounds=bounds, method="highs")
    if res.status != 0:
        return None
    step = np.zeros(n)
    step[idx] = res.x[:m] - res.x[m:2 * m]
    return step, float(res.x[-1])


# --------------------------------------------------------------------------- phase 3

def _grid_neighbour(grid: np.ndarray, x: float, direction: int) -> float | None:
    """The next admissible grid position from x in the given direction."""
    if direction > 0:
        ahead = grid[grid > x + 1e-6]
        return float(ahead[0]) if ahead.size else None
    behind = grid[grid < x - 1e-6]
    return float(behind[-1]) if behind.size else None


def _local_moves(x: np.ndarray, u: np.ndarray, grid: np.ndarray, free: np.ndarray, lo: float,
                 hi: float, s_min: float, s_max: float, hot_tolerance: float = 0.02) -> list:
    """Neighbouring grid layouts for the local search: every governing block (within 2 % of
    the largest utilisation, at most three) and its neighbours are moved by one or two grid
    positions, alone or together with the next block; only admissible layouts are returned."""
    n = x.size
    peak = float(np.max(u))
    hot = [int(i) for i in np.argsort(-u)[:3] if u[i] >= (1.0 - hot_tolerance) * peak]
    movable = sorted({k for i in hot for k in (i - 1, i, i + 1) if 0 <= k < n and free[k]})
    trials, seen = [], set()
    for k in movable:
        for direction in (-1, 1):
            for steps in (1, 2):
                for group in ((k,), (k, k + 1)):
                    if group[-1] >= n or not all(free[g] for g in group):
                        continue
                    z = x.copy()
                    for g in group:
                        v = z[g]
                        for _ in range(steps):
                            v = _grid_neighbour(grid, v, direction)
                            if v is None:
                                break
                        if v is None:
                            break
                        z[g] = v
                    else:
                        gaps = np.diff(z)
                        key = tuple(np.round(z, 6))
                        if key in seen or z[0] < lo - _TOL or z[-1] > hi + _TOL or \
                                np.any(gaps < s_min - _TOL) or np.any(gaps > s_max + _TOL):
                            continue
                        seen.add(key)
                        trials.append(z)
    return trials


def _round_to_grid(x: np.ndarray, grid: np.ndarray, fixed_ends: bool, lo: float, hi: float,
                   s_min: float, s_max: float) -> np.ndarray | None:
    """Snap the free end blocks to the nearest grid position within the extent, then the
    interior blocks in one admissible sweep (see _admissible)."""
    y = x.copy()
    if not fixed_ends:
        inside = grid[(grid >= lo - _TOL) & (grid <= hi + _TOL)]
        if inside.size == 0:
            return None
        for k in (0, -1):
            y[k] = inside[np.argmin(np.abs(inside - y[k]))]
    return _admissible(y, s_min, s_max, grid)


# --------------------------------------------------------------------------- driver

def design_spacing(blocks: Sequence[Block], evaluate: Evaluate, design: SpacingDesign,
                   x_bounds: tuple[float, float], structure: np.ndarray | None = None,
                   equidistribution_iterations: int = 15, fd_step: float = 0.05,
                   local_search_evaluations: int | None = None) -> SpacingResult:
    """Design the positions of the given blocks (count and order kept).

    ``x_bounds`` is the allowed extent of the block line (used unless the design sets
    x_min / x_max); ``structure`` gives admissible final positions if
    ``design.snap_to_structure``; otherwise ``design.grid`` (if set) is used.
    """
    current = relabel(blocks)
    n = len(current)
    strength = np.array([np.nan if b.effective_strength is None else b.effective_strength
                         for b in current], dtype=float)
    if np.any(~np.isfinite(strength)) or np.any(strength <= 0):
        missing = ", ".join(b.label for b, s in zip(current, strength, strict=True)
                            if not np.isfinite(s) or s <= 0)
        raise InputError(f"Spacing design needs a strength for every block; missing: {missing}")
    target = design.target_utilisation
    x0 = np.array([b.x for b in current])
    lo = design.x_min if design.x_min is not None else x_bounds[0]
    hi = design.x_max if design.x_max is not None else x_bounds[1]
    if design.fix_end_blocks:
        lo, hi = x0[0], x0[-1]
    s_min = design.min_spacing
    # default maximum spacing: twice the mean spacing of the given line (an unlimited gap
    # would let the minimax optimum leave long unsupported spans amidships)
    s_max = design.max_spacing if design.max_spacing is not None else \
        max(2.0 * (x0[-1] - x0[0]) / max(n - 1, 1), s_min)
    if (n - 1) * s_min > hi - lo + 1e-6:
        raise InputError(f"Spacing design: {n} blocks at a minimum spacing of {s_min:g} m need "
                         f"{(n - 1) * s_min:.2f} m, but the block line is only {hi - lo:.2f} m "
                         f"long ({lo:g}–{hi:g} m)")
    if np.isfinite(s_max) and (n - 1) * s_max < hi - lo - 1e-6:
        if design.fix_end_blocks:
            raise InputError(f"Spacing design: {n} blocks at a maximum spacing of {s_max:g} m "
                             f"cannot span the fixed block line {lo:g}–{hi:g} m")
        hi = lo + (n - 1) * s_max      # free ends: shorten the line to what the blocks can span
    run = _Counter(evaluate, current, strength)
    history: list[SpacingStep] = []

    def record(phase, it, u):
        history.append(SpacingStep(phase, it, _peak(u), float(np.min(u)), run.count))

    first_eval, u_first = run(x0)
    record("start", 0, u_first)

    # starting layout: ends at the extent limits, interior mapped proportionally
    x = x0.copy()
    if not design.fix_end_blocks:
        x = lo + (x0 - x0[0]) * (hi - lo) / max(x0[-1] - x0[0], _TOL)
        x[0], x[-1] = lo, hi
    x = _admissible(x, s_min, s_max)
    if x is None:
        raise InputError("Spacing design: no admissible start layout within the spacing limits")
    evaluation, u = run(x)
    best = (x, evaluation, u)

    # phase 1 — strength-weighted equidistribution (fixed point, under-relaxed)
    for it in range(1, equidistribution_iterations + 1):
        y = equidistribution_step(x, evaluation.reactions, strength)
        y = _admissible(x + 0.5 * (y - x), s_min, s_max)
        if y is None or np.max(np.abs(y - x)) < 1e-3:
            break
        x = y
        evaluation, u = run(x)
        record("equidistribution", it, u)
        if _score(u) < _score(best[2]):
            best = (x, evaluation, u)
    x, evaluation, u = best

    # phase 2 — sequential linear programming on the minimax problem
    free = np.ones(n, dtype=bool)
    if design.fix_end_blocks:
        free[[0, -1]] = False
    mean_gap = (hi - lo) / max(n - 1, 1)
    radius, radius_max, radius_min = 0.25 * mean_gap, mean_gap, 0.005
    smoothing = 0.002 / mean_gap        # utilisation traded for 1 m of spacing variation
    J = None
    for it in range(1, design.max_iterations + 1):
        if J is None:
            J = _sensitivities(run, x, u, free, lo, hi, s_min, s_max, fd_step)
        step = _slp_step(x, u, J, free, lo, hi, s_min, s_max, radius, smoothing)
        if step is None:
            break
        delta, predicted = step
        if _peak(u) - predicted < 1e-5:
            break                                          # linearised optimum reached
        y = x + delta
        y[0], y[-1] = max(y[0], lo), min(y[-1], hi)
        if np.any(np.diff(y) < s_min - 1e-6):
            y = _admissible(y, s_min, s_max)
            if y is None:
                radius *= 0.5
                continue
        ev, uy = run(y)
        if _score(uy) < _score(u):
            x, evaluation, u = y, ev, uy
            record("slp", it, u)
            radius = min(2.0 * radius, radius_max)
            J = None
        else:
            radius *= 0.5
            if radius < radius_min:
                break

    # phase 3 — admissible positions (structure / grid) and 1-opt local search
    grid = None
    if design.snap_to_structure and structure is not None and np.size(structure):
        grid = np.sort(np.asarray(structure, dtype=float))
    elif design.grid:
        grid = np.arange(np.ceil((lo - _TOL) / design.grid), np.floor((hi + _TOL) / design.grid) + 1) \
            * design.grid
    if grid is not None:
        y = _round_to_grid(x, grid, design.fix_end_blocks, lo, hi, s_min, s_max)
        if y is not None:
            x, (evaluation, u) = y, run(y)
            budget = local_search_evaluations if local_search_evaluations is not None else 25 * n
            start = run.count
            while run.count - start < budget:
                trials = _local_moves(x, u, grid, free, lo, hi, s_min, s_max)
                best_trial = None
                for z in trials:
                    if run.count - start >= budget:
                        break
                    ev, uz = run(z)
                    if best_trial is None or _score(uz) < _score(best_trial[2]):
                        best_trial = (z, ev, uz)
                if best_trial is None or not _score(best_trial[2]) < _score(u):
                    break                                  # local optimum on the grid
                x, evaluation, u = best_trial
            record("grid", 0, u)

    final = run.blocks_at(x)
    peak = _peak(u)
    total_strength = float(np.sum(strength))
    total_load = float(np.sum(np.maximum(evaluation.reactions, 0.0)))
    if peak <= target * (1 + 1e-9):
        status = "fits"
        message = (f"all {n} blocks within strength: max reaction/strength {peak:.3f} "
                   f"(target {target:g}) after {run.count} analyses")
        suggested = None
    else:
        status = "does_not_fit"
        suggested = max(n + 1, int(np.ceil(n * peak / target)))
        if total_strength * target < total_load:
            reason = (f"the total strength Σ S = {total_strength / 9806.65:.1f} t × {target:g} "
                      f"is less than the load {total_load / 9806.65:.1f} t")
        else:
            reason = f"the best spacing gives max reaction/strength {peak:.3f}"
        worst = int(np.argmax(u))
        if worst in (0, n - 1) and not design.fix_end_blocks:
            side = "aft" if worst == 0 else "forward"
            limit = "x_min" if worst == 0 else "x_max"
            reason += (f", governed by the {side} end block at the limit of the block line — "
                       f"allow a longer line ({limit}) to shorten the overhang")
        open_blocks = int(np.count_nonzero(np.asarray(evaluation.reactions) <= 0.0))
        if open_blocks:
            message = (f"{n} blocks with these strengths cannot be spaced within the target: "
                       f"{reason}; {open_blocks} block(s) carry no load because of gaps under "
                       f"the keel (e.g. welding deformation), so the gaps rather than the "
                       f"spacing govern and more blocks may not help — adjust the block "
                       f"heights to the keel line, or use stronger blocks where the keel "
                       f"bears (about {suggested} blocks would be needed if the load spread)")
        else:
            message = (f"{n} blocks with these strengths cannot be spaced within the target: "
                       f"{reason}; or try about {suggested} blocks or stronger blocks")
    return SpacingResult(status, message, current, final, first_eval, evaluation, history,
                         run.count, target, suggested, total_strength, total_load)
