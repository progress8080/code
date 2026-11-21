"""Compression-only (optionally crushing) block springs and the contact-state
iteration shared by the rigid and FEM solvers.

Each block force is piecewise linear in its compression δ = v − gap:

    open     δ ≤ 0          R = 0
    elastic  0 < δ ≤ δy     R = k·δ
    crushed  δ > δy         R = k₂·δ + (k − k₂)·δy

For a fixed set of states the problem is linear, so we iterate on the states
until they are consistent with the solution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from ..model import Block

OPEN, ELASTIC, CRUSHED = 0, 1, 2
STATE_NAMES = {OPEN: "open", ELASTIC: "elastic", CRUSHED: "crushed"}


class ContactError(RuntimeError):
    """The contact iteration failed."""


class UnstableSupportError(ContactError):
    """Too few blocks in contact to support the ship."""


@dataclass(frozen=True, eq=False)
class BlockSprings:
    x: np.ndarray
    k: np.ndarray
    gap: np.ndarray
    crush_deflection: np.ndarray   # δy, inf where the block does not crush
    k_crushed: np.ndarray

    @classmethod
    def from_blocks(cls, blocks: Sequence[Block]) -> "BlockSprings":
        x = np.array([b.x for b in blocks], dtype=float)
        k = np.array([b.type.stiffness for b in blocks], dtype=float)
        gap = np.array([b.gap for b in blocks], dtype=float)
        fy = np.array([np.inf if b.type.crush_force is None else b.type.crush_force for b in blocks])
        k2 = np.array([b.type.post_crush_ratio for b in blocks]) * k
        return cls(x, k, gap, fy / k, k2)

    def __len__(self) -> int:
        return self.x.size

    def states(self, delta: np.ndarray) -> np.ndarray:
        s = np.full(delta.shape, ELASTIC)
        s[delta > self.crush_deflection] = CRUSHED
        s[delta <= 0.0] = OPEN
        return s

    def linearize(self, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(kt, c) such that R = kt·δ + c for the given states."""
        dy = np.where(np.isfinite(self.crush_deflection), self.crush_deflection, 0.0)
        kt = np.select([states == ELASTIC, states == CRUSHED], [self.k, self.k_crushed], 0.0)
        c = np.where(states == CRUSHED, (self.k - self.k_crushed) * dy, 0.0)
        return kt, c

    def force(self, delta: np.ndarray) -> np.ndarray:
        kt, c = self.linearize(self.states(delta))
        return kt * delta + c

    def violation(self, delta: np.ndarray, states: np.ndarray) -> np.ndarray:
        """How far each block's δ lies outside the range of its assumed state."""
        dy = self.crush_deflection
        v = np.zeros(delta.shape)
        v = np.where(states == OPEN, np.maximum(delta, 0.0), v)
        v = np.where(states == ELASTIC, np.maximum(-delta, 0.0) + np.maximum(delta - dy, 0.0), v)
        v = np.where(states == CRUSHED, np.maximum(dy - delta, 0.0), v)
        return v


# linear_solve(kt, c, gap) -> (δ at each block, extra data)
LinearSolve = Callable[[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, object]]


def is_supported(x: np.ndarray, kt: np.ndarray) -> bool:
    stiff = kt > 0.0
    return np.count_nonzero(stiff) >= 2 and np.ptp(x[stiff]) >= 1e-9


def check_support(x: np.ndarray, kt: np.ndarray) -> None:
    if not is_supported(x, kt):
        stiff = kt > 0.0
        raise UnstableSupportError(
            f"Only {np.count_nonzero(stiff)} block position(s) carry load with stiffness; at least "
            f"two separated blocks are needed. Check block positions, gaps and the LCG.")


def _iterate(springs: BlockSprings, linear_solve: LinearSolve, gap: np.ndarray,
             states: np.ndarray, max_iterations: int):
    """State iteration for one gap vector, starting from the given states.

    All states are updated at once; if a set of states repeats (cycling), or a
    full update would leave the hull without support, only the worst violator
    is changed per iteration.
    """
    seen: set[bytes] = set()
    one_at_a_time = False
    for iteration in range(1, max_iterations + 1):
        kt, c = springs.linearize(states)
        check_support(springs.x, kt)
        delta, extra = linear_solve(kt, c, gap)
        new_states = springs.states(delta)
        changed = new_states != states
        if not changed.any():
            return delta, states, iteration, extra

        seen.add(states.tobytes())
        unsupported = not is_supported(springs.x, springs.linearize(new_states)[0])
        if one_at_a_time or unsupported or new_states.tobytes() in seen:
            one_at_a_time = True
            worst = np.where(changed, springs.violation(delta, states), -np.inf)
            i = int(np.argmax(worst))
            new_states = states.copy()
            new_states[i] = springs.states(delta)[i]
        states = new_states
    raise ContactError(f"Block contact iteration did not converge in {max_iterations} iterations")


def solve_contact(springs: BlockSprings, linear_solve: LinearSolve,
                  max_iterations: int = 200) -> tuple[np.ndarray, np.ndarray, int, object]:
    """Iterate block states to consistency; returns (δ, states, iterations, extra).

    The state iteration is tried directly. If it fails — large, uneven gaps such as
    a welding end lift can make it cycle or lose all support on the way — a gap
    continuation is used: the problem is first solved with all gaps equal to the
    smallest gap (a rigid shift, so every block starts in contact), and the gaps
    are then increased in steps to their true values, each step starting from the
    previous contact states. Failed steps are halved.
    """
    initial = np.full(len(springs), ELASTIC)
    try:
        return _iterate(springs, linear_solve, springs.gap, initial, max_iterations)
    except ContactError as direct_error:
        first_error = direct_error

    g0 = float(springs.gap.min())
    target = springs.gap
    try:
        delta, states, total, extra = _iterate(springs, linear_solve, np.full_like(target, g0),
                                               initial, max_iterations)
    except ContactError:
        raise first_error
    lam, step = 0.0, 0.25
    while lam < 1.0:
        nxt = min(1.0, lam + step)
        gap = g0 + nxt * (target - g0)
        try:
            delta, states_next, iterations, extra = _iterate(springs, linear_solve, gap, states,
                                                             max_iterations)
        except ContactError:
            step /= 2.0
            if step < 1.0 / 4096:
                raise first_error
            continue
        total += iterations
        lam, states = nxt, states_next
        step = min(2.0 * step, 0.5)
    return delta, states, total, extra
