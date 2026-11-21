from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .contact import CRUSHED, OPEN


@dataclass(eq=False)
class SupportResult:
    method: str
    block_x: np.ndarray
    block_nodes: np.ndarray
    reactions: np.ndarray     # N, upward on the hull
    delta: np.ndarray         # m, block compression (negative = gap still open)
    states: np.ndarray
    iterations: int
    nodes: np.ndarray
    v: np.ndarray             # m, hull displacement at nodes, downward positive
    theta: np.ndarray         # rad, dv/dx at nodes
    total_load: float         # N
    load_moment: float        # N·m about x = 0

    @property
    def force_residual(self) -> float:
        """(ΣR − W) / W."""
        return (self.reactions.sum() - self.total_load) / self.total_load

    @property
    def moment_residual(self) -> float:
        """(ΣR·x − W·LCG) / (W·span)."""
        span = max(np.ptp(self.nodes), 1e-12)
        return (np.dot(self.reactions, self.block_x) - self.load_moment) / (self.total_load * span)

    @property
    def n_open(self) -> int:
        return int(np.count_nonzero(self.states == OPEN))

    @property
    def n_crushed(self) -> int:
        return int(np.count_nonzero(self.states == CRUSHED))
