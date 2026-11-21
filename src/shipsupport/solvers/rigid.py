"""Method A: rigid hull on elastic blocks.

The hull displacement is linear, v(x) = a + b·(x − x0), so equilibrium alone
gives a 2×2 system for each set of block states:

    [ Σkt      Σkt·xc  ] [a]   [ W  + Σ(kt·g − c)    ]
    [ Σkt·xc   Σkt·xc² ] [b] = [ M0 + Σ(kt·g − c)·xc ]

with xc = x − x0 (x0 = block centroid, for conditioning) and M0 the load
moment about x0.
"""

from __future__ import annotations

import numpy as np

from ..weights import LoadModel, node_index
from .contact import BlockSprings, UnstableSupportError, solve_contact
from .result import SupportResult


def solve_rigid(load: LoadModel, springs: BlockSprings, max_iterations: int = 200) -> SupportResult:
    x0 = float(springs.x.mean())
    xc = springs.x - x0
    W = load.total
    M0 = load.moment - W * x0

    def linear_solve(kt, c, gap):
        r = kt * gap - c
        A = np.array([[kt.sum(), np.dot(kt, xc)],
                      [np.dot(kt, xc), np.dot(kt, xc * xc)]])
        try:
            a, b = np.linalg.solve(A, [W + r.sum(), M0 + np.dot(r, xc)])
        except np.linalg.LinAlgError as exc:
            raise UnstableSupportError("Rigid-hull system is singular") from exc
        return a + b * xc - gap, (a, b)

    delta, states, iterations, (a, b) = solve_contact(springs, linear_solve, max_iterations)
    kt, c = springs.linearize(states)
    nodes = load.nodes
    return SupportResult(
        method="rigid",
        block_x=springs.x,
        block_nodes=node_index(nodes, springs.x),
        reactions=kt * delta + c,
        delta=delta,
        states=states,
        iterations=iterations,
        nodes=nodes,
        v=a + b * (nodes - x0),
        theta=np.full(nodes.size, b),
        total_load=W,
        load_moment=load.moment,
    )
