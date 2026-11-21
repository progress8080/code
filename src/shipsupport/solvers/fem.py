"""Method B: elastic hull girder (Euler–Bernoulli / Timoshenko beam elements)
on compression-only block springs.

DOFs per node: v (downward) and θ = dv/dx. The beam is free–free; it is held
only by the springs, so at least two separated blocks must be in contact.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..model import InputError, Section
from ..weights import LoadModel, node_index
from .contact import BlockSprings, UnstableSupportError, solve_contact
from .result import SupportResult


def element_properties(sections: Sequence[Section], nodes: np.ndarray,
                       shear_deformation: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """EI and kGA for each element (kGA = inf means no shear deformation)."""
    xm = 0.5 * (nodes[:-1] + nodes[1:])
    EI = np.full(xm.size, np.nan)
    kGA = np.full(xm.size, np.inf)
    for s in sections:
        m = np.isnan(EI) & (xm >= s.xa) & (xm <= s.xb)
        EI[m] = s.EI
        if shear_deformation and s.kGA is not None:
            kGA[m] = s.kGA
    missing = np.isnan(EI)
    if missing.any():
        raise InputError(f"No hull section covers x = {xm[missing][0]:.3f} m")
    return EI, kGA


def element_stiffness(L: np.ndarray, EI: np.ndarray, kGA: np.ndarray) -> np.ndarray:
    """Timoshenko beam element stiffness, shape (n, 4, 4); Φ = 0 gives Euler–Bernoulli."""
    L = np.asarray(L, dtype=float)
    phi = 12.0 * EI / (kGA * L**2)
    c = EI / (L**3 * (1.0 + phi))
    one = np.ones_like(L)
    k = np.array([
        [12 * one,  6 * L,                  -12 * one,  6 * L],
        [6 * L,     (4 + phi) * L**2,       -6 * L,     (2 - phi) * L**2],
        [-12 * one, -6 * L,                 12 * one,   -6 * L],
        [6 * L,     (2 - phi) * L**2,       -6 * L,     (4 + phi) * L**2],
    ])
    return np.moveaxis(k, -1, 0) * c[:, None, None]


def element_loads(L: np.ndarray, q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Consistent nodal loads for a linear load q1 → q2, shape (n, 4)."""
    return np.column_stack([
        L * (7 * q1 + 3 * q2) / 20,
        L**2 * (3 * q1 + 2 * q2) / 60,
        L * (3 * q1 + 7 * q2) / 20,
        -L**2 * (2 * q1 + 3 * q2) / 60,
    ])


def assemble(load: LoadModel, EI: np.ndarray, kGA: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    nodes = load.nodes
    ndof = 2 * nodes.size
    L = load.lengths
    e = np.arange(L.size)
    dofs = np.column_stack([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])

    ke = element_stiffness(L, EI, kGA)
    rows = np.repeat(dofs, 4, axis=1).ravel()
    cols = np.tile(dofs, (1, 4)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()

    F = np.zeros(ndof)
    np.add.at(F, dofs, element_loads(L, load.q1, load.q2))
    F[0::2] += load.point
    return K, F


def solve_fem(load: LoadModel, sections: Sequence[Section], springs: BlockSprings,
              shear_deformation: bool = True, max_iterations: int = 200) -> SupportResult:
    nodes = load.nodes
    ndof = 2 * nodes.size
    EI, kGA = element_properties(sections, nodes, shear_deformation)
    K0, F0 = assemble(load, EI, kGA)
    block_nodes = node_index(nodes, springs.x)
    vdofs = 2 * block_nodes

    def linear_solve(kt, c, gap):
        kdiag = np.zeros(ndof)
        np.add.at(kdiag, vdofs, kt)
        F = F0.copy()
        np.add.at(F, vdofs, kt * gap - c)
        u = spla.spsolve((K0 + sp.diags(kdiag)).tocsc(), F)
        if not np.all(np.isfinite(u)):
            raise UnstableSupportError("FEM system is singular")
        return u[vdofs] - gap, u

    delta, states, iterations, u = solve_contact(springs, linear_solve, max_iterations)
    kt, c = springs.linearize(states)
    return SupportResult(
        method="fem",
        block_x=springs.x,
        block_nodes=block_nodes,
        reactions=kt * delta + c,
        delta=delta,
        states=states,
        iterations=iterations,
        nodes=nodes,
        v=u[0::2],
        theta=u[1::2],
        total_load=load.total,
        load_moment=load.moment,
    )
