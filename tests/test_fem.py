import numpy as np
import pytest

from shipsupport.model import Section, WeightItem
from shipsupport.solvers import OPEN, BlockSprings, solve_fem, solve_rigid
from shipsupport.weights import build_load_model, build_nodes

from helpers import interp_v, make_blocks, uniform_setup

RIGID_K = 1e15


def test_simply_supported_midspan_deflection_euler_bernoulli():
    L, q, EI = 10.0, 1e3, 1e8
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=RIGID_K, EI=EI, max_len=0.5)
    r = solve_fem(load, sections, springs)
    sag = interp_v(r, L / 2) - 0.5 * (interp_v(r, 0.0) + interp_v(r, L))
    assert sag == pytest.approx(5 * q * L**4 / (384 * EI), rel=1e-9)


def test_simply_supported_midspan_deflection_timoshenko():
    L, q, EI, kGA = 10.0, 1e3, 1e8, 1e6
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=RIGID_K, EI=EI, kGA=kGA,
                                               max_len=0.5)
    r = solve_fem(load, sections, springs, shear_deformation=True)
    sag = interp_v(r, L / 2) - 0.5 * (interp_v(r, 0.0) + interp_v(r, L))
    expected = 5 * q * L**4 / (384 * EI) + q * L**2 / (8 * kGA)
    assert sag == pytest.approx(expected, rel=1e-9)


def test_shear_flag_off_ignores_kGA():
    load, sections, _, springs = uniform_setup(10.0, 1e3, [0.0, 10.0], k=RIGID_K, EI=1e8,
                                               kGA=1e6)
    r = solve_fem(load, sections, springs, shear_deformation=False)
    sag = interp_v(r, 5.0) - 0.5 * (interp_v(r, 0.0) + interp_v(r, 10.0))
    assert sag == pytest.approx(5 * 1e3 * 1e4 / (384 * 1e8), rel=1e-9)


def test_two_span_continuous_beam_reactions():
    L, q = 20.0, 1e3
    load, sections, _, springs = uniform_setup(L, q, [0.0, 10.0, 20.0], k=RIGID_K, EI=1e8)
    r = solve_fem(load, sections, springs)
    W = q * L
    assert np.allclose(r.reactions / W, [3 / 16, 10 / 16, 3 / 16], rtol=1e-6)


def test_stiff_hull_converges_to_rigid_method():
    # The FEM–rigid difference is proportional to 1/EI (until round-off dominates)
    xs = np.arange(3.0, 98.0, 2.5)
    extra = (WeightItem("engine", 2e6, 25.0, 20.0, 30.0), WeightItem("p", 5e5, 80.0, 80.0, 80.0))
    errors = []
    for EI in (1e16, 1e17, 1e18):
        load, sections, _, springs = uniform_setup(100.0, 4e4, list(xs), k=2e9, EI=EI,
                                                   extra_items=extra, max_len=0.5)
        fem = solve_fem(load, sections, springs).reactions
        rigid = solve_rigid(load, springs).reactions
        errors.append(np.max(np.abs(fem - rigid)) / np.max(rigid))
    assert errors[0] / errors[1] == pytest.approx(10.0, rel=0.1)
    assert errors[1] / errors[2] == pytest.approx(10.0, rel=0.1)
    assert errors[2] < 5e-5


def test_flexible_hull_concentrates_load_under_weight():
    xs = np.arange(3.0, 98.0, 2.5)
    extra = (WeightItem("engine", 2e6, 25.0, 20.0, 30.0),)
    load, sections, _, springs = uniform_setup(100.0, 4e4, list(xs), k=2e9, EI=4e12,
                                               extra_items=extra)
    fem = solve_fem(load, sections, springs)
    rigid = solve_rigid(load, springs)
    i = int(np.argmin(np.abs(xs - 25.0)))
    assert fem.reactions[i] > 1.05 * rigid.reactions[i]
    assert abs(fem.force_residual) < 1e-9 and abs(fem.moment_residual) < 1e-9


def test_large_gap_equals_removed_block():
    xs = [2.0, 6.0, 10.0, 14.0, 18.0]
    load, sections, _, springs = uniform_setup(20.0, 1e4, xs, k=1e9, EI=1e9,
                                               gaps=[0, 0, 0.05, 0, 0])
    r = solve_fem(load, sections, springs)
    assert r.states[2] == OPEN and r.reactions[2] == 0.0

    reduced = BlockSprings.from_blocks(make_blocks([2.0, 6.0, 14.0, 18.0], k=1e9))
    r2 = solve_fem(load, sections, reduced)
    assert np.allclose(np.delete(r.reactions, 2), r2.reactions, rtol=1e-9)


def test_small_gap_partially_closes():
    xs = [2.0, 6.0, 10.0, 14.0, 18.0]
    load, sections, _, springs = uniform_setup(20.0, 1e5, xs, k=1e9, EI=1e9,
                                               gaps=[0, 0, 2e-4, 0, 0])
    r = solve_fem(load, sections, springs)
    assert 0.0 < r.reactions[2] < r.reactions[1]
    assert r.delta[2] == pytest.approx(r.reactions[2] / 1e9, rel=1e-9)


def test_mesh_independent_reactions():
    # Hermite elements with consistent loads are nodally exact for piecewise-linear
    # loads, so reactions must not depend on the mesh density.
    xs = list(np.arange(2.0, 59.0, 3.0))
    extra = (WeightItem("w", 1e6, 15.0, 10.0, 20.0),)

    def solve(max_len):
        load, sections, _, springs = uniform_setup(60.0, 2e4, xs, k=1e9, EI=2e11,
                                                   extra_items=extra, max_len=max_len)
        return solve_fem(load, sections, springs).reactions

    reference = solve(0.1)
    for max_len in (3.0, 1.5, 0.5):
        assert np.max(np.abs(solve(max_len) - reference)) < 1e-7 * np.max(reference)


def test_varying_sections_must_cover_hull():
    items = [WeightItem("u", 1e5, 5.0, 0.0, 10.0)]
    nodes = build_nodes([0.0, 10.0, 2.0, 8.0], 0.0, 10.0, 0.5)
    load = build_load_model(items, nodes)
    springs = BlockSprings.from_blocks(make_blocks([2.0, 8.0]))
    with pytest.raises(ValueError):
        solve_fem(load, [Section(0.0, 6.0, 1e9)], springs)
