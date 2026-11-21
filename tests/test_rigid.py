import numpy as np
import pytest

from shipsupport.model import WeightItem
from shipsupport.solvers import OPEN, BlockSprings, UnstableSupportError, solve_rigid

from helpers import make_blocks, uniform_setup


def test_two_blocks_statically_determinate():
    # Uniform 1 kN/m over 20 m, blocks at 2 and 15: R2 = W·(10 − 2)/13
    load, _, _, springs = uniform_setup(20.0, 1e3, [2.0, 15.0], k=3e8)
    r = solve_rigid(load, springs)
    W = 2e4
    assert r.reactions[1] == pytest.approx(W * 8 / 13, rel=1e-12)
    assert r.reactions[0] == pytest.approx(W * 5 / 13, rel=1e-12)


def test_equal_blocks_symmetric_load():
    xs = np.linspace(1, 19, 10)
    load, _, _, springs = uniform_setup(20.0, 5e3, xs)
    r = solve_rigid(load, springs)
    assert np.allclose(r.reactions, 1e5 / 10, rtol=1e-12)


def test_eccentric_closed_form():
    xs = np.linspace(2, 18, 9)
    extra = (WeightItem("point", 4e4, 14.0, 14.0, 14.0),)
    load, _, _, springs = uniform_setup(20.0, 1e3, xs, extra_items=extra)
    r = solve_rigid(load, springs)
    W, xbar = load.total, xs.mean()
    e = load.lcg - xbar
    expected = W / xs.size + W * e * (xs - xbar) / np.sum((xs - xbar) ** 2)
    assert np.all(expected > 0)
    assert np.allclose(r.reactions, expected, rtol=1e-10)


def test_lift_off_matches_reduced_block_set():
    # Heavy point load near the forward end lifts the aft block
    xs = [1.0, 5.0, 9.0, 13.0, 17.0, 19.0]
    extra = (WeightItem("p", 4e5, 18.5, 18.5, 18.5),)
    load, _, _, springs = uniform_setup(20.0, 1e3, xs, extra_items=extra)
    r = solve_rigid(load, springs)
    assert r.states[0] == OPEN and r.reactions[0] == 0.0
    assert np.all(r.reactions >= 0)

    reduced = BlockSprings.from_blocks(make_blocks(xs[1:]))
    r2 = solve_rigid(load, reduced)
    assert np.allclose(r.reactions[1:], r2.reactions, rtol=1e-10)
    assert abs(r.force_residual) < 1e-12 and abs(r.moment_residual) < 1e-12


def test_crushing_two_symmetric_blocks():
    k, fy, ratio = 1e9, 3e5, 0.1
    load, _, _, springs = uniform_setup(10.0, 1e5, [2.0, 8.0], k=k, crush_stress=fy,
                                        post_crush_ratio=ratio)
    r = solve_rigid(load, springs)
    dy = fy / k
    expected_delta = dy + (5e5 - fy) / (ratio * k)
    assert np.allclose(r.delta, expected_delta, rtol=1e-12)
    assert np.allclose(r.reactions, 5e5, rtol=1e-12)
    assert r.n_crushed == 2


def test_lcg_outside_support_is_unstable():
    load, _, _, springs = uniform_setup(20.0, 1e3, [1.0, 4.0])
    with pytest.raises(UnstableSupportError):
        solve_rigid(load, springs)
