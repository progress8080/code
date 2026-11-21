import numpy as np
import pytest

from shipsupport.model import BendingEnvelope, Section, Settings, ShearEnvelope, Ship, WeightItem
from shipsupport.solvers import solve_fem
from shipsupport.strength import compute_strength, shear_moment_curves

from helpers import uniform_setup


def test_simply_supported_shear_and_moment():
    L, q = 10.0, 2e3
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=1e15, EI=1e9)
    r = solve_fem(load, sections, springs)
    x, V, M = shear_moment_curves(load, r.block_nodes, r.reactions)
    assert M.max() == pytest.approx(q * L**2 / 8, rel=1e-9)        # sagging positive
    assert x[np.argmax(M)] == pytest.approx(L / 2)
    assert V.max() == pytest.approx(q * L / 2, rel=1e-9)
    assert abs(V[-1]) < 1e-6 and abs(M[-1]) < 1e-6


def test_overhangs_give_hogging_over_supports():
    L, q, a = 10.0, 1e3, 2.0
    load, sections, _, springs = uniform_setup(L, q, [a, L - a], k=1e15, EI=1e9)
    r = solve_fem(load, sections, springs)
    x, V, M = shear_moment_curves(load, r.block_nodes, r.reactions)
    at_support = M[np.argmin(np.abs(x - a))]
    assert at_support == pytest.approx(-q * a**2 / 2, rel=1e-9)
    mid = M[np.argmin(np.abs(x - L / 2))]
    assert mid == pytest.approx(q * L * (L / 2 - a) / 2 - q * (L / 2)**2 / 2, rel=1e-9)


def test_checks_and_utilisation():
    L, q = 10.0, 2e3
    load, _, blocks, springs = uniform_setup(L, q, [0.0, L], k=1e15, EI=1e9)
    sections = [Section(0.0, L, 1e9, None, z_deck=1e-4, z_keel=2e-4, shear_area=1e-3)]
    ship = Ship("t", L, [WeightItem("u", q * L, L / 2, 0.0, L)], sections,
                permissible_bm=BendingEnvelope([0.0, L], [5e4, 5e4], [5e4, 5e4]),
                permissible_sf=ShearEnvelope([0.0, L], [2e4, 2e4]))
    settings = Settings(allowable_bending_stress=200e6, allowable_shear_stress=100e6)
    r = solve_fem(load, sections, springs)
    s = compute_strength(load, r, blocks, ship, settings)
    checks = {c.name: c for c in s.checks}

    bm = checks["Hull bending moment"]
    assert bm.utilisation == pytest.approx(q * L**2 / 8 / 5e4, rel=1e-9)
    deck = checks["Deck bending stress"]
    assert deck.demand == pytest.approx(-(q * L**2 / 8) / 1e-4, rel=1e-9)   # compression
    assert deck.passed is False
    assert abs(checks["Hull shear stress"].demand) == pytest.approx(1e4 / 1e-3, rel=1e-9)
    assert s.passed is False
