"""Block plan input and the spacing design for a given number of blocks (spacing.py)."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, load_case, run_case
from shipsupport.model import Block, BlockType, InputError, SpacingDesign
from shipsupport.spacing import design_spacing, equidistribution_step

ROOT = Path(__file__).resolve().parents[1]
TYPE = BlockType("t", 1.0, 1.0, stiffness_override=1e9)


@dataclass
class FakeEvaluation:
    reactions: np.ndarray


def tributary_load(q, a=0.0, b=30.0):
    """Substitute analysis of a hull from a to b: every block carries the load q(x) of its
    tributary cell, and the end blocks also carry the overhangs (the load is conserved)."""
    def evaluate(blocks):
        xs = np.array([b_.x for b_ in blocks])
        edges = np.concatenate([[a], 0.5 * (xs[1:] + xs[:-1]), [b]])
        grid = np.linspace(a, b, 40001)
        mid = 0.5 * (grid[1:] + grid[:-1])
        cell = np.clip(np.searchsorted(edges, mid, side="right") - 1, 0, xs.size - 1)
        return FakeEvaluation(np.bincount(cell, q(mid) * np.diff(grid), minlength=xs.size))
    return evaluate


def blocks_at(xs, strengths):
    return [Block(float(x), TYPE, strength=float(s)) for x, s in zip(xs, strengths, strict=True)]


def free_design(**kw):
    defaults = dict(min_spacing=0.5, max_spacing=10.0, snap_to_structure=False, max_iterations=30)
    defaults.update(kw)
    return SpacingDesign(**defaults)


# --------------------------------------------------------------------------- the method

def test_uniform_load_reaches_the_analytic_optimum():
    """Uniform load q on a hull of length L carried by N equal blocks of strength S: every
    block can carry exactly q L / N, at the cell centres x_i = (i − ½) L / N, so the
    smallest possible largest utilisation is U* = q L / (N S) and all utilisations are equal."""
    n, length, q, s = 11, 30.0, 950.0 / 30.0, 100.0
    xs = np.sort(np.random.default_rng(1).uniform(0.5, 29.5, n))
    result = design_spacing(blocks_at(xs, [s] * n), tributary_load(lambda x: np.full_like(x, q)),
                            free_design(), (0.0, length))
    u_star = q * length / (n * s)
    u = result.final_utilisation
    assert result.fits
    assert np.max(u) == pytest.approx(u_star, rel=0.01)
    assert np.max(u) - np.min(u) < 0.02 * u_star                 # every block equally used
    assert np.max(result.initial_utilisation) > 1.3 * u_star


def test_weak_blocks_get_a_smaller_spacing():
    n = 16
    strengths = [50.0] * 8 + [150.0] * 8                   # weak blocks aft, strong forward
    result = design_spacing(blocks_at(np.linspace(1, 29, n), strengths),
                            tributary_load(lambda x: np.full_like(x, 40.0)),
                            free_design(), (0.0, 30.0))
    gaps = result.spacing
    assert result.fits
    assert np.mean(gaps[1:7]) < 0.6 * np.mean(gaps[8:-1])
    before = result.initial_utilisation
    assert np.max(result.final_utilisation) < 0.75 * np.max(before)
    assert [b.effective_strength for b in result.final_blocks] == strengths   # order kept


def test_equidistribution_step_moves_blocks_towards_a_load_peak():
    xs = np.linspace(0.0, 20.0, 11)
    reactions = np.where(xs < 5.0, 300.0, 50.0)
    y = equidistribution_step(xs, reactions, np.full(11, 100.0))
    assert y[0] == 0.0 and y[-1] == 20.0
    assert np.all(np.diff(y) > 0)
    assert y[2] < xs[2]                                    # denser where the load is high


def test_spacing_limits_extent_and_grid_are_respected():
    q = lambda x: np.where(x < 6.0, 120.0, 15.0)            # noqa: E731
    design = SpacingDesign(min_spacing=0.8, max_spacing=3.0, x_min=-2.0, x_max=34.0,
                           grid=0.4, snap_to_structure=False)
    result = design_spacing(blocks_at(np.linspace(0, 30, 16), [100.0] * 16), tributary_load(q),
                            design, (0.0, 30.0))
    xs = np.array([b.x for b in result.final_blocks])
    gaps = np.diff(xs)
    assert np.all(gaps >= 0.8 - 1e-6) and np.all(gaps <= 3.0 + 1e-6)
    assert xs[0] >= -2.0 - 1e-6 and xs[-1] <= 34.0 + 1e-6
    assert np.allclose(xs / 0.4, np.round(xs / 0.4), atol=1e-6)    # on the 0.4 m grid
    assert [h.phase for h in result.history][0] == "start" and result.history[-1].phase == "grid"


def test_fixed_end_blocks_stay():
    result = design_spacing(blocks_at([1.0, 4.0, 5.0, 9.0, 13.0], [80.0] * 5),
                            tributary_load(lambda x: np.full_like(x, 20.0)),
                            free_design(fix_end_blocks=True), (0.0, 20.0))
    assert result.final_blocks[0].x == 1.0 and result.final_blocks[-1].x == 13.0


def test_reports_when_the_strength_is_not_enough():
    n = 6
    result = design_spacing(blocks_at(np.linspace(0, 20, n), [30.0] * n),
                            tributary_load(lambda x: np.full_like(x, 20.0)),   # 400 > 6 × 30
                            free_design(), (0.0, 20.0))
    assert result.status == "does_not_fit" and not result.fits
    assert result.suggested_count is not None and result.suggested_count > n
    assert "total strength" in result.message


def test_invalid_spacing_problems_are_rejected():
    load = tributary_load(lambda x: np.full_like(x, 1.0))
    with pytest.raises(InputError, match="strength"):
        design_spacing([Block(0.0, TYPE), Block(5.0, TYPE)], load, free_design(), (0.0, 5.0))
    with pytest.raises(InputError, match="minimum spacing"):
        design_spacing(blocks_at(np.linspace(0, 4, 6), [10.0] * 6), load,
                       free_design(min_spacing=1.0), (0.0, 4.0))
    with pytest.raises(InputError):
        SpacingDesign(min_spacing=2.0, max_spacing=1.0)


# --------------------------------------------------------------------------- integration

def test_spacing_example_fits_every_block_within_strength(tmp_path):
    result = run_case(load_case(ROOT / "examples" / "spacing_case.yaml"), ["fem"])
    sp = result.spacing
    assert sp.fits and len(sp.final_blocks) == 44
    assert np.max(sp.initial_utilisation) > 1.0 >= np.max(sp.final_utilisation)
    assert result.strength["fem+spacing"].passed is True
    assert list(result.block_supports) == ["fem"]
    gaps = sp.spacing
    assert gaps.min() >= 0.8 - 1e-6 and gaps.max() <= 4.8 + 1e-6
    structure = np.arange(0.0, 120.01, 0.8)
    for b in sp.final_blocks:                                          # under floors
        assert np.min(np.abs(structure - b.x)) < 1e-6
    assert [b.effective_strength for b in sp.final_blocks] == \
        [b.effective_strength for b in sp.initial_blocks]               # strengths kept in order
    final = result.supports["fem+spacing"]
    assert abs(final.force_residual) < 1e-9 and abs(final.moment_residual) < 1e-9
    from shipsupport.cli import main
    assert main(["run", str(ROOT / "examples" / "spacing_case.yaml"), "-o", str(tmp_path),
                 "--methods", "fem"]) == 0
    for name in ("spacing.csv", "spacing_history.csv", "spacing_design.png"):
        assert (tmp_path / name).exists()
    rows = (tmp_path / "spacing.csv").read_text(encoding="utf-8").splitlines()
    assert len(rows) == 45 and rows[0].startswith("block,x_m,spacing_to_next_m")
    assert "BLOCK SPACING DESIGN" in (tmp_path / "summary.txt").read_text(encoding="utf-8")


def test_initial_gaps_move_with_their_blocks():
    data = yaml.safe_load((ROOT / "examples" / "spacing_case.yaml").read_text(encoding="utf-8"))
    gap = [0.0] * 44
    gap[30] = 0.4
    data["blocks"]["plan"]["gap_mm"] = gap
    sp = run_case(case_from_dict(data), ["fem"]).spacing
    assert sp.final_blocks[30].gap == pytest.approx(0.4e-3)
    assert sp.final_evaluation.springs.gap[30] == pytest.approx(0.4e-3)
    assert np.count_nonzero(sp.final_evaluation.springs.gap) == 1


@pytest.fixture(scope="module")
def welded():
    """The welded hull of welding_case.yaml with the spacing design switched on."""
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    data["spacing_design"]["enabled"] = True
    return run_case(case_from_dict(data))


def test_spacing_design_on_the_welded_hull_recomputes_the_welding_gaps(welded):
    """Welding deformation is a constraint of the spacing design: every candidate layout is
    analysed with the welding gaps at its own block positions."""
    from shipsupport.deformation import gaps_from_lift, keel_lift
    sp = welded.spacing
    x0 = np.array([b.x for b in sp.initial_blocks])
    x1 = np.array([b.x for b in sp.final_blocks])
    assert not np.allclose(x0, x1)
    expected = gaps_from_lift(keel_lift(welded.deformation, x1), x1)
    assert np.allclose(sp.final_evaluation.springs.gap, expected, rtol=0, atol=1e-12)
    assert np.allclose(sp.initial_evaluation.springs.gap, welded.weld_gap, rtol=0, atol=1e-12)
    assert not np.allclose(expected, welded.weld_gap)                  # the gaps moved with x
    # the designed layout analysed from scratch gives the same reactions
    import dataclasses
    fresh = run_case(dataclasses.replace(welded.case, blocks=sp.final_blocks, spacing=None),
                     ["fem"])
    assert np.allclose(fresh.supports["fem"].reactions, welded.supports["fem+spacing"].reactions,
                       rtol=1e-12, atol=1e-6)
    assert np.allclose(fresh.weld_gap, expected, atol=1e-12)


def test_spacing_design_on_the_welded_hull_reports_its_outcome(welded):
    sp = welded.spacing
    u0, u1 = np.max(sp.initial_utilisation), np.max(sp.final_utilisation)
    assert u1 < u0                                                     # better than as built
    assert sp.total_load == pytest.approx(welded.load.total, rel=1e-9)
    if sp.fits:
        assert u1 <= 1.0
    else:
        assert sp.status == "does_not_fit" and sp.suggested_count > len(sp.final_blocks)
        assert any(w.startswith("Spacing design:") for w in welded.warnings)
    from shipsupport.report import format_summary
    text = format_summary(welded)
    assert "recomputed at the designed positions" in text


def scaled_welding(factor):
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    for joint in data["welding_deformation"]["joints"]:
        joint["angle_mrad"] *= factor
    for zone in data["welding_deformation"]["zones"]:
        zone["deck_shrinkage_mm_per_m"] *= factor
        zone["bottom_shrinkage_mm_per_m"] *= factor
    data["spacing_design"]["enabled"] = True
    return run_case(case_from_dict(data), ["fem"])


def test_moderate_welding_deformation_is_absorbed_by_the_spacing_design():
    """Half the example's welding deformation (end lift about 21 mm): the as-built spacing
    overloads blocks, the designed spacing brings every block within strength."""
    result = scaled_welding(0.5)
    sp = result.spacing
    assert 18e-3 < np.max(result.weld_gap) < 23e-3
    assert np.max(sp.initial_utilisation) > 1.4
    assert sp.fits and result.strength["fem+spacing"].passed is True


def test_severe_welding_deformation_is_diagnosed_as_gap_governed(welded):
    """The full welding deformation (end lift 41 mm) cannot be compensated by spacing alone:
    blocks are kept out of contact and the message says that the gaps govern."""
    sp = welded.spacing
    assert not sp.fits
    assert np.count_nonzero(sp.final_evaluation.reactions <= 0) > 0
    assert "gaps rather than the spacing govern" in sp.message
