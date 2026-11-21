from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, load_case, run_case
from shipsupport.cli import main
from shipsupport.report import format_summary
from shipsupport.units import TONNE_FORCE

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "example_case.yaml"


@pytest.fixture(scope="module")
def result():
    return run_case(load_case(EXAMPLE))


def test_example_equilibrium(result):
    assert result.load.total == pytest.approx(4500 * TONNE_FORCE, rel=1e-12)
    for r in result.supports.values():
        assert abs(r.force_residual) < 1e-9
        assert abs(r.moment_residual) < 1e-9
        assert np.all(r.reactions >= 0)
    for s in result.strength.values():
        assert abs(s.closure_shear) < 1e-6 and abs(s.closure_moment) < 1e-6


def test_summary_and_cli(result, tmp_path):
    text = format_summary(result)
    assert "SUPPORT SOLUTIONS" in text and "LANDING" not in text and "SHIM" not in text
    assert main(["run", str(EXAMPLE), "-o", str(tmp_path)]) == 0
    for name in ("summary.txt", "blocks.csv", "curves_fem.csv", "curves_rigid.csv",
                 "results.png", "block_forces.png"):
        assert (tmp_path / name).exists()
    for name in ("landing.csv", "landing.png", "layout.csv", "welding_and_shims.png",
                 "spacing.csv"):
        assert not (tmp_path / name).exists()
    header = (tmp_path / "blocks.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "shim" not in header and "strength_t" in header


def test_cli_survey_case_outputs(tmp_path):
    assert main(["run", str(EXAMPLE.with_name("welding_survey_case.yaml")), "-o", str(tmp_path),
                 "--methods", "fem"]) == 0
    assert (tmp_path / "welding_data.png").exists()
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "identified angles" in summary


def test_cli_welding_case_outputs(tmp_path):
    assert main(["run", str(EXAMPLE.with_name("welding_case.yaml")), "-o", str(tmp_path)]) == 0
    assert (tmp_path / "welding_data.png").exists()
    assert not (tmp_path / "spacing.csv").exists()               # spacing_design.enabled: false
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "WELDING DEFORMATION" in summary and "SHIM" not in summary


def test_weight_curve_case_runs(tmp_path):
    case = load_case(EXAMPLE.with_name("weight_curve_case.yaml"))
    ship = case.ship
    assert ship.total_weight == pytest.approx(4643.0 * TONNE_FORCE, rel=1e-12)
    assert ship.lcg == pytest.approx(56.639, abs=5e-4)
    curve = ship.weight_curve
    assert curve.x.size == 18
    assert sum(w.weight for w in curve.items()) == pytest.approx(4343.0 * TONNE_FORCE, rel=1e-12)
    assert main(["run", str(EXAMPLE.with_name("weight_curve_case.yaml")), "-o", str(tmp_path),
                 "--methods", "fem", "--no-plots"]) == 0
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "of which curve" in summary and "BLOCK SPACING DESIGN" in summary


def items_as_curve(items):
    """The piecewise-linear curve (with jumps) that is the sum of the trapezoidal items."""
    from shipsupport.weights import trapezoid_ends
    parts = [(w.xa, w.xb, *trapezoid_ends(w)) for w in items]
    breaks = sorted({p[0] for p in parts} | {p[1] for p in parts})

    def value(x, side):
        q = 0.0
        for xa, xb, wa, wb in parts:
            inside = xa < x <= xb if side < 0 else xa <= x < xb
            if inside:
                q += wa + (wb - wa) * (x - xa) / (xb - xa)
        return q

    points = []
    for x in breaks:
        left, right = value(x, -1), value(x, +1)
        points.append({"x": x, "q_t_per_m": left / TONNE_FORCE})
        if right != left:
            points.append({"x": x, "q_t_per_m": right / TONNE_FORCE})
    return points


def test_weight_curve_equivalent_to_the_example_items_gives_identical_reactions(result):
    """The example's distributed items as one weight curve (sum of trapezoids = a piecewise-
    linear curve with jumps) give the same load and the same block reactions."""
    data = yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))
    case = load_case(EXAMPLE)
    distributed = [w for w in case.ship.weights if not w.is_point]
    data["weight_curve"] = {"name": "items as curve", "points": items_as_curve(distributed)}
    data["weights"] = [w for w in data["weights"] if "x" in w]           # the point load stays
    curved = run_case(case_from_dict(data))
    assert curved.load.total == pytest.approx(result.load.total, rel=1e-12)
    assert curved.load.lcg == pytest.approx(result.load.lcg, rel=1e-12)
    for m in result.supports:
        assert np.allclose(curved.supports[m].reactions, result.supports[m].reactions,
                           rtol=1e-9, atol=1e-9 * result.load.total)


def test_blocks_built_in_code_get_labels():
    # regression: blocks without a label were shown with empty names
    import dataclasses
    from shipsupport.model import Block
    case = load_case(EXAMPLE)
    case.blocks = [Block(b.x, b.type, b.gap) for b in reversed(case.blocks)]
    case.blocks[0] = dataclasses.replace(case.blocks[0], gap=0.05)      # forward end block lifts off
    result = run_case(case)
    labels = result.strength["fem"].blocks.label
    assert labels[0] == "B01" and labels[-1] == f"B{len(labels):02d}"
    assert any("carry no load" in w and f"B{len(labels):02d}" in w for w in result.warnings)


def test_negative_weight_tolerance_rejected():
    from shipsupport.model import InputError, Settings
    with pytest.raises(InputError):
        Settings(weight_tolerance=-0.1)
