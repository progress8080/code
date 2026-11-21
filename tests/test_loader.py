"""Reading a case: structure positions, block plan, weight items (LCG or moment) and the
weight curve."""

from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, run_case
from shipsupport.loader import structure_positions
from shipsupport.model import InputError
from shipsupport.units import MM, TONNE_FORCE

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "example_case.yaml"


def example_dict():
    return yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- structure

def test_structure_start_end_spacing():
    assert np.allclose(structure_positions({"start": 0.0, "end": 120.0, "spacing": 0.8}),
                       0.8 * np.arange(151))


def test_structure_start_spacing_count_and_start_end_count():
    x = structure_positions({"start": 6.4, "spacing": 2.4, "count": 44})
    assert x.size == 44 and x[-1] == pytest.approx(109.6)
    x = structure_positions({"start": 10.0, "end": 100.0, "count": 31})
    assert x.size == 31 and np.allclose(np.diff(x), 3.0)


def test_structure_explicit_list_wins():
    assert structure_positions({"x": [3, 1, 2], "start": 0, "count": 9, "spacing": 1}).tolist() \
        == [3, 1, 2]


@pytest.mark.parametrize("spec, message", [
    ({"start": 0, "end": 10}, "exactly two"),
    ({"start": 0, "end": 10, "spacing": 1, "count": 11}, "exactly two"),
    ({"start": 0, "spacing": 1, "count": 2.5}, "whole number"),
    ({"start": 0, "spacing": -1, "count": 3}, "spacing"),
    ({"start": 5, "end": 1, "count": 3}, "end must be"),
])
def test_invalid_structure_specs(spec, message):
    with pytest.raises(InputError, match=message):
        structure_positions(spec)


# --------------------------------------------------------------------------- block plan

def test_block_plan_defines_count_spacing_strength_and_gap_of_every_block():
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "count": 4,
                              "spacing": [2.0, 3.0, 4.0], "strength_t": [150, None, 170, 180],
                              "gap_mm": [0.0, 0.5, 0.0, 1.0]}
    case = case_from_dict(data)
    assert [b.x for b in case.blocks] == [10.0, 12.0, 15.0, 19.0]
    assert [b.label for b in case.blocks] == ["B01", "B02", "B03", "B04"]
    # an empty strength is the type capacity (200 t)
    assert [round(b.effective_strength / TONNE_FORCE) for b in case.blocks] == [150, 200, 170, 180]
    assert [b.gap for b in case.blocks] == pytest.approx([0.0, 0.5 * MM, 0.0, 1.0 * MM])
    data["blocks"]["plan"] = {"type": "keel", "start_x": 6.4, "count": 44, "spacing": 2.4,
                              "strength_t": 190, "gap_mm": 0.2}
    case = case_from_dict(data)
    assert len(case.blocks) == 44 and case.blocks[-1].x == pytest.approx(109.6)
    assert all(b.strength == pytest.approx(190 * TONNE_FORCE) for b in case.blocks)
    assert all(b.gap == pytest.approx(0.2 * MM) for b in case.blocks)


def test_block_count_from_the_lists():
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "spacing": 3.0,
                              "strength_t": [100, 120, 140]}
    assert len(case_from_dict(data).blocks) == 3
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "spacing": [3.0, 2.0, 1.0, 4.0]}
    assert len(case_from_dict(data).blocks) == 5


@pytest.mark.parametrize("plan, message", [
    ({"count": 4, "spacing": [2.0, 3.0]}, "3 expected"),
    ({"count": 3, "spacing": [2.0, -1.0]}, "> 0"),
    ({"count": 2.5, "spacing": 2.0}, "whole number"),
    ({"count": 1, "spacing": 2.0}, "whole number"),
    ({"count": 3, "spacing": 2.0, "strength_t": [1, 2]}, "3 expected"),
    ({"count": 3, "spacing": 2.0, "gap_mm": [0, 1, 2, 3]}, "3 expected"),
    ({"count": 3, "spacing": "wide"}, "expected a number"),
    ({"spacing": 2.0}, "count"),
    ({"count": 3, "spacing": 2.0, "type": "pier"}, "unknown block type"),
])
def test_invalid_block_plans_are_rejected(plan, message):
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, **plan}
    with pytest.raises(InputError, match=message):
        case_from_dict(data)


def test_block_groups_are_no_longer_accepted():
    data = example_dict()
    data["blocks"].pop("plan")
    data["blocks"]["groups"] = [{"type": "keel", "start": 6.4, "spacing": 2.4, "count": 44}]
    with pytest.raises(InputError, match="blocks.plan"):
        case_from_dict(data)


def test_changing_count_and_spacing_changes_reactions():
    data = example_dict()
    data["structure"] = None
    data["blocks"]["plan"].update(count=33, spacing=3.2)
    result = run_case(case_from_dict(data))
    fem = result.supports["fem"]
    assert fem.reactions.size == 33
    assert np.allclose(np.diff(fem.block_x), 3.2)
    assert fem.reactions.sum() == pytest.approx(result.load.total, rel=1e-12)


# --------------------------------------------------------------------------- weights

def test_moment_gives_the_lcg():
    data = example_dict()
    data["weights"] = [{"name": "a", "weight_t": 400, "moment_tm": 22000, "xa": 40, "xb": 70},
                       {"name": "b", "weight_t": 100, "moment_tm": 3000}]          # point at 30 m
    items = case_from_dict(data).ship.weights
    assert items[0].lcg == pytest.approx(55.0) and (items[0].xa, items[0].xb) == (40.0, 70.0)
    assert items[1].is_point and items[1].lcg == pytest.approx(30.0)
    same = example_dict()
    same["weights"] = [{"name": "a", "weight_t": 400, "lcg": 55.0, "xa": 40, "xb": 70},
                       {"name": "b", "weight_t": 100, "x": 30.0}]
    assert [(w.weight, w.lcg, w.xa, w.xb) for w in case_from_dict(same).ship.weights] == \
        pytest.approx([(w.weight, w.lcg, w.xa, w.xb) for w in items])


def test_moment_and_lcg_must_agree():
    data = example_dict()
    data["weights"][1] = {"name": "main engine", "weight_t": 380, "lcg": 22.0,
                          "moment_tm": 380 * 22.0, "xa": 16.0, "xb": 28.0}
    assert case_from_dict(data).ship.weights[1].lcg == pytest.approx(22.0)
    data["weights"][1]["moment_tm"] = 380 * 23.0
    with pytest.raises(InputError, match="main engine.*does not agree"):
        case_from_dict(data)


@pytest.mark.parametrize("item, message", [
    ({"weight_t": 10, "xa": 0, "xb": 10}, "give lcg, moment_tm"),
    ({"weight_t": 10, "lcg": 5, "xa": 0}, "both xa and xb"),
    ({"weight_t": 10, "x": 5, "xa": 0, "xb": 10}, "no extent"),
    ({"weight_t": 0, "moment_tm": 5}, "must not be zero"),
])
def test_invalid_weight_items(item, message):
    data = example_dict()
    data["weights"] = [item]
    with pytest.raises(InputError, match=message):
        case_from_dict(data)


def test_weight_curve_is_read_as_items():
    data = example_dict()
    data["weight_curve"] = {"name": "lightship", "points": [
        {"x": 0.0, "q_t_per_m": 10.0}, {"x": 60.0, "q_t_per_m": 20.0},
        {"x": 60.0, "q_t_per_m": 30.0}, {"x": 120.0, "q_t_per_m": 0.0}]}
    ship = case_from_dict(data).ship
    curve_items = [w for w in ship.weights if w.name.startswith("lightship")]
    assert [w.name for w in curve_items] == ["lightship 1/3", "lightship 3/3"]
    assert ship.weight_curve.name == "lightship"
    added = (15.0 * 60 + 15.0 * 60) * TONNE_FORCE
    assert ship.total_weight == pytest.approx(4500 * TONNE_FORCE + added)
    as_list = example_dict()
    as_list["weight_curve"] = data["weight_curve"]["points"]          # a plain list of points
    assert case_from_dict(as_list).ship.total_weight == pytest.approx(ship.total_weight)
    only_curve = example_dict()
    only_curve.pop("weights")
    only_curve["weight_curve"] = data["weight_curve"]
    assert case_from_dict(only_curve).ship.total_weight == pytest.approx(added)


@pytest.mark.parametrize("points, message", [
    ([{"x": 0, "q_t_per_m": 1}, {"x": "aft", "q_t_per_m": 1}], r"weight_curve\.points\[1\]\.x"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 10}], r"weight_curve\.points\[1\].*q_t_per_m"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 10, "q_t_per_m": -1}], r"point 2 \(x = 10 m\)"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 130, "q_t_per_m": 1}], "outside the hull"),
])
def test_invalid_weight_curve_points_are_named(points, message):
    data = example_dict()
    data["weight_curve"] = {"points": points}
    with pytest.raises(InputError, match=message):
        run_case(case_from_dict(data))


def test_no_weights_at_all_is_rejected():
    data = example_dict()
    data.pop("weights")
    with pytest.raises(InputError, match="weight_curve"):
        case_from_dict(data)
