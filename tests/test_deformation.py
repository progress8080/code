from pathlib import Path

import dataclasses

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, run_case
from shipsupport.deformation import (gaps_from_lift, identify_joint_angles, interpolate_profile,
                                     keel_lift, kink_shape, zone_shape)
from shipsupport.model import (InputError, MeasuredProfile, ShrinkageZone, WeldingDeformation,
                               WeldJoint)
from shipsupport.solvers import ELASTIC, OPEN, solve_fem, solve_rigid

from helpers import uniform_setup

ROOT = Path(__file__).resolve().parents[1]


def test_kink_changes_slope_by_angle():
    x = np.array([9.0, 10.0, 11.0])
    z = 2e-3 * kink_shape(x, 10.0)
    slope_left, slope_right = z[1] - z[0], z[2] - z[1]
    assert slope_right - slope_left == pytest.approx(2e-3)


def test_zone_has_unit_curvature_inside_and_is_straight_outside():
    x = np.linspace(-5, 25, 3001)
    z = zone_shape(x, 0.0, 20.0)
    curvature = np.gradient(np.gradient(z, x), x)
    inside = (x > 1) & (x < 19)
    outside = (x < -1) | (x > 21)
    assert np.allclose(curvature[inside], 1.0, atol=1e-6)
    assert np.allclose(curvature[outside], 0.0, atol=1e-6)
    assert zone_shape(np.array([20.0]), 0, 20)[0] == pytest.approx(200.0)


def test_gaps_ignore_straight_line_part():
    x = np.linspace(5, 95, 20)
    lift = 1e-3 * (x - 50) ** 2 / 100
    tilted = lift + 0.4 + 3e-3 * x
    assert np.allclose(gaps_from_lift(lift, x), gaps_from_lift(tilted, x), atol=1e-15)
    assert gaps_from_lift(lift, x).min() == 0.0


def test_reactions_unchanged_by_rigid_tilt_of_keel_line():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    gap = 1e-3 * ((x - 50) / 50) ** 2 * 10          # concave-up keel, 10 mm at the ends
    a = solve_fem(load, sections, dataclasses.replace(springs, gap=gap))
    b = solve_fem(load, sections, dataclasses.replace(springs, gap=gap + 0.01 + 2e-4 * x))
    assert np.allclose(a.reactions, b.reactions, rtol=1e-8, atol=1e-6 * load.total)


def test_large_end_lift_unloads_first_and_last_blocks():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    defo = WeldingDeformation(joints=tuple(WeldJoint(xj, 1.5e-3) for xj in (20, 40, 60, 80)))
    gap = gaps_from_lift(keel_lift(defo, x), x)
    fem = solve_fem(load, sections, dataclasses.replace(springs, gap=gap))
    assert fem.states[0] == OPEN and fem.states[-1] == OPEN
    assert fem.reactions[0] == 0.0 and fem.reactions[-1] == 0.0
    # gap terms k·g are ~40·W here, so round-off in the residual is larger than usual
    assert abs(fem.force_residual) < 1e-7 and abs(fem.moment_residual) < 1e-7
    # complementarity: open blocks have no compression, loaded blocks have positive compression
    assert np.all(fem.delta[fem.states == OPEN] <= 0)
    assert np.all(fem.delta[fem.states == ELASTIC] > 0)


def test_continuation_solves_rigid_case_with_large_gaps():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    gap = 0.04 * ((x - 50) / 47) ** 2                # 40 mm end lift: few blocks touch a rigid hull
    rigid = solve_rigid(load, dataclasses.replace(springs, gap=gap))
    assert rigid.n_open > 20
    assert abs(rigid.force_residual) < 1e-12 and abs(rigid.moment_residual) < 1e-12
    assert np.all(rigid.delta[rigid.states == OPEN] <= 0)
    assert np.all(rigid.delta[rigid.states == ELASTIC] > 0)


def test_identification_recovers_joint_angles_despite_datum():
    joints = (20.0, 45.0, 70.0, 95.0)
    angles = np.array([0.8e-3, 0.4e-3, 0.5e-3, 0.7e-3])
    xs = np.linspace(0, 120, 25)
    z = sum(a * kink_shape(xs, xj) for a, xj in zip(angles, joints)) + 0.123 - 4e-4 * xs
    ident = identify_joint_angles(MeasuredProfile(xs, z), joints)
    assert np.allclose([j.angle for j in ident.joints], angles, rtol=1e-9, atol=1e-12)
    assert ident.rms_residual < 1e-12


def test_identification_with_noise_and_regularisation():
    rng = np.random.default_rng(3)
    joints = (20.0, 45.0, 70.0, 95.0)
    angles = np.array([0.8e-3, 0.4e-3, 0.5e-3, 0.7e-3])
    xs = np.linspace(0, 120, 49)
    z = sum(a * kink_shape(xs, xj) for a, xj in zip(angles, joints))
    noisy = z + rng.normal(0, 0.5e-3, xs.size)                 # 0.5 mm survey noise
    ident = identify_joint_angles(MeasuredProfile(xs, noisy), joints, regularization=1e-6)
    assert np.allclose([j.angle for j in ident.joints], angles, atol=0.15e-3)
    assert ident.rms_residual == pytest.approx(0.5e-3, rel=0.35)


@pytest.mark.parametrize("method", ["pchip", "linear", "spline"])
def test_profile_interpolation_passes_through_survey_points(method):
    profile = MeasuredProfile([0, 10, 25, 40, 60], [5e-3, 2e-3, 0.0, 1e-3, 6e-3],
                              interpolation=method)
    assert np.allclose(interpolate_profile(profile, profile.x), profile.lift, atol=1e-15)


def test_loader_converts_shrinkage_to_curvature():
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    case = case_from_dict(data)
    zone = case.deformation.zones[0]
    assert zone.curvature == pytest.approx((0.24e-3 - 0.20e-3) / 10.0)
    assert case.deformation.joints[0].angle == pytest.approx(0.9e-3)
    assert case.spacing is None                                   # spacing_design.enabled: false


def test_invalid_deformation_input():
    with pytest.raises(InputError):
        WeldingDeformation()
    with pytest.raises(InputError):
        ShrinkageZone(10.0, 5.0, 1e-5)
    with pytest.raises(InputError):
        WeldingDeformation(measured=MeasuredProfile([0, 1], [0, 0], identify_joints=True))


def test_survey_csv_and_identification_with_known_zones(tmp_path):
    from shipsupport import load_case
    from shipsupport.deformation import effective_model
    from shipsupport.loader import read_survey_csv

    xs, zs = read_survey_csv(ROOT / "examples" / "welding_survey.csv")
    assert len(xs) == 25 and xs[0] == 0.0 and xs[-1] == 120.0
    case = load_case(ROOT / "examples" / "welding_survey_case.yaml")
    assert case.deformation.measured.identify_joints
    used, ident = effective_model(case.deformation)
    true = [0.9e-3, 0.5e-3, 0.4e-3, 0.5e-3, 0.6e-3]          # used to generate the sample survey
    assert np.allclose([j.angle for j in used.joints], true, atol=0.1e-3)
    assert ident.rms_residual < 0.6e-3                         # survey noise was 0.5 mm
    assert used.zones == case.deformation.zones                # zones kept, not double counted


def test_identification_subtracts_known_zone_part_exactly():
    from shipsupport.deformation import effective_model
    joints = (WeldJoint(30.0, 0.0), WeldJoint(70.0, 0.0))
    zones = (ShrinkageZone(0.0, 50.0, 8e-6), ShrinkageZone(50.0, 100.0, 3e-6))
    xs = np.linspace(0, 100, 21)
    truth = WeldingDeformation(joints=(WeldJoint(30.0, 0.6e-3), WeldJoint(70.0, 0.3e-3)),
                               zones=zones)
    survey = MeasuredProfile(xs, keel_lift(truth, xs) + 0.02, identify_joints=True)
    used, ident = effective_model(WeldingDeformation(joints=joints, zones=zones, measured=survey))
    assert np.allclose([j.angle for j in used.joints], [0.6e-3, 0.3e-3], rtol=1e-8)
    assert ident.rms_residual < 1e-12


def test_survey_without_identification_is_the_total_shape():
    from shipsupport.deformation import effective_model
    survey = MeasuredProfile([0, 50, 100], [5e-3, 0.0, 5e-3])
    used, ident = effective_model(WeldingDeformation(joints=(WeldJoint(50.0, 1e-3),),
                                                     measured=survey))
    assert ident is None and not used.joints and used.measured is survey


def test_bad_survey_file(tmp_path):
    from shipsupport.loader import read_survey_csv
    path = tmp_path / "s.csv"
    path.write_text("station,height\n0,1\n5,2\n", encoding="utf-8")
    with pytest.raises(InputError, match="header"):
        read_survey_csv(path)


def test_welding_example_reproduces_end_lift():
    result = run_case(case_from_dict(yaml.safe_load(
        (ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))))
    fem = result.supports["fem"]
    assert fem.reactions[0] == 0.0 and fem.reactions[-1] == 0.0      # first and last blocks lifted
    assert result.weld_gap[0] > 0.03 and result.weld_gap[-1] > 0.03
    assert np.allclose(result.springs.gap, result.weld_gap)          # no initial block gaps
    assert set(result.supports) == {"rigid", "fem"}
