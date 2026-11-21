import os
import time
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("PyQt5")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import QEventLoop, QTimer  # noqa: E402
from PyQt5.QtWidgets import QApplication, QGroupBox  # noqa: E402

from shipsupport import case_from_dict, load_case, run_case  # noqa: E402
from shipsupport.units import TONNE_FORCE  # noqa: E402
from shipsupport_gui.input_panels import CaseEditor, format_layers, parse_layers  # noqa: E402
from shipsupport_gui.main_window import EXAMPLES, MainWindow  # noqa: E402
from shipsupport_gui.plan_panel import COL_GAP, COL_SPACING, COL_STRENGTH, COL_U, COL_X  # noqa: E402
from shipsupport_gui.worker import AnalysisThread, run_analysis  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "example_case.yaml"
SPACING = ROOT / "examples" / "spacing_case.yaml"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def reference():
    return run_case(load_case(EXAMPLE))


def example_dict(path=EXAMPLE):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def assert_same_result(a, b):
    assert a.supports.keys() == b.supports.keys()
    for m in a.supports:
        assert np.allclose(a.supports[m].reactions, b.supports[m].reactions, rtol=1e-12)


def wait_for(window, app, seconds=300):
    end = time.time() + seconds
    while window.thread is not None and window.thread.isRunning() and time.time() < end:
        app.processEvents()
        time.sleep(0.02)
    for _ in range(30):
        app.processEvents()


def tab_titles(tabs):
    return [tabs.tabText(i).replace("&&", "&") for i in range(tabs.count())]


def test_layers_text_round_trip():
    layers = parse_layers("0.9@30000; 0.2@1000;0.1@400", "t")
    assert layers[2] == {"thickness": 0.1, "E_MPa": 400.0}
    assert parse_layers(format_layers(layers), "t") == layers
    with pytest.raises(ValueError):
        parse_layers("0.9-30000", "t")


def test_tabs_as_built_and_no_removed_widgets(app):
    window = MainWindow()
    assert tab_titles(window.editor) == ["Ship", "Weights", "Hull sections", "Blocks",
                                         "Structure & limits", "Welding deformation", "Settings"]
    assert tab_titles(window.results) == ["Summary", "Block forces", "Plots", "Block reactions",
                                          "Checks", "Welding data", "Spacing design", "Warnings"]
    for name in ("landing_box", "landing_form", "hydro", "optimization_box", "optimization_form",
                 "layout_box", "layout_form", "strengths", "block_groups", "gaps",
                 "layout_preview", "strength_tab"):
        assert not hasattr(window.editor, name), name
    for name in ("landing_plot", "layout_plot", "weld_plot"):
        assert not hasattr(window.results, name), name
    plan = window.editor.plan_editor
    assert not hasattr(plan, "form") and not hasattr(plan, "copy_requested")
    assert not hasattr(window, "sample_layout_action")
    titles = [b.title() for b in window.findChildren(QGroupBox)]
    assert not any("ptimisation" in t or "landing" in t.lower() or "layout" in t.lower()
                   for t in titles)
    menu_texts = [a.text() for a in window.findChildren(type(window.run_action))]
    assert not any("strength &&" in t or "layout" in t.lower() for t in menu_texts)
    assert [path.name for _, path in EXAMPLES] == ["example_case.yaml", "spacing_case.yaml",
                                                   "weight_curve_case.yaml", "welding_case.yaml",
                                                   "welding_survey_case.yaml"]
    window.modified = False
    window.close()


@pytest.mark.parametrize("name", ["example_case", "welding_case", "welding_survey_case",
                                  "spacing_case", "weight_curve_case"])
def test_every_example_round_trips_through_the_editor(app, name):
    from shipsupport.loader import inline_data_files
    path = ROOT / "examples" / f"{name}.yaml"
    data = inline_data_files(example_dict(path), path.parent)
    editor = CaseEditor()
    editor.load_dict(data)
    via_gui = case_from_dict(editor.to_dict())
    direct = load_case(path)
    assert [b.x for b in via_gui.blocks] == [b.x for b in direct.blocks]
    assert [b.effective_strength for b in via_gui.blocks] == \
        [b.effective_strength for b in direct.blocks]
    assert [(w.name, w.weight, w.lcg, w.xa, w.xb) for w in via_gui.ship.weights] == \
        [(w.name, w.weight, w.lcg, w.xa, w.xb) for w in direct.ship.weights]
    assert (via_gui.spacing is None) == (direct.spacing is None)


def test_editor_round_trip_gives_identical_results(app, reference):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    assert_same_result(run_analysis(editor.to_dict()), reference)


def test_weights_tab_items_with_moment_and_weight_curve(app):
    editor = CaseEditor()
    editor.load_dict(example_dict(ROOT / "examples" / "weight_curve_case.yaml"))
    assert editor.curve.table.rowCount() == 18 and editor.curve_name.text() == "lightship curve"
    assert "Σ weight = 4643.0 t, LCG = 56.639 m" in editor.weight_summary.text()
    assert "curve 4343.0 t" in editor.weight_summary.text()
    data = editor.to_dict()
    assert data["weight_curve"]["name"] == "lightship curve"
    assert data["weights"][0]["moment_tm"] == 14000 and "lcg" not in data["weights"][0]
    assert data["weights"][2] == {"name": "stores", "weight_t": 40.0, "moment_tm": 1200.0}
    # an item given by its LCG on the example: the moment column gives the same result
    editor = CaseEditor()
    editor.load_dict(example_dict())
    table = editor.weights.table
    headers = [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]
    assert headers == ["Name", "Weight [t]", "LCG [m]", "Moment [t·m]", "xa [m]", "xb [m]"]
    table.item(1, 2).setText("")
    table.item(1, 3).setText(str(380 * 22.0))
    assert "Σ weight = 4500.0 t, LCG = 50.469 m" in editor.weight_summary.text()
    case = case_from_dict(editor.to_dict())
    assert case.ship.weights[1].lcg == pytest.approx(22.0)
    table.item(1, 2).setText("23.0")                               # disagrees with the moment
    with pytest.raises(Exception, match="does not agree"):
        case_from_dict(editor.to_dict())


def test_blocks_tab_types_plan_and_count_spacing_strength_gap(app):
    editor = CaseEditor()
    editor.load_dict(example_dict(SPACING))
    assert editor.indexOf(editor.blocks_tab) == 3
    plan = editor.plan_editor
    assert plan.count.value() == 44 and plan.block_type.currentText() == "keel"
    assert [plan.block_type.itemText(i) for i in range(plan.block_type.count())] == ["keel"]
    plan.count.setValue(46)                                  # rows added, copying the last row
    d = plan.plan()
    assert d["count"] == 46 and d["spacing"] == 2.4 and d["strength_t"][-1] == 190.0
    plan.count.setValue(10)
    plan.all_spacing.setText("3.0")
    plan._set_all(COL_SPACING, plan.all_spacing)
    plan.all_strength.setText("250")
    plan._set_all(COL_STRENGTH, plan.all_strength)
    plan.table.item(4, COL_SPACING).setText("1.6")           # edit one spacing
    plan.table.item(2, COL_GAP).setText("0.5")               # and one initial gap
    plan.table.item(7, COL_STRENGTH).setText("")             # empty = type capacity
    case = case_from_dict(editor.to_dict())
    xs = [b.x for b in case.blocks]
    assert len(xs) == 10 and np.allclose(np.diff(xs), [3, 3, 3, 3, 1.6, 3, 3, 3, 3])
    assert [round(b.effective_strength / TONNE_FORCE) for b in case.blocks] == \
        [250] * 7 + [200] + [250] * 2
    assert case.blocks[2].gap == pytest.approx(0.5e-3)
    assert float(plan.table.item(9, COL_X).text()) == pytest.approx(xs[-1])   # x column
    assert not plan.table.item(0, COL_X).flags() & 2                          # x read-only
    plan.table.item(3, COL_SPACING).setText("abc")
    with pytest.raises(ValueError, match="B04 spacing"):
        editor.to_dict()


def test_hidden_spacing_design_round_trips(app, tmp_path):
    window = MainWindow()
    assert window.open_path(SPACING)
    loaded = example_dict(SPACING)["spacing_design"]
    assert window.editor.to_dict()["spacing_design"] == loaded
    target = tmp_path / "saved.yaml"
    assert window._write(target)
    assert example_dict(target)["spacing_design"] == loaded
    assert window.editor.auto_spacing_design() == {**loaded, "enabled": True}
    window.open_path(EXAMPLE)                               # no spacing_design in the file
    assert "spacing_design" not in window.editor.to_dict()
    assert window.editor.auto_spacing_design() == {"enabled": True, "min_spacing": 0.8}
    window.editor.structure_box.setChecked(False)
    assert window.editor.auto_spacing_design() == {"enabled": True}
    window.modified = False
    window.close()


def test_f5_does_not_run_the_spacing_design(app):
    window = MainWindow()
    assert window.open_path(SPACING)
    window.run_analysis()
    wait_for(window, app)
    assert window.result is not None and window.result.spacing is None
    plan = window.editor.plan_editor
    assert plan.table.item(0, COL_U).text() != ""           # reactions shown in the plan
    window.modified = False
    window.close()


def test_auto_calculate_spacing_writes_the_design_into_the_plan(app):
    window = MainWindow()
    window.open_path(SPACING)
    window.editor.plan_editor.auto_button.click()
    wait_for(window, app)
    sp = window.result.spacing
    assert sp is not None and sp.fits
    case = case_from_dict(window.editor.to_dict())           # the plan now holds the design
    assert np.allclose([b.x for b in case.blocks], [b.x for b in sp.final_blocks], atol=1e-6)
    assert [b.effective_strength for b in case.blocks] == \
        pytest.approx([b.effective_strength for b in sp.final_blocks])
    plan = window.editor.plan_editor
    assert "Fits" in plan.status.text()
    u = [float(plan.table.item(r, COL_U).text()) for r in range(44)]
    assert max(u) <= 1.0 and max(u) == pytest.approx(np.max(sp.final_utilisation), abs=1e-3)
    assert not window.results.spacing_plot.action.isVisibleTo(window.results.spacing_plot)
    assert window.editor.to_dict()["spacing_design"] == example_dict(SPACING)["spacing_design"]
    window.modified = False
    window.close()


def test_auto_calculate_spacing_without_settings_uses_the_floor_spacing(app):
    window = MainWindow()
    window.open_path(EXAMPLE)
    plan = window.editor.plan_editor
    plan.count.setValue(30)
    plan.all_spacing.setText("3.2")
    plan._set_all(COL_SPACING, plan.all_spacing)
    window.auto_spacing()
    wait_for(window, app)
    sp = window.result.spacing
    assert sp is not None and window.result.case.spacing.min_spacing == pytest.approx(0.8)
    assert np.all(sp.spacing >= 0.8 - 1e-6)
    window.modified = False
    window.close()


def test_survey_case_opens_with_points_in_table(app):
    window = MainWindow()
    survey_case = EXAMPLE.with_name("welding_survey_case.yaml")
    assert window.open_path(survey_case)
    assert window.editor.weld_survey.table.rowCount() == 25
    via_gui = run_analysis(window.editor.to_dict())
    direct = run_case(load_case(survey_case))
    assert np.allclose(via_gui.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-12)
    window.show_result(via_gui)
    assert not window.results.weld_data_plot.message.isVisible()
    window.modified = False
    window.close()


def test_sample_welding_button_fills_and_shows_welding_tab(app):
    window = MainWindow()
    assert window.open_path(EXAMPLE)                  # example without welding data
    window.show_result(run_analysis(window.editor.to_dict()))
    assert window.results.weld_data_plot.action.isVisibleTo(window.results.weld_data_plot)
    window.editor.load_sample_welding(0.0, 120.0)     # what the button does before running
    result = run_analysis(window.editor.to_dict())
    window.show_result(result)
    assert result.deformation is not None
    assert not window.results.weld_data_plot.message.isVisible()
    assert not window.results.weld_data_plot.action.isVisibleTo(window.results.weld_data_plot)
    direct = run_case(load_case(EXAMPLE.with_name("welding_case.yaml")))
    assert np.allclose(result.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-12, atol=1e-6)
    window.modified = False
    window.close()


def test_sample_welding_is_mapped_onto_the_hull(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.load_sample_welding(-5.0, 115.0)
    joints = [row["x"] for row in editor.weld_joints.rows()]
    zones = editor.weld_zones.rows()
    assert joints == pytest.approx([15.0, 35.0, 55.0, 75.0, 95.0])
    assert zones[0]["xa"] == pytest.approx(-5.0) and zones[-1]["xb"] == pytest.approx(115.0)


def test_unresolved_survey_file_is_reported(app):
    data = example_dict()
    data["welding_deformation"] = {"measured": {"file": "survey.csv"}}
    with pytest.raises(ValueError):
        CaseEditor().load_dict(data)


def test_bad_cell_reports_location(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.weights.table.item(1, 1).setText("abc")
    with pytest.raises(ValueError, match="row 2"):
        editor.to_dict()


def test_point_load_needs_both_or_neither_extent(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.weights.table.item(0, 5).setText("")
    with pytest.raises(ValueError, match="both xa and xb"):
        editor.to_dict()


def test_main_window_show_result_and_save(app, reference, tmp_path):
    window = MainWindow()
    assert window.open_path(EXAMPLE)
    window.show_result(run_analysis(window.editor.to_dict()))
    assert window.results.blocks.rowCount() == 44
    assert window.results.checks.rowCount() == 14
    assert window.export_action.isEnabled()
    assert "All checks passed" in window.state_label.text() or \
        "CHECKS FAILED" in window.state_label.text()
    assert "shims" not in window.state_label.text()

    target = tmp_path / "saved.yaml"
    assert window._write(target)
    assert not window.modified
    assert_same_result(run_case(load_case(target)), reference)
    window.close()


def test_open_every_example_and_run(app):
    window = MainWindow()
    for _, path in EXAMPLES:
        window.modified = False
        window.open_example(path)
        window.run_analysis()
        wait_for(window, app)
        assert window.result is not None, path
        assert window.result.spacing is None
        assert window.results.blocks.rowCount() == len(window.result.case.blocks)
    window.modified = False
    window.close()


def test_analysis_thread(app, reference):
    got = {}
    loop = QEventLoop()
    thread = AnalysisThread(example_dict())
    thread.succeeded.connect(lambda r: got.setdefault("result", r))
    thread.failed.connect(lambda msg: got.setdefault("error", msg))
    thread.finished.connect(loop.quit)
    QTimer.singleShot(30000, loop.quit)
    thread.start()
    loop.exec_()
    thread.wait()
    assert "error" not in got
    assert_same_result(got["result"], reference)


def test_analysis_thread_reports_input_error(app):
    data = example_dict()
    data["weights"][0]["lcg"] = 1.0   # outside the middle third
    got = {}
    loop = QEventLoop()
    thread = AnalysisThread(data)
    thread.failed.connect(lambda msg: got.setdefault("error", msg))
    thread.finished.connect(loop.quit)
    QTimer.singleShot(30000, loop.quit)
    thread.start()
    loop.exec_()
    thread.wait()
    assert "middle third" in got["error"]


# --------------------------------------------------------------------------- weight curve sample

def test_sample_weight_curve_matches_the_example_file(app):
    from shipsupport_gui.input_panels import SAMPLE_WEIGHT_CURVE
    data = yaml.safe_load((ROOT / "examples" / "weight_curve_case.yaml").read_text(encoding="utf-8"))
    points = [(p["x"], p["q_t_per_m"]) for p in data["weight_curve"]["points"]]
    assert SAMPLE_WEIGHT_CURVE["points"] == points
    assert [(i["name"], i["weight_t"], i["moment_tm"]) for i in SAMPLE_WEIGHT_CURVE["items"]] == \
        [(i["name"], i["weight_t"], i["moment_tm"]) for i in data["weights"]]


def test_sample_weight_curve_loads_and_runs(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    total, lcg = editor.load_sample_weight_curve()
    assert total == pytest.approx(4643.0) and lcg == pytest.approx(56.639, abs=1e-3)
    result = run_analysis(editor.to_dict())
    assert result.load.total == pytest.approx(4643.0 * TONNE_FORCE, rel=1e-9)
    assert result.case.ship.weight_curve is not None and not result.warnings
    direct = run_case(load_case(ROOT / "examples" / "weight_curve_case.yaml"))
    assert np.allclose(result.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-9)


def test_sample_weight_curve_is_mapped_onto_another_hull(app):
    editor = CaseEditor()
    data = example_dict()
    data["ship"].update(lpp=150.0, x_fwd=150.0)
    editor.load_dict(data)
    total, lcg = editor.load_sample_weight_curve()
    assert total == pytest.approx(4643.0, rel=1e-6)                 # weights kept
    assert lcg == pytest.approx(56.639 * 150.0 / 120.0, abs=2e-3)    # positions scaled
    assert editor.curve.rows()[-1]["x"] == pytest.approx(150.0)


def test_block_forces_and_spacing_charts_show_the_weight_distribution():
    from matplotlib.figure import Figure

    from shipsupport.report import draw_block_forces, draw_spacing
    result = run_case(load_case(ROOT / "examples" / "weight_curve_case.yaml"), ["fem"])
    for draw in (draw_block_forces, draw_spacing):
        fig = Figure()
        draw(fig, result)
        top = fig.axes[0]
        assert top.get_ylabel() == "load [t/m]"
        labels = [line.get_label() for line in top.get_lines()]
        assert "total distributed load w(x)" in labels
        assert any(label.startswith("weight curve") for label in labels)
