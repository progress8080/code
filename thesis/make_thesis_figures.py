"""Create the thesis figures.

Two diagrams are drawn here (method flowchart, software architecture); all
result figures are produced by the program itself, through
make_result_figures.py and the welding, weight-curve and spacing functions below.

Run from the project root:
    python thesis/make_thesis_figures.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

THESIS = Path(__file__).resolve().parent
ROOT = THESIS.parent
OUT = THESIS / "figures"
sys.path.insert(0, str(ROOT / "src"))     # use the repository's package, installed or not

RESULT_FIGURES = [
    "fig_trapezoid.png", "fig_weight_distribution.png", "fig_block_spring_law.png",
    "fig_hermite_shapes.png", "fig_results_example.png", "fig_block_forces_example.png",
    "fig_spacing_study.png", "fig_ei_study.png",
    "fig_gui_weights.png", "fig_gui_blocks.png", "fig_gui_block_forces.png", "fig_gui_checks.png",
]

BLUE, BLUE_EDGE = "#dfe9f5", "#1f4e79"
GREEN, GREEN_EDGE = "#e2f0e4", "#2f6b3a"
ORANGE, ORANGE_EDGE = "#fbecd5", "#a0611a"
GREY, GREY_EDGE = "#eeeeee", "#555555"


def _box(ax, x, y, w, h, text, face, edge, size=9.5, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                facecolor=face, edgecolor=edge, lw=1.4))
    ax.text(x, y, text, ha="center", va="center", fontsize=size,
            fontweight="bold" if bold else "normal", wrap=True)


def _arrow(ax, x0, y0, x1, y1, text=None, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=14,
                                 lw=1.3, color="#333333", linestyle=ls))
    if text:
        ax.text((x0 + x1) / 2 + 0.08, (y0 + y1) / 2, text, fontsize=8, color="#333333", va="center")


def fig_method_flow():
    fig = Figure(figsize=(9.5, 11))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")

    _box(ax, 5, 11.3, 8.6, 0.9,
         "Input: weight items (LCG or moment) and weight curve, hull sections (EI, κGA, Z),\n"
         "block types, block plan (spacing, strength, gap), floors, permissible BM/SF, welding",
         GREY, GREY_EDGE)
    _box(ax, 5, 10.0, 6.0, 0.75, "Validation and unit conversion to SI\n(loader, data model)",
         BLUE, BLUE_EDGE)
    _box(ax, 5, 8.85, 9.0, 0.75, "Mesh and exact weight distribution q(x)   +   welding deformation "
         "of the assembled hull\n(joint kinks, shrinkage curvature, survey / identified angles) "
         "→ initial block gaps g_i", BLUE, BLUE_EDGE)
    _box(ax, 2.6, 7.4, 4.0, 0.95, "Method A – rigid hull\nv(x) = a + b·x  (2×2 system)",
         GREEN, GREEN_EDGE, bold=True)
    _box(ax, 7.4, 7.4, 4.0, 0.95, "Method B – elastic hull FEM\nTimoshenko beam on block springs",
         GREEN, GREEN_EDGE, bold=True)
    _box(ax, 5, 5.85, 8.0, 0.95,
         "Contact state iteration for every block: open / elastic / crushed\n"
         "R = k_t·δ + c; gap continuation for large gaps (cycle-safe)",
         ORANGE, ORANGE_EDGE)
    _box(ax, 5, 4.45, 8.0, 0.85,
         "Block reactions R_i; spacing design: R_i ≤ S_i by equidistribution → SLP\n"
         "minimax (HiGHS) → floors + 1-opt; each candidate re-analysed with its welding gaps",
         ORANGE, ORANGE_EDGE)
    _box(ax, 2.6, 3.0, 4.2, 1.05,
         "Hull girder strength\nV(x), M(x), σ deck/keel, τ\npermissible envelopes",
         BLUE, BLUE_EDGE)
    _box(ax, 7.4, 3.0, 4.2, 1.05,
         "Block checks\nbearing pressure, capacity,\nposition vs structure",
         BLUE, BLUE_EDGE)
    _box(ax, 5, 1.55, 8.0, 0.8,
         "Verdict: hull and block checks; designed spacing fits / does not fit", BLUE, BLUE_EDGE)
    _box(ax, 5, 0.45, 8.0, 0.6, "Report: summary, CSV tables, plots, GUI result tabs",
         GREY, GREY_EDGE)

    _arrow(ax, 5, 10.85, 5, 10.38)
    _arrow(ax, 5, 9.62, 5, 9.23)
    _arrow(ax, 4.2, 8.47, 2.9, 7.88)
    _arrow(ax, 5.8, 8.47, 7.1, 7.88)
    _arrow(ax, 2.6, 6.92, 3.4, 6.33)
    _arrow(ax, 7.4, 6.92, 6.6, 6.33)
    _arrow(ax, 5, 5.37, 5, 4.88)
    _arrow(ax, 4.0, 4.02, 2.9, 3.53)
    _arrow(ax, 6.0, 4.02, 7.1, 3.53)
    _arrow(ax, 8.9, 4.45, 9.55, 4.45, style="-")
    _arrow(ax, 9.55, 4.45, 9.55, 1.55, style="-")
    _arrow(ax, 9.55, 1.55, 9.0, 1.55)
    _arrow(ax, 5, 1.15, 5, 0.75)
    _arrow(ax, 2.6, 2.47, 2.6, 1.95)
    _arrow(ax, 7.4, 2.47, 7.4, 1.95)
    return fig


def fig_architecture():
    """Layered architecture: every layer only uses the layers below it."""
    fig = Figure(figsize=(11, 7.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 7.6)
    ax.axis("off")

    layers = [
        # (y, title, face, edge, [(x, width, label), ...])
        (6.55, "Layer 4 – user interfaces", GREEN, GREEN_EDGE, [
            (2.15, 2.6, "cli.py\ncommand line"),
            (6.9, 6.6, "shipsupport_gui (PyQt5)\napp · main_window · input_panels · "
                       "results_panels · worker · widgets")]),
        (4.85, "Layer 3 – orchestration and output", BLUE, BLUE_EDGE, [
            (3.4, 4.6, "analysis.py\nrun_case(): runs all steps, collects warnings"),
            (8.3, 3.8, "report.py\nsummary · CSV · plots")]),
        (3.15, "Layer 2 – numerical core", BLUE, BLUE_EDGE, [
            (1.35, 1.6, "weights.py\nload model"),
            (3.3, 1.9, "deformation.py\nwelding gaps"),
            (5.45, 2.1, "solvers/\ncontact·rigid·fem"),
            (7.6, 1.9, "spacing.py\nspacing design"),
            (9.6, 1.7, "strength.py\nV, M, checks")]),
        (1.45, "Layer 1 – input and data model", ORANGE, ORANGE_EDGE, [
            (2.6, 3.4, "loader.py\nYAML / dict → Case"),
            (6.0, 2.8, "model.py\nvalidated dataclasses"),
            (8.9, 2.4, "units.py\nSI constants")]),
    ]
    for y, title, face, edge, boxes in layers:
        ax.add_patch(FancyBboxPatch((0.3, y - 0.62), 10.4, 1.24, boxstyle="round,pad=0.02",
                                    facecolor="white", edgecolor=edge, lw=1.0, ls="--"))
        ax.text(0.45, y + 0.5, title, fontsize=9, fontweight="bold", color=edge, va="center")
        for x, w, label in boxes:
            _box(ax, x, y - 0.1, w, 0.78, label, face, edge, size=8.8)

    for x in (2.15, 6.9):          # interfaces → orchestration
        _arrow(ax, x, 5.93, x if x < 5 else 4.6, 5.53)
    _arrow(ax, 8.3, 5.93, 8.3, 5.53)
    for x in (1.35, 3.3, 5.45, 7.6, 9.6):  # orchestration → core
        _arrow(ax, 3.4, 4.23, x, 3.83)
    for x in (2.6, 6.0, 8.9):        # core → data model
        _arrow(ax, x, 2.53, x, 2.13)
    ax.text(0.45, 0.3, "Arrows point from a module to the modules it uses. The calculation "
            "package (layers 1–3 and cli.py) has no GUI dependency;\nthe GUI builds the same "
            "input dictionary as the YAML file and calls loader.case_from_dict → "
            "analysis.run_case → report drawing functions.",
            fontsize=8.3, color="#333333", va="center")
    return fig


# --------------------------------------------------------------------------- welding, weights, spacing

def _welding_data():
    import yaml
    return yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))


def _run(data, methods=None):
    from shipsupport import case_from_dict, run_case
    return run_case(case_from_dict(data), methods)


def fig_welding_data():
    """Welding data figure of the survey-based sample case (identified joint angles)."""
    from shipsupport import load_case, run_case
    from shipsupport.report import draw_welding_data
    fig = Figure(figsize=(13, 11.5))
    draw_welding_data(fig, run_case(load_case(ROOT / "examples" / "welding_survey_case.yaml")))
    fig.savefig(OUT / "fig_welding_data.png", dpi=130)


def fig_identification():
    import numpy as np
    from shipsupport.deformation import identify_joint_angles, kink_shape
    from shipsupport.model import MeasuredProfile
    data = _welding_data()
    joints = [j["x"] for j in data["welding_deformation"]["joints"]]
    true = np.array([j["angle_mrad"] for j in data["welding_deformation"]["joints"]]) * 1e-3
    xs = np.arange(0.0, 120.01, 5.0)
    rng = np.random.default_rng(2024)
    clean = sum(a * kink_shape(xs, xj) for a, xj in zip(true, joints, strict=True)) + 0.05 - 3e-4 * xs
    survey = clean + rng.normal(0.0, 1.0e-3, xs.size)                 # 1 mm survey noise
    ident = identify_joint_angles(MeasuredProfile(xs, survey), joints)
    est = np.array([j.angle for j in ident.joints])
    xc = np.linspace(0, 120, 500)
    fitted_shape = sum(a * kink_shape(xc, xj) for a, xj in zip(est, joints, strict=True))
    # align for plotting: remove the fitted straight line from the survey
    A = np.column_stack([np.ones_like(xs), xs])
    line = np.linalg.lstsq(A, survey - sum(a * kink_shape(xs, xj) for a, xj in zip(est, joints, strict=True)),
                           rcond=None)[0]
    fig = Figure(figsize=(11, 4.6))
    ax1, ax2 = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.6, 1.0]})
    ax1.plot(xs, (survey - line[0] - line[1] * xs) * 1e3, "s", mfc="none", color="#2f6b3a",
             label="synthetic survey (1 mm noise)")
    ax1.plot(xc, fitted_shape * 1e3, color="#1f4e79", lw=2, label="identified kink model")
    for xj in joints:
        ax1.axvline(xj, color="0.6", ls=":", lw=1)
    ax1.set_xlabel("x from AP [m]")
    ax1.set_ylabel("keel lift [mm]")
    ax1.set_title(f"Survey fit (rms residual {ident.rms_residual * 1e3:.2f} mm)")
    ax1.legend(fontsize=8)
    width = 0.35
    pos = np.arange(len(joints))
    ax2.bar(pos - width / 2, true * 1e3, width, label="true", color="#b5b5b5")
    ax2.bar(pos + width / 2, est * 1e3, width, label="identified", color="#1f4e79")
    ax2.set_xticks(pos, [f"{x:g} m" for x in joints])
    ax2.set_ylabel("joint angle [mrad]")
    ax2.set_title("Joint angles")
    ax2.legend(fontsize=8)
    for a in (ax1, ax2):
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig_identification.png", dpi=130)


def fig_gui_welding():
    import os
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    if os.name == "nt":
        os.environ.setdefault("QT_QPA_FONTDIR", r"C:\Windows\Fonts")
    from PyQt5.QtWidgets import QApplication
    from shipsupport_gui.main_window import MainWindow
    from shipsupport_gui.style import apply_style
    from shipsupport_gui.worker import run_analysis
    app = QApplication.instance() or QApplication([])
    apply_style(app)
    w = MainWindow()
    w.resize(1600, 980)
    w.open_path(ROOT / "examples" / "welding_case.yaml")
    w.show()
    app.processEvents()
    w.show_result(run_analysis(w.editor.to_dict()))
    w.editor.setCurrentIndex(w.editor.indexOf(w.editor.welding_tab))
    w.results.setCurrentIndex(w.results.indexOf(w.results.weld_data_plot))
    app.processEvents()
    w.grab().save(str(OUT / "fig_gui_welding.png"))
    w.modified = False
    w.close()


# Run times measured for Section 5.8 (median of repeated runs; Python 3.14.6, NumPy 2.4.6,
# SciPy 1.18.0, Intel64 Family 6 Model 158, 8 logical CPUs, Windows 11). Re-measure on the
# current machine with:  python thesis/make_thesis_figures.py --measure-scaling
SCALING = {
    "n_e": [96, 147, 255, 504, 1200, 2400],
    "fem_ms": [1.32, 1.34, 1.68, 2.22, 3.85, 6.32],
    "strength_ms": [2.84, 4.03, 7.63, 12.68, 29.26, 60.03],
    "run_ms": [5.36, 6.81, 10.02, 16.13, 34.56, 67.54],
}


def _measure_scaling():
    """Re-measure the run times of Section 5.8 versus the mesh size (takes about a minute)."""
    import copy
    import statistics
    import time

    import yaml
    from shipsupport import case_from_dict
    from shipsupport.analysis import run_case
    from shipsupport.solvers import BlockSprings, solve_fem
    from shipsupport.strength import compute_strength
    from shipsupport.weights import build_load_model, build_nodes

    def timed(f, reps):
        times = []
        for _ in range(reps):
            t = time.perf_counter()
            out = f()
            times.append(time.perf_counter() - t)
        return statistics.median(times), out

    def prepare(case):
        ship = case.ship
        blocks = sorted(case.blocks, key=lambda b: b.x)
        keys = ([w.xa for w in ship.weights] + [w.xb for w in ship.weights]
                + [s.xa for s in ship.sections] + [s.xb for s in ship.sections] + [b.x for b in blocks])
        load = build_load_model(ship.weights, build_nodes(keys, ship.x_aft, ship.x_fwd,
                                                          case.settings.max_element_length))
        return blocks, load, BlockSprings.from_blocks(blocks)

    data = {key: [] for key in SCALING}
    base = yaml.safe_load((ROOT / "examples" / "example_case.yaml").read_text(encoding="utf-8"))
    for h in (2.0, 1.0, 0.5, 0.25, 0.1, 0.05):
        d = copy.deepcopy(base)
        d["settings"].update(max_element_length=h, methods=["fem"])
        case = case_from_dict(d)
        blocks, load, springs = prepare(case)
        # loop variables are bound as default arguments so each lambda is self-contained
        t_fem, res = timed(lambda c=case, ld=load, sp=springs: solve_fem(ld, c.ship.sections, sp),
                           7)
        data["n_e"].append(load.nodes.size - 1)
        data["fem_ms"].append(1e3 * t_fem)
        data["strength_ms"].append(1e3 * timed(
            lambda c=case, ld=load, r=res, b=blocks: compute_strength(ld, r, b, c.ship, c.settings),
            5)[0])
        data["run_ms"].append(1e3 * timed(lambda c=case: run_case(c), 5)[0])
    return data


def fig_scaling(measure: bool = False):
    """Log–log run times versus the number of elements (Section 5.8)."""
    import numpy as np

    data = _measure_scaling() if measure else SCALING
    fig = Figure(figsize=(7, 4.6))
    ax = fig.subplots()
    ne = np.asarray(data["n_e"], dtype=float)
    ax.loglog(ne, data["fem_ms"], "o-", label="contact solution (solve_fem)")
    ax.loglog(ne, data["strength_ms"], "s-", label="strength post-processing")
    ax.loglog(ne, data["run_ms"], "^-", label="complete analysis (run_case)")
    ax.plot(ne, data["strength_ms"][0] * ne / ne[0], color="0.6", ls=":", lw=1.2,
            label="slope 1 (linear)")
    ax.set_xlabel("number of elements $n_e$")
    ax.set_ylabel("run time [ms]")
    ax.set_title("Mesh refinement, example vessel")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_scaling.png", dpi=130)
    return data


def fig_weight_curve():
    """Distributed load of examples/weight_curve_case.yaml: the lightship weight curve
    (with its jumps), the deadweight items given by their moments, and the total."""
    import numpy as np
    from shipsupport import load_case
    from shipsupport.units import to_tonnes
    from shipsupport.weights import build_load_model, build_nodes

    case = load_case(ROOT / "examples" / "weight_curve_case.yaml")
    ship = case.ship
    curve = ship.weight_curve
    segments = curve.items()
    names = {w.name for w in segments}
    items = [w for w in ship.weights if w.name not in names]
    keys = [w.xa for w in ship.weights] + [w.xb for w in ship.weights]
    nodes = build_nodes(keys, ship.x_aft, ship.x_fwd, 0.25)
    total = build_load_model(ship.weights, nodes)
    curve_load = build_load_model(segments, nodes)

    fig = Figure(figsize=(12, 5.2))
    ax = fig.subplots()
    x, w = total.intensity_curve()
    xc, wc = curve_load.intensity_curve()
    ax.fill_between(x, to_tonnes(w), color="#9fb7d0", alpha=0.6, label="total w(x)")
    ax.plot(xc, to_tonnes(wc), color="#1f4e79", lw=1.6,
            label=f"{curve.name}: {to_tonnes(curve_load.total):.0f} t, LCG {curve_load.lcg:.2f} m")
    ax.plot(curve.x, to_tonnes(curve.q), "o", ms=4, color="#1f4e79", label="curve points")
    colours = ["#e8a33d", "#4c9a5f", "#7b5ea7", "#d64545"]
    for item, colour in zip(items, colours * len(items)):
        moment = to_tonnes(item.weight * item.lcg)
        label = (f"{item.name}: {to_tonnes(item.weight):.0f} t, moment {moment:.0f} t·m "
                 f"→ LCG {item.lcg:.2f} m")
        if item.is_point:
            top = float(np.max(to_tonnes(w))) * 0.55
            ax.annotate("", xy=(item.lcg, 0), xytext=(item.lcg, top),
                        arrowprops=dict(arrowstyle="-|>", color=colour, lw=1.8))
            ax.plot([], [], color=colour, lw=1.8, label=label + " (point)")
        else:
            xi, wi = build_load_model([item], nodes).intensity_curve()
            ax.plot(xi, to_tonnes(wi), color=colour, lw=1.4, ls="--", label=label)
    ax.axvline(total.lcg, color="k", ls=":", lw=1)
    ax.text(total.lcg + 0.8, float(np.max(to_tonnes(w))) * 1.02,
            f"W = {to_tonnes(total.total):.0f} t, LCG = {total.lcg:.2f} m", fontsize=9)
    ax.set_xlim(ship.x_aft, ship.x_fwd)
    ax.set_ylim(0, float(np.max(to_tonnes(w))) * 1.15)
    ax.set_xlabel("x from AP [m]")
    ax.set_ylabel("w(x) [t/m]")
    ax.set_title("Weight distribution from a lightship weight curve and items given by moment")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig_weight_curve.png", dpi=130)


def fig_spacing_design():
    """Spacing design for the block plan of examples/spacing_case.yaml."""
    from shipsupport import load_case, run_case
    from shipsupport.report import draw_spacing
    fig = Figure(figsize=(13, 11))
    draw_spacing(fig, run_case(load_case(ROOT / "examples" / "spacing_case.yaml"), ["fem"]))
    fig.savefig(OUT / "fig_spacing_design.png", dpi=130)


def fig_spacing_welding():
    """Spacing design of the welded hull (examples/welding_case.yaml with the spacing design
    switched on): every candidate is analysed with the welding gaps at its block positions."""
    from shipsupport.report import draw_spacing
    data = _welding_data()
    data["spacing_design"]["enabled"] = True
    fig = Figure(figsize=(13, 11))
    draw_spacing(fig, _run(data, ["fem"]))
    fig.savefig(OUT / "fig_spacing_welding.png", dpi=130)


def fig_welding_spacing_study():
    """Largest reaction/strength at the given and at the designed spacing for welding
    deformations scaled from 0 to 1 × examples/welding_case.yaml."""
    import copy

    import numpy as np
    base = _welding_data()
    factors = [0.0, 0.1, 0.25, 0.4, 0.5, 0.55, 0.6, 0.75, 1.0]
    lift, given, designed = [], [], []
    for f in factors:
        data = copy.deepcopy(base)
        for joint in data["welding_deformation"]["joints"]:
            joint["angle_mrad"] *= f
        for zone in data["welding_deformation"]["zones"]:
            zone["deck_shrinkage_mm_per_m"] *= f
            zone["bottom_shrinkage_mm_per_m"] *= f
        data["spacing_design"]["enabled"] = True
        result = _run(data, ["fem"])
        lift.append(1000 * float(np.max(result.weld_gap)) if result.weld_gap is not None else 0.0)
        given.append(float(np.max(result.spacing.initial_utilisation)))
        designed.append(float(np.max(result.spacing.final_utilisation)))
    fig = Figure(figsize=(9, 5))
    ax = fig.add_subplot(111)
    ax.plot(lift, given, "D-", color="#222222", label="given spacing (2.4 m)")
    ax.plot(lift, designed, "o-", color="#1f4e79", label="designed spacing")
    ax.axhline(1.0, color="#d64545", ls="--", lw=1, label="block strength")
    ax.set_xlabel("welding end lift of the keel [mm]")
    ax.set_ylabel("max reaction / strength [–]")
    ax.set_title("Spacing design under welding deformation (44 blocks, elastic hull)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig_welding_spacing_study.png", dpi=130)
    return list(zip(factors, lift, given, designed, strict=True))


def main(argv):
    OUT.mkdir(parents=True, exist_ok=True)
    missing = [f for f in RESULT_FIGURES if not (OUT / f).exists()]
    if missing or "--regenerate" in argv:
        args = [sys.executable, str(THESIS / "make_result_figures.py")]
        subprocess.run(args + (["--no-gui"] if "--no-gui" in argv else []), check=True)
    fig_method_flow().savefig(OUT / "fig_method_flow.png", dpi=140)
    fig_architecture().savefig(OUT / "fig_architecture.png", dpi=140)
    fig_welding_data()
    fig_identification()
    fig_weight_curve()
    fig_spacing_design()
    fig_spacing_welding()
    print("welding spacing study:", fig_welding_spacing_study())
    measured = fig_scaling(measure="--measure-scaling" in argv)
    if "--measure-scaling" in argv:
        print("measured scaling data:", measured)
    if "--no-gui" not in argv:
        fig_gui_welding()
    print(f"thesis figures ready in {OUT}")


if __name__ == "__main__":
    main(sys.argv[1:])
