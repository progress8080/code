"""Result figures of the thesis, computed with the program itself.

Called by make_thesis_figures.py; can also be run on its own from the project root:
    python thesis/make_result_figures.py        (add --no-gui to skip GUI screenshots)

All results are computed with the shipsupport package; nothing is hand-drawn.
"""

from __future__ import annotations

import copy
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = Path(__file__).resolve().parent / "figures"
EXAMPLE = ROOT / "examples" / "example_case.yaml"

import numpy as np  # noqa: E402
import yaml  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from shipsupport import case_from_dict, load_case, run_case  # noqa: E402
from shipsupport.model import WeightItem  # noqa: E402
from shipsupport.report import draw_results, plot_block_forces  # noqa: E402
from shipsupport.units import MM, to_tonnes  # noqa: E402
from shipsupport.weights import build_load_model, build_nodes, trapezoid_ends  # noqa: E402


def example_dict() -> dict:
    return yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))


def save(fig: Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=130)
    print("wrote", OUT / name)


# --------------------------------------------------------------------------- theory figures

def fig_trapezoid():
    fig = Figure(figsize=(12, 3.8))
    axes = fig.subplots(1, 3, sharey=True)
    W, xa, xb = 100.0, 0.0, 12.0
    for ax, lcg, title in zip(axes, (4.0, 5.0, 8.0),
                              ("LCG at L/3 → triangle", "LCG inside middle third → trapezoid",
                               "LCG at 2L/3 → triangle"), strict=True):
        wa, wb = trapezoid_ends(WeightItem("item", W, lcg, xa, xb))
        ax.fill_between([xa, xb], [wa, wb], color="#9fb7d0", edgecolor="#1f4e79", lw=1.5)
        ax.axvline(lcg, color="#d64545", ls="--")
        ax.text(lcg, max(wa, wb) * 1.02, f" LCG = {lcg:g} m", color="#d64545", va="bottom")
        ax.axvspan(xa + (xb - xa) / 3, xa + 2 * (xb - xa) / 3, color="0.9", zorder=0)
        ax.set_title(title, fontsize=10)
        ax.text(xa + 0.2, wa + 0.4, f"w_a = {wa:.2f} t/m", fontsize=9)
        ax.text(xb - 0.2, wb + 0.4, f"w_b = {wb:.2f} t/m", fontsize=9, ha="right")
        ax.set_xlabel("x [m]")
        ax.set_ylim(0, 20.5)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("w(x) [t/m]")
    fig.suptitle("Trapezoidal spreading of a 100 t item over 12 m (grey band = middle third)")
    fig.tight_layout()
    save(fig, "fig_trapezoid.png")


def fig_weight_distribution():
    case = load_case(EXAMPLE)
    ship = case.ship
    keys = [w.xa for w in ship.weights] + [w.xb for w in ship.weights]
    nodes = build_nodes(keys, ship.x_aft, ship.x_fwd, 0.5)
    fig = Figure(figsize=(12, 5))
    ax = fig.subplots()
    x = None
    layers, labels = [], []
    for item in ship.weights:
        load = build_load_model([item], nodes)
        if item.is_point:
            continue
        x, w = load.intensity_curve()
        layers.append(to_tonnes(w))
        labels.append(f"{item.name} ({to_tonnes(item.weight):.0f} t)")
    ax.stackplot(x, *layers, labels=labels, alpha=0.85)
    for item in ship.weights:
        if item.is_point:
            ax.annotate(f"{item.name}\n{to_tonnes(item.weight):.0f} t (point)",
                        xy=(item.lcg, 0), xytext=(item.lcg + 3, 70),
                        arrowprops=dict(arrowstyle="->"), fontsize=9)
    total = build_load_model(ship.weights, nodes)
    ax.axvline(total.lcg, color="k", ls="--")
    ax.text(total.lcg + 0.8, 82, f"LCG = {total.lcg:.2f} m\nW = {to_tonnes(total.total):.0f} t")
    ax.set_xlabel("x from AP [m]")
    ax.set_ylabel("w(x) [t/m]")
    ax.set_title("Example vessel: longitudinal weight distribution built from weight items")
    ax.legend(fontsize=8, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_weight_distribution.png")


def fig_block_spring_law():
    case = load_case(EXAMPLE)
    bt = case.blocks[0].type
    k, fy, ratio = bt.stiffness, bt.crush_force, bt.post_crush_ratio
    dy = fy / k
    d = np.linspace(-0.5e-3, 3.5e-3, 400)
    R = np.where(d <= 0, 0.0, np.where(d <= dy, k * d, k * dy + ratio * k * (d - dy)))
    fig = Figure(figsize=(8, 4.8))
    ax = fig.subplots()
    ax.plot(d / MM, to_tonnes(R), lw=2.2, color="#1f4e79")
    ax.axvline(0, color="0.4", lw=0.8)
    ax.axvline(dy / MM, color="#d64545", ls="--", lw=1)
    ax.text(-0.45, 30, "open\n(gap not closed)\nR = 0", fontsize=9)
    ax.text(0.25, to_tonnes(fy) * 0.45, f"elastic\nR = k·δ\nk = {k / 1e6:.0f} kN/mm", fontsize=9)
    ax.text(dy / MM + 0.1, to_tonnes(fy) * 0.75,
            f"crushed\nR = k₂·δ + (k − k₂)·δy\nk₂ = {ratio:g}·k", fontsize=9)
    ax.annotate(f"crush force F_y = {to_tonnes(fy):.1f} t\nat δy = {dy / MM:.2f} mm",
                xy=(dy / MM, to_tonnes(fy)), xytext=(2.4, to_tonnes(fy) * 1.18),
                arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax.set_ylim(-20, to_tonnes(fy) * 1.4)
    ax.set_xlabel("block compression δ = v − g [mm]")
    ax.set_ylabel("block reaction R [t]")
    ax.set_title("Bilinear compression-only block law (example keel block)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_block_spring_law.png")


def fig_hermite_shapes():
    xi = np.linspace(0, 1, 200)
    shapes = {
        "N₁ = 1 − 3ξ² + 2ξ³  (v₁)": 1 - 3 * xi**2 + 2 * xi**3,
        "N₂ / L = ξ − 2ξ² + ξ³  (θ₁)": xi - 2 * xi**2 + xi**3,
        "N₃ = 3ξ² − 2ξ³  (v₂)": 3 * xi**2 - 2 * xi**3,
        "N₄ / L = −ξ² + ξ³  (θ₂)": -xi**2 + xi**3,
    }
    fig = Figure(figsize=(8, 4.5))
    ax = fig.subplots()
    for label, n in shapes.items():
        ax.plot(xi, n, lw=2, label=label)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xlabel("ξ = s / L")
    ax.set_ylabel("shape function value")
    ax.set_title("Hermite cubic shape functions of the beam element")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_hermite_shapes.png")


# --------------------------------------------------------------------------- example results

def fig_example_results():
    result = run_case(load_case(EXAMPLE))
    fig = Figure(figsize=(12, 13))
    draw_results(fig, result)
    save(fig, "fig_results_example.png")
    plot_block_forces(result, OUT / "fig_block_forces_example.png")
    print("wrote", OUT / "fig_block_forces_example.png")


def _studies_base() -> dict:
    data = example_dict()
    data["blocks"]["plan"].pop("gap_mm", None)   # initial gaps belong to the original plan
    data["settings"]["methods"] = ["rigid", "fem"]
    return data


def fig_spacing_study():
    spacings = np.array([1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.8])
    rows = []
    for s in spacings:
        data = _studies_base()
        count = int(np.floor((109.6 - 6.4) / s + 1e-9)) + 1         # blocks from 6.4 to ≤ 109.6 m
        data["blocks"]["plan"].update(start_x=6.4, count=count, spacing=float(s))
        res = run_case(case_from_dict(data))
        fem, rigid = res.supports["fem"].reactions, res.supports["rigid"].reactions
        rows.append((s, fem.size, to_tonnes(fem.max()), to_tonnes(fem[0]), to_tonnes(fem[-1]),
                     to_tonnes(rigid.max()), to_tonnes(res.load.total) / fem.size))
    rows = np.array(rows)
    fig = Figure(figsize=(10, 5.5))
    ax = fig.subplots()
    ax.plot(rows[:, 0], rows[:, 2], "o-", lw=2, label="FEM – maximum block")
    ax.plot(rows[:, 0], rows[:, 3], "s--", label="FEM – aft end block")
    ax.plot(rows[:, 0], rows[:, 4], "^--", label="FEM – forward end block")
    ax.plot(rows[:, 0], rows[:, 5], "D:", color="k", label="rigid hull – maximum block")
    ax.plot(rows[:, 0], rows[:, 6], color="0.5", label="average W / N")
    ax.axhline(200, color="#d64545", ls="--", lw=1, label="block capacity (example)")
    for s, n, top in zip(rows[:, 0], rows[:, 1], rows[:, 2], strict=True):
        ax.text(s, top + 8, f"N={int(n)}", ha="center", fontsize=8)
    ax.set_xlabel("block spacing [m]  (block line 6.4 → ≤109.6 m)")
    ax.set_ylabel("block reaction [t]")
    ax.set_title("Effect of block spacing on block reactions (example vessel)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "fig_spacing_study.png")


def fig_ei_study():
    multipliers = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000.0])
    base = _studies_base()
    rigid = None
    rows = []
    for m in multipliers:
        data = copy.deepcopy(base)
        for sec in data["sections"]:
            sec["I_m4"] = sec["I_m4"] * float(m)
        res = run_case(case_from_dict(data))
        fem = res.supports["fem"].reactions
        rigid = res.supports["rigid"].reactions
        m_fem = np.max(np.abs(res.strength["fem"].M))
        rows.append((m, to_tonnes(fem[0]), to_tonnes(fem[-1]), to_tonnes(fem.max()),
                     to_tonnes(m_fem)))
    rows = np.array(rows)
    fig = Figure(figsize=(11, 5))
    ax1, ax2 = fig.subplots(1, 2)
    ax1.semilogx(rows[:, 0], rows[:, 1], "s-", label="FEM aft end block")
    ax1.semilogx(rows[:, 0], rows[:, 2], "^-", label="FEM forward end block")
    ax1.semilogx(rows[:, 0], rows[:, 3], "o-", label="FEM maximum")
    ax1.axhline(to_tonnes(rigid[0]), color="C0", ls=":", label="rigid aft end")
    ax1.axhline(to_tonnes(rigid[-1]), color="C1", ls=":", label="rigid forward end")
    ax1.set_xlabel("EI multiplier (1 = example hull)")
    ax1.set_ylabel("reaction [t]")
    ax1.set_title("End-block reactions vs hull stiffness")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3, which="both")
    ax2.semilogx(rows[:, 0], rows[:, 4], "o-", color="C3")
    ax2.set_xlabel("EI multiplier")
    ax2.set_ylabel("max |M| [t·m]")
    ax2.set_title("Hull bending moment vs hull stiffness")
    ax2.grid(alpha=0.3, which="both")
    fig.tight_layout()
    save(fig, "fig_ei_study.png")


# --------------------------------------------------------------------------- GUI screenshots

def fig_gui():
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
    w.open_path(EXAMPLE)
    w.show()
    app.processEvents()
    w.show_result(run_analysis(w.editor.to_dict()))
    shots = {"fig_gui_weights.png": (1, 0), "fig_gui_blocks.png": (3, 3),
             "fig_gui_block_forces.png": (3, 1), "fig_gui_checks.png": (6, 4)}
    for name, (editor_tab, result_tab) in shots.items():
        w.editor.setCurrentIndex(editor_tab)
        w.results.setCurrentIndex(result_tab)
        app.processEvents()
        w.grab().save(str(OUT / name))
        print("wrote", OUT / name)
    w.modified = False
    w.close()


def main(argv):
    fig_trapezoid()
    fig_weight_distribution()
    fig_block_spring_law()
    fig_hermite_shapes()
    fig_example_results()
    fig_spacing_study()
    fig_ei_study()
    if "--no-gui" not in argv:
        fig_gui()


if __name__ == "__main__":
    main(sys.argv[1:])
