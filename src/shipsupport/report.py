"""Text summary, CSV tables and plots for an AnalysisResult."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .analysis import AnalysisResult
from .strength import Check
from .units import KN_PER_MM, MM, MPA, TONNE_FORCE, to_tonnes

try:                       # one chart style for the GUI, the exported PNGs and the thesis
    import matplotlib

    matplotlib.rcParams.update({
        "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 9.5, "axes.labelcolor": "#1d2733", "axes.edgecolor": "#9aa5b1",
        "axes.spines.top": False, "axes.spines.right": False, "axes.titlecolor": "#1d2733",
        "xtick.color": "#3d4a57", "ytick.color": "#3d4a57", "grid.color": "#d9dee4",
        "grid.linewidth": 0.7, "legend.frameon": True, "legend.framealpha": 0.92,
        "legend.edgecolor": "#cfd6df", "legend.fancybox": True, "figure.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight", "savefig.pad_inches": 0.15,
    })
except ImportError:        # matplotlib is optional for the calculation core
    pass

_DISPLAY = {"N": (TONNE_FORCE, "t"), "N·m": (TONNE_FORCE, "t·m"), "Pa": (MPA, "MPa")}


def _fmt_check(c: Check) -> str:
    scale, unit = _DISPLAY[c.unit]
    capacity = "—" if c.capacity is None else f"{c.capacity / scale:12.2f}"
    util = "—" if np.isnan(c.utilisation) else f"{c.utilisation:6.3f}"
    status = {True: "OK", False: "FAIL", None: "no limit"}[c.passed]
    return (f"  {c.name:<26} {c.demand / scale:12.2f} {capacity:>12} {unit:<4} "
            f"{util:>6}  x = {c.x:8.3f} m  {status}")


def format_summary(result: AnalysisResult) -> str:
    case, load = result.case, result.load
    ship = case.ship
    out = []
    add = out.append

    add("=" * 100)
    add(f"SHIP SUPPORT FORCE ANALYSIS — {ship.name}")
    add("=" * 100)
    add(f"Total weight        : {to_tonnes(load.total):12.1f} t")
    add(f"LCG                 : {load.lcg:12.3f} m from AP")
    curve = ship.weight_curve
    if curve is not None:
        segments = curve.items()
        w_curve = sum(w.weight for w in segments)
        lcg_curve = sum(w.weight * w.lcg for w in segments) / w_curve
        add(f"  of which curve    : {to_tonnes(w_curve):12.1f} t, LCG {lcg_curve:.3f} m "
            f"('{curve.name}', {curve.x.size} points, {curve.x[0]:g} … {curve.x[-1]:g} m)")
        add(f"  weight items      : {to_tonnes(ship.total_weight - w_curve):12.1f} t "
            f"({len(ship.weights) - len(segments)} items)")
    if ship.declared_weight is not None:
        add(f"Declared weight     : {to_tonnes(ship.declared_weight):12.1f} t")
    if ship.declared_lcg is not None:
        add(f"Declared LCG        : {ship.declared_lcg:12.3f} m")
    add(f"Mesh                : {load.nodes.size} nodes, {load.nodes.size - 1} elements")
    x = result.springs.x
    add(f"Blocks              : {x.size}, x = {x.min():.3f} … {x.max():.3f} m")

    add("")
    add("SUPPORT SOLUTIONS")
    add(f"  {'method':<11}{'ΣR [t]':>12}{'max R [t]':>12}{'at x [m]':>10}{'min R [t]':>12}"
        f"{'open':>6}{'crushed':>9}{'iter':>6}{'force res':>12}{'moment res':>12}")
    for m, r in result.supports.items():
        i = int(np.argmax(r.reactions))
        add(f"  {m:<11}{to_tonnes(r.reactions.sum()):12.1f}{to_tonnes(r.reactions[i]):12.2f}"
            f"{r.block_x[i]:10.3f}{to_tonnes(r.reactions.min()):12.2f}{r.n_open:6d}{r.n_crushed:9d}"
            f"{r.iterations:6d}{r.force_residual:12.1e}{r.moment_residual:12.1e}")

    design = result.design_method
    table = result.strength[design].blocks
    add("")
    add(f"BLOCK REACTIONS (pressure and state from '{design}')")
    head = f"  {'block':<6}{'x [m]':>9}{'gap [mm]':>10}{'k [kN/mm]':>11}"
    for m in result.block_supports:
        head += f"{'R ' + m + ' [t]':>14}"
    head += f"{'p [MPa]':>10}{'δ [mm]':>9}  state"
    add(head)
    springs = result.springs
    for i in range(x.size):
        line = (f"  {table.label[i]:<6}{x[i]:9.3f}{springs.gap[i] / MM:10.2f}"
                f"{springs.k[i] / KN_PER_MM:11.1f}")
        for r in result.block_supports.values():
            line += f"{to_tonnes(r.reactions[i]):14.2f}"
        delta = result.supports[design].delta[i]
        line += f"{table.pressure[i] / MPA:10.3f}{delta / MM:9.3f}  {table.state[i]}"
        if table.off_structure[i]:
            line += f"  (off structure {table.structure_offset[i]:.2f} m)"
        add(line)

    for m, s in result.strength.items():
        add("")
        add(f"HULL GIRDER & BLOCK CHECKS — {m}")
        tol = 1e-6 * max(np.max(np.abs(s.M)), 1e-12)   # ignore round-off at the free ends
        sag = (f"{to_tonnes(s.M.max()):10.1f} t·m at x = {s.x[np.argmax(s.M)]:.3f} m"
               if s.M.max() > tol else "      none")
        hog = (f"{to_tonnes(s.M.min()):10.1f} t·m at x = {s.x[np.argmin(s.M)]:.3f} m"
               if s.M.min() < -tol else "      none")
        add(f"  max sagging BM {sag},  max hogging BM {hog}")
        add(f"  {'check':<26} {'demand':>12} {'limit':>12} {'unit':<4} {'util':>6}")
        for c in s.checks:
            add(_fmt_check(c))
        verdict = {True: "ALL CHECKS PASSED", False: "ONE OR MORE CHECKS FAILED",
                   None: "NO LIMITS DEFINED"}[s.passed]
        add(f"  => {verdict}")

    _welding(result, add)
    _spacing_design(result, add)

    add("")
    add("WARNINGS" if result.warnings else "WARNINGS: none")
    for w in result.warnings:
        add(f"  - {w}")
    return "\n".join(out)


def _labels(result: AnalysisResult, indices) -> str:
    table = result.strength[result.design_method].blocks
    names = [table.label[i] for i in indices]
    return ", ".join(names) if names else "none"


def _spacing_design(result: AnalysisResult, add) -> None:
    sp = result.spacing
    if sp is None:
        return
    design = result.case.spacing
    xs = np.array([b.x for b in sp.final_blocks])
    add("")
    add(f"BLOCK SPACING DESIGN ({len(xs)} blocks, strengths kept, {design.method}, target "
        f"utilisation {design.target_utilisation:g})")
    add(f"  result                : {sp.status.replace('_', ' ').upper()} — {sp.message}")
    add(f"  block line            : {xs[0]:.2f} – {xs[-1]:.2f} m (given "
        f"{sp.initial_blocks[0].x:.2f} – {sp.initial_blocks[-1].x:.2f} m)")
    add(f"  spacing               : {sp.spacing.min():.2f} – {sp.spacing.max():.2f} m (given "
        f"{sp.initial_spacing.min():.2f} – {sp.initial_spacing.max():.2f} m)")
    u0, u1 = sp.initial_utilisation, sp.final_utilisation
    add(f"  max reaction/strength : {np.max(u0):.3f} → {np.max(u1):.3f}   (mean "
        f"{np.mean(u0):.3f} → {np.mean(u1):.3f})")
    add(f"  total strength        : {to_tonnes(sp.total_strength):.1f} t for a load of "
        f"{to_tonnes(sp.total_load):.1f} t")
    phases = {}
    for h in sp.history:
        phases[h.phase] = h
    add("  phases                : " + ", ".join(
        f"{name} {h.max_utilisation:.3f}" for name, h in phases.items())
        + f"  ({sp.evaluations} analyses)")
    add(f"  {'block':>7}{'x [m]':>10}{'spacing to next [m]':>21}{'strength [t]':>14}"
        f"{'R [t]':>10}{'R/S':>8}")
    final = sp.final_evaluation.reactions
    gaps = list(sp.spacing) + [None]
    for b, r, u, g in zip(sp.final_blocks, final, u1, gaps, strict=True):
        add(f"  {b.label:>7}{b.x:>10.2f}{'' if g is None else f'{g:.2f}':>21}"
            f"{to_tonnes(b.effective_strength):>14.1f}{to_tonnes(r):>10.1f}{u:>8.3f}"
            + ("  > target" if u > sp.target * (1 + 1e-9) else ""))
    if result.deformation is not None:
        gap = sp.final_evaluation.springs.gap
        lifted = [b.label for b, r in zip(sp.final_blocks, final, strict=True) if r <= 0]
        add(f"  welding gaps          : recomputed at the designed positions, up to "
            f"{gap.max() / MM:.1f} mm; blocks without contact: {', '.join(lifted) or 'none'}")


def _welding(result: AnalysisResult, add) -> None:
    defo, gap = result.deformation, result.weld_gap
    if defo is None or gap is None:
        return
    add("")
    add("WELDING DEFORMATION OF THE ASSEMBLED HULL")
    parts = []
    if defo.joints:
        parts.append(f"{len(defo.joints)} joint kinks")
    if defo.zones:
        parts.append(f"{len(defo.zones)} shrinkage zones")
    if defo.measured is not None:
        parts.append(f"surveyed profile ({defo.measured.x.size} points, "
                     f"{defo.measured.interpolation})")
    add(f"  model                 : {', '.join(parts)}")
    ident = result.identification
    if ident is not None:
        angles = ", ".join(f"{j.angle * 1e3:.3f}" for j in ident.joints)
        add(f"  identified angles     : [{angles}] mrad at joints "
            f"x = {', '.join(f'{j.x:g}' for j in ident.joints)} m")
        add(f"  survey fit residual   : rms {ident.rms_residual / MM:.2f} mm, "
            f"max {ident.max_residual / MM:.2f} mm")
    add(f"  keel lift at blocks   : aft end {gap[0] / MM:.1f} mm, forward end "
        f"{gap[-1] / MM:.1f} mm, max {gap.max() / MM:.1f} mm (relative to the block line)")
    design = result.supports[result.design_method]
    add(f"  blocks without contact ({result.design_method}, as built): "
        f"{_labels(result, np.flatnonzero(design.reactions <= 0))}")


def write_csv(result: AnalysisResult, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    design = result.design_method
    table = result.strength[design].blocks
    path = out_dir / "blocks.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        weld = result.weld_gap if result.weld_gap is not None else np.zeros(table.x.size)
        w.writerow(["block", "x_m", "gap_mm", "weld_gap_mm", "k_kN_per_mm", "strength_t"]
                   + [f"R_{m}_t" for m in result.block_supports]
                   + ["pressure_MPa", "pressure_util", "capacity_util", "delta_mm", "state",
                      "structure_offset_m"])
        blocks = sorted(result.case.blocks, key=lambda b: b.x)
        for i in range(table.x.size):
            s = blocks[i].effective_strength
            w.writerow([table.label[i], f"{table.x[i]:.4f}", f"{result.springs.gap[i] / MM:.3f}",
                        f"{weld[i] / MM:.3f}", f"{result.springs.k[i] / KN_PER_MM:.3f}",
                        "" if s is None else f"{to_tonnes(s):.3f}"]
                       + [f"{to_tonnes(r.reactions[i]):.4f}" for r in result.block_supports.values()]
                       + [f"{table.pressure[i] / MPA:.4f}", f"{table.pressure_util[i]:.4f}",
                          f"{table.capacity_util[i]:.4f}",
                          f"{result.supports[design].delta[i] / MM:.4f}", table.state[i],
                          f"{table.structure_offset[i]:.4f}"])
    written.append(path)

    for m, s in result.strength.items():
        path = out_dir / f"curves_{m.replace('+', '_')}.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["x_m", "V_t", "M_tm", "sigma_deck_MPa", "sigma_keel_MPa", "tau_MPa",
                        "bm_util", "sf_util"])
            for row in zip(s.x, to_tonnes(s.V), to_tonnes(s.M), s.sigma_deck / MPA,
                           s.sigma_keel / MPA, s.tau / MPA, s.bm_util, s.sf_util):
                w.writerow([f"{v:.6g}" for v in row])
        written.append(path)

    sp = result.spacing
    if sp is not None:
        path = out_dir / "spacing.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["block", "x_m", "spacing_to_next_m", "given_x_m", "given_spacing_m",
                        "strength_t", "R_t", "utilisation", "given_R_t", "given_utilisation"])
            gaps, gaps0 = list(sp.spacing) + [None], list(sp.initial_spacing) + [None]
            for i, b in enumerate(sp.final_blocks):
                w.writerow([b.label, f"{b.x:.4f}", "" if gaps[i] is None else f"{gaps[i]:.4f}",
                            f"{sp.initial_blocks[i].x:.4f}",
                            "" if gaps0[i] is None else f"{gaps0[i]:.4f}",
                            f"{to_tonnes(b.effective_strength):.3f}",
                            f"{to_tonnes(sp.final_evaluation.reactions[i]):.4f}",
                            f"{sp.final_utilisation[i]:.4f}",
                            f"{to_tonnes(sp.initial_evaluation.reactions[i]):.4f}",
                            f"{sp.initial_utilisation[i]:.4f}"])
        written.append(path)
        path = out_dir / "spacing_history.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["phase", "iteration", "analyses", "max_utilisation", "min_utilisation"])
            for h in sp.history:
                w.writerow([h.phase, h.iteration, h.evaluations, f"{h.max_utilisation:.6g}",
                            f"{h.min_utilisation:.6g}"])
        written.append(path)
    return written


def plot_results(result: AnalysisResult, path: Path) -> Path:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(12, 13))
    draw_results(fig, result)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    return path


def _layout(fig, h_pad: float = 0.06) -> None:
    """Constrained layout: recomputed on every draw, so the chart also fits after the GUI
    resizes the canvas (tight_layout would only fit the size at the time of the call)."""
    fig.set_layout_engine("constrained")
    fig.get_layout_engine().set(h_pad=h_pad, w_pad=0.06, hspace=0.04, wspace=0.04)


def draw_results(fig, result: AnalysisResult) -> None:
    """Draw load/reactions, settlement, V and M into a matplotlib Figure."""
    ship = result.case.ship
    fig.clear()
    _layout(fig)
    axes = fig.subplots(4, 1, sharex=True)

    ax = axes[0]
    xw, w = result.load.intensity_curve()
    ax.fill_between(xw, to_tonnes(w), color="0.75", label="weight w(x)")
    ax.set_ylabel("w [t/m]")
    ax.set_title(f"{ship.name} — docking analysis")
    bars = ax.twinx()
    bars.spines["right"].set_visible(True)
    n = len(result.block_supports)
    spacing = np.min(np.diff(np.unique(result.springs.x))) if result.springs.x.size > 1 else 1.0
    width = 0.8 * spacing / n
    for j, (m, r) in enumerate(result.block_supports.items()):
        bars.bar(r.block_x + (j - (n - 1) / 2) * width, to_tonnes(r.reactions), width,
                 label=f"R {m}", alpha=0.8)
    bars.set_ylabel("block reaction [t]")
    lines = ax.get_legend_handles_labels()
    more = bars.get_legend_handles_labels()
    ax.legend(lines[0] + more[0], lines[1] + more[1], loc="upper right")

    ax = axes[1]
    for m, r in result.supports.items():
        ax.plot(r.nodes, r.v / MM, label=m)
    ax.invert_yaxis()
    ax.set_ylabel("settlement v [mm] (down +)")
    ax.legend()

    ax = axes[2]
    for m, s in result.strength.items():
        ax.plot(s.x, to_tonnes(s.V), label=m)
    if ship.permissible_sf is not None:
        e = ship.permissible_sf
        ax.plot(e.x, to_tonnes(e.limit), "k--", lw=1, label="permissible")
        ax.plot(e.x, -to_tonnes(e.limit), "k--", lw=1)
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_ylabel("shear force V [t]")
    ax.legend()

    ax = axes[3]
    for m, s in result.strength.items():
        ax.plot(s.x, to_tonnes(s.M), label=m)
    if ship.permissible_bm is not None:
        e = ship.permissible_bm
        ax.plot(e.x, to_tonnes(e.sag), "k--", lw=1, label="permissible (sag +, hog −)")
        ax.plot(e.x, -to_tonnes(e.hog), "k--", lw=1)
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_ylabel("bending moment M [t·m] (sag +)")
    ax.set_xlabel("x from AP [m]")
    ax.legend()

    for ax in axes:
        ax.grid(True, alpha=0.3)


METHOD_NAMES = {"fem": "elastic hull FEM", "rigid": "rigid hull",
                "fem+spacing": "elastic hull FEM, designed block spacing",
                "rigid+spacing": "rigid hull, designed block spacing"}
_OK, _NEAR, _OVER, _OPEN = "#4c9a5f", "#e8a33d", "#d64545", "#b5b5b5"


def draw_weight_panel(ax, result: AnalysisResult, block_x=None, other_x=None,
                      other_label: str = "", block_label: str = "blocks") -> None:
    """Weight distribution along the hull: the total distributed load w(x) [t/m] (items and
    weight curve), the weight curve on its own, point loads, W and LCG, and the blocks."""
    ship = result.case.ship
    load = result.load
    xw, w = load.intensity_curve()
    w_t = to_tonnes(w)
    ax.fill_between(xw, 0.0, w_t, color="#9fb7d0", alpha=0.55, lw=0, zorder=1)
    ax.plot(xw, w_t, color="#1f4e79", lw=1.4, zorder=3, label="total distributed load w(x)")
    curve = ship.weight_curve
    top = float(np.max(w_t)) if w_t.size else 1.0
    if curve is not None:
        ax.plot(curve.x, to_tonnes(curve.q), color="#b5651d", lw=1.3, ls="--", zorder=4,
                label=f"weight curve ({curve.name})" if curve.name else "weight curve")
        top = max(top, float(np.max(to_tonnes(curve.q))))
    nodes, point = load.nodes, to_tonnes(load.point)
    for xp, p_t in zip(nodes[point > 0], point[point > 0], strict=True):
        ax.annotate(f"{p_t:.0f} t", xy=(xp, 0.0), xytext=(xp, 0.82 * top), ha="center",
                    fontsize=7.5, color="#6a3d9a",
                    arrowprops=dict(arrowstyle="-|>", color="#6a3d9a", lw=1.0))
    if np.any(point > 0):
        ax.plot([], [], color="#6a3d9a", marker=r"$\downarrow$", ls="", label="point load")
    ax.axvline(load.lcg, color="#d64545", lw=1.0, ls=":", zorder=2)
    ax.text(load.lcg, 0.96, f"W = {to_tonnes(load.total):.1f} t   LCG = {load.lcg:.2f} m",
            transform=ax.get_xaxis_transform(), fontsize=8.5, color="#a11d1d", ha="center",
            va="top", zorder=6,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e3b3b3", lw=0.8))
    if other_x is not None:
        ax.plot(other_x, np.full(len(other_x), -0.07 * top), "^", ms=4, color="#9fb7d0",
                mec="#5a6a7a", mew=0.4, clip_on=False, zorder=5, label=other_label)
    if block_x is not None:
        ax.plot(block_x, np.full(len(block_x), -0.07 * top if other_x is None else -0.15 * top),
                "^", ms=4.5, color="#1f4e79", clip_on=False, zorder=5, label=block_label)
    ax.set_ylim(-0.2 * top, 1.18 * top)
    ax.set_xlim(ship.x_aft, ship.x_fwd)
    ax.axhline(0.0, color="#333333", lw=0.8)
    ax.set_ylabel("load [t/m]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right", ncol=2, framealpha=0.92)


def plot_block_forces(result: AnalysisResult, path: Path, method: str | None = None) -> Path:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(18, 12.5))
    draw_block_forces(fig, result, method)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    return path


def draw_block_forces(fig, result: AnalysisResult, method: str | None = None) -> None:
    """One chart with the force on every block: hull schematic plus labelled bars,
    coloured by utilisation (the worse of load/capacity and pressure/allowable)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Polygon, Rectangle

    method = method or result.design_method
    support = result.supports[method]
    table = result.strength[method].blocks
    ship = result.case.ship
    blocks = sorted(result.case.blocks, key=lambda b: b.x)

    x = support.block_x
    R = to_tonnes(support.reactions)
    n = x.size
    with np.errstate(invalid="ignore"):
        util = np.fmax(table.capacity_util, table.pressure_util)
    colours = np.where(R <= 0, _OPEN, np.where(util > 1.0, _OVER, np.where(util > 0.9, _NEAR, _OK)))
    font = max(5.0, min(8.5, 360.0 / n))
    spacing = np.min(np.diff(np.unique(x))) if np.unique(x).size > 1 else 1.0
    r_max = max(R.max(), 1e-9)

    fig.clear()
    _layout(fig)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.35, 1.0, 3.6], hspace=0.06)
    axw = fig.add_subplot(gs[0])
    ax0 = fig.add_subplot(gs[1], sharex=axw)
    ax = fig.add_subplot(gs[2], sharex=axw)
    draw_weight_panel(axw, result, x)
    axw.tick_params(labelbottom=False)

    # Hull schematic with the blocks underneath
    x0, x1 = ship.x_aft, ship.x_fwd
    L = x1 - x0
    hull = [(x0, 1.0), (x1, 1.0), (x1 - 0.05 * L, 0.0), (x0 + 0.07 * L, 0.0), (x0, 0.45)]
    ax0.add_patch(Polygon(hull, closed=True, facecolor="#dfe6ee", edgecolor="#4a5a6a", lw=1.2))
    for b, xi, c in zip(blocks, x, colours):
        width = min(b.type.contact_length, 0.9 * spacing)
        ax0.add_patch(Rectangle((xi - width / 2, -0.42), width, 0.40 - 0.02, facecolor=c,
                                edgecolor="#333333", lw=0.5))
    ax0.axhline(-0.42, color="#555555", lw=2)
    ax0.text(x0 + 0.01 * L, 0.72, "AFT", fontsize=9, color="#4a5a6a", va="center")
    ax0.text(x1 - 0.01 * L, 0.72, "FWD", fontsize=9, color="#4a5a6a", va="center", ha="right")
    ax0.annotate("", xy=(result.load.lcg, 0.0), xytext=(result.load.lcg, 0.95),
                 arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=1.8))
    ax0.text(result.load.lcg + 0.006 * L, 0.55, f"W = {to_tonnes(result.load.total):.1f} t\n"
             f"LCG = {result.load.lcg:.2f} m", fontsize=9, color="#1f4e79", va="center")
    ax0.set_ylim(-0.55, 1.15)
    ax0.set_yticks([])
    ax0.tick_params(labelbottom=False)
    for side in ("top", "right", "left"):
        ax0.spines[side].set_visible(False)
    axw.set_title(f"{ship.name} — weight distribution and block forces\n"
                  f"{METHOD_NAMES.get(method, method)}: "
                  f"{n} blocks, ΣR = {R.sum():.1f} t, max {R.max():.1f} t at "
                  f"{table.label[int(np.argmax(R))]}, min {R.min():.1f} t at "
                  f"{table.label[int(np.argmin(R))]}", fontsize=12)

    # Labelled bars
    others = [m for m in result.block_supports if m != method]
    label_base = np.max([R] + [to_tonnes(result.supports[m].reactions) for m in others], axis=0)
    ax.bar(x, R, width=0.7 * spacing, color=colours, edgecolor="#333333", lw=0.5, zorder=3)
    for xi, ri, yi, label in zip(x, R, label_base, table.label):
        ax.text(xi, yi + 0.035 * r_max, f"{ri:.1f}" if ri > 0 else "0 (open)", rotation=90,
                ha="center", va="bottom", fontsize=font, zorder=4)
        ax.text(xi, -0.015 * r_max, label, rotation=90, ha="center", va="top", fontsize=font,
                color="#333333")

    handles = [Patch(color=_OK, label="≤ 90 % of limit"), Patch(color=_NEAR, label="90–100 % of limit"),
               Patch(color=_OVER, label="> 100 % of limit"), Patch(color=_OPEN, label="no contact")]
    for other in others:
        ax.plot(x, to_tonnes(result.supports[other].reactions), "D", ms=4, color="#222222",
                zorder=5)
        handles.append(Line2D([], [], marker="D", ls="", color="#222222", ms=4,
                              label=f"{METHOD_NAMES.get(other, other)} (comparison)"))
    capacity = np.array([np.nan if b.effective_strength is None else to_tonnes(b.effective_strength)
                         for b in blocks])
    if np.isfinite(capacity).any():
        ax.step(x, capacity, where="mid", color=_OVER, ls="--", lw=1.2, zorder=2)
        handles.append(Line2D([], [], color=_OVER, ls="--", label="block capacity"))
    mean = R.sum() / n
    ax.axhline(mean, color="#1f4e79", ls=":", lw=1.2, zorder=2)
    handles.append(Line2D([], [], color="#1f4e79", ls=":", label=f"average W/N = {mean:.1f} t"))

    top = max(label_base.max(), np.nanmax(capacity) if np.isfinite(capacity).any() else 0.0)
    ax.set_ylim(-0.11 * top, 1.42 * top)   # headroom for value labels and the legend
    ax.axhline(0.0, color="#333333", lw=0.8)
    ax.set_xlim(x0 - 0.01 * L, x1 + 0.01 * L)
    ax.set_ylabel("block reaction [t]")
    ax.set_xlabel("x from AP [m]")
    ax.grid(True, axis="y", alpha=0.3, zorder=0)
    ax.legend(handles=handles, loc="upper center", ncol=min(len(handles), 4), fontsize=9,
              framealpha=0.95)


def plot_welding_data(result: AnalysisResult, path: Path) -> Path | None:
    if result.deformation is None:
        return None
    from matplotlib.figure import Figure

    fig = Figure(figsize=(13, 11.5))
    draw_welding_data(fig, result)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    return path


def draw_welding_data(fig, result: AnalysisResult) -> bool:
    """The welding deformation input and its effect: keel shape (components, total,
    survey points), joint angles and zone curvatures, and the resulting block gaps."""
    from .deformation import detrended_lift, lift_components

    fig.clear()
    _layout(fig)
    defo = result.deformation
    if defo is None or result.weld_gap is None:
        return False
    ship = result.case.ship
    x = result.springs.x
    labels = result.strength[result.design_method].blocks.label
    given = result.case.deformation
    survey = given.measured if given is not None else None
    ident = result.identification
    xc = np.linspace(ship.x_aft, ship.x_fwd, 600)
    axes = fig.subplots(3, 1, sharex=True, gridspec_kw={"height_ratios": [1.5, 0.9, 1.0]})

    # 1 — keel shape after welding
    ax = axes[0]
    total = detrended_lift(defo, xc, x)
    colours = {"joint kinks": "#e8a33d", "zone shrinkage": "#7b5ea7",
               "survey (interpolated)": "#2f6b3a"}
    parts = lift_components(defo, xc, x)
    if len(parts) > 1:
        for name, curve in parts.items():
            ax.plot(xc, curve / MM, ls="--", lw=1.4, color=colours.get(name, "0.5"),
                    label=f"{name} (straight line removed)")
    ax.plot(xc, total / MM, color="#1f4e79", lw=2.4, label="keel lift used for the gaps")
    if survey is not None:
        model = detrended_lift(defo, survey.x, x)
        c1, c0 = np.polyfit(survey.x, survey.lift - model, 1)
        ax.plot(survey.x, (survey.lift - (c0 + c1 * survey.x)) / MM, "s", ms=5, mfc="white",
                color="#2f6b3a", label=f"survey points ({survey.x.size})", zorder=5)
    ax.plot(x, result.weld_gap / MM, "o", ms=3.5, color="#d64545", label="block positions",
            zorder=6)
    for k, j in enumerate(defo.joints, start=1):
        ax.axvline(j.x, color="0.65", ls=":", lw=1)
        ax.text(j.x, 1.0, f"J{k}", transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=8, color="0.35")
    if ident is not None:
        source = (f"joint angles identified from the survey (rms residual "
                  f"{ident.rms_residual / MM:.2f} mm)")
    elif survey is not None and not defo.joints and not defo.zones:
        source = "surveyed keel line used directly"
    else:
        source = "welding model (joint kinks and shrinkage zones)"
    ax.set_title(f"{ship.name} — welding deformation data: {source}", fontsize=11, pad=16)
    ax.set_ylabel("keel lift [mm]")
    ax.legend(fontsize=8, loc="upper center", ncol=3)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.25)

    # 2 — input data: joint angles and zone curvature
    ax = axes[1]
    if defo.joints or defo.zones:
        jx = np.array([j.x for j in defo.joints])
        width = 0.03 * (ship.x_fwd - ship.x_aft)
        if ident is not None and given is not None and any(j.angle for j in given.joints):
            ax.bar(jx - width / 2, [j.angle * 1e3 for j in given.joints], width, color="#b5b5b5",
                   label="given angle")
            ax.bar(jx + width / 2, [j.angle * 1e3 for j in defo.joints], width, color="#e8a33d",
                   label="identified angle")
        elif defo.joints:
            ax.bar(jx, [j.angle * 1e3 for j in defo.joints], width, color="#e8a33d",
                   label="identified joint angle" if ident is not None else "joint angle")
        for xj, j in zip(jx, defo.joints):
            ax.text(xj, j.angle * 1e3, f"{j.angle * 1e3:.2f}", ha="center", va="bottom",
                    fontsize=8)
        ax.set_ylabel("joint angle [mrad]")
        if defo.joints:
            ax.set_ylim(0, max(max(j.angle for j in defo.joints) * 1e3 * 1.3, 0.1))
        if defo.zones:
            ax2 = ax.twinx()
            ax2.spines["right"].set_visible(True)
            for z in defo.zones:
                ax2.plot([z.xa, z.xb], [z.curvature * 1e6] * 2, color="#7b5ea7", lw=3,
                         solid_capstyle="butt")
            ax2.set_ylabel("shrinkage curvature [10⁻⁶ /m]", color="#7b5ea7")
            top = max(abs(z.curvature) for z in defo.zones) * 1e6
            ax2.set_ylim(0, max(top * 1.6, 1e-3))
            ax.plot([], [], color="#7b5ea7", lw=3, label="zone curvature (right axis)")
        ax.legend(fontsize=8, loc="upper center", ncol=3)
    else:
        ax.text(0.5, 0.5, "surveyed keel line used directly — no joint or zone parameters",
                transform=ax.transAxes, ha="center", va="center", color="0.4")
        ax.set_yticks([])

    # 3 — resulting gaps at the blocks
    ax = axes[2]
    design = result.supports[result.design_method]
    open_ = design.reactions <= 0
    spacing = np.min(np.diff(np.unique(x))) if np.unique(x).size > 1 else 1.0
    ax.bar(x, result.weld_gap / MM, 0.7 * spacing,
           color=np.where(open_, "#d64545", "#9fb7d0"), edgecolor="#333333", lw=0.4)
    for i in np.flatnonzero(open_):
        ax.text(x[i], result.weld_gap[i] / MM, labels[i], rotation=90, ha="center", va="bottom",
                fontsize=7, color="#d64545")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#9fb7d0", label="block in contact (as built)"),
                       Patch(color="#d64545", label="block without contact (as built)")],
              fontsize=8, loc="upper center", ncol=2)
    ax.set_ylabel("gap from welding [mm]")
    ax.set_xlabel("x from AP [m]")
    ax.set_ylim(0, max(result.weld_gap.max() / MM * 1.3, 1.0))
    for a in axes:
        a.grid(True, alpha=0.3)
    return True


def plot_spacing(result: AnalysisResult, path: Path) -> Path | None:
    if result.spacing is None:
        return None
    from matplotlib.figure import Figure

    fig = Figure(figsize=(13, 14))
    draw_spacing(fig, result)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    return path


def draw_spacing(fig, result: AnalysisResult) -> bool:
    """Designed spacing against the given spacing, reactions against strength, convergence."""
    fig.clear()
    _layout(fig)
    sp = result.spacing
    if sp is None:
        return False
    ship = result.case.ship
    target = sp.target
    x0 = np.array([b.x for b in sp.initial_blocks])
    x1 = np.array([b.x for b in sp.final_blocks])
    axes = fig.subplots(4, 1, gridspec_kw={"height_ratios": [1.0, 1.1, 1.3, 0.8]})

    # 0 — weight distribution with the given and the designed blocks
    draw_weight_panel(axes[0], result, x1, x0, "given blocks", "designed blocks")
    axes[0].set_title(f"{ship.name} — block spacing design: {sp.status.replace('_', ' ')}\n"
                      f"{len(x1)} blocks (strengths kept), max reaction/strength "
                      f"{np.max(sp.initial_utilisation):.3f} → {np.max(sp.final_utilisation):.3f}",
                      fontsize=11)
    axes = axes[1:]

    # 1 — spacing along the block line
    ax = axes[0]
    for xs, colour, label, z in ((x0, "#9fb7d0", "given spacing", 1),
                                 (x1, "#1f4e79", "designed spacing", 2)):
        gaps = np.diff(xs)
        ax.bar(xs[:-1], gaps, gaps, align="edge", color=colour, alpha=0.55 if z == 1 else 0.85,
               edgecolor="white", lw=0.6, label=label, zorder=z)
    gaps1 = sp.spacing
    start = 0
    for i in range(1, gaps1.size + 1):                     # label runs of equal spacing that
        if i == gaps1.size or abs(gaps1[i] - gaps1[start]) > 1e-6:   # are long enough
            if x1[i] - x1[start] >= 0.04 * (ship.x_fwd - ship.x_aft):
                ax.text(0.5 * (x1[start] + x1[i]), gaps1[start], f"{gaps1[start]:.1f} m",
                        ha="center", va="bottom", fontsize=7.5, color="#1f4e79")
            start = i
    ax.plot(x1, np.zeros_like(x1), "^", ms=5, color="#1f4e79", clip_on=False, zorder=3)
    ax.set_xlim(ship.x_aft, ship.x_fwd)
    ax.set_ylabel("spacing to next block [m]")
    ax.set_title("spacing to the next block", fontsize=9.5, loc="left")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    # 2 — reactions against strength
    ax = axes[1]
    strength = np.array([b.effective_strength for b in sp.final_blocks])
    u = sp.final_utilisation
    gaps = np.diff(x1)
    local = np.minimum(np.concatenate([[gaps[0]], gaps]), np.concatenate([gaps, [gaps[-1]]]))
    colours = np.where(sp.final_evaluation.reactions <= 0, _OPEN,
                       np.where(u > target, _OVER, np.where(u > 0.9 * target, _NEAR, _OK)))
    ax.bar(x1, to_tonnes(sp.final_evaluation.reactions), 0.7 * local, color=colours,
           edgecolor="#333333", lw=0.4, label="designed spacing")
    ax.plot(x0, to_tonnes(sp.initial_evaluation.reactions), "D", ms=3.5, color="#222222",
            label="given spacing")
    ax.step(x1, to_tonnes(strength * target), where="mid", color=_OVER, ls="--", lw=1.2,
            label="block strength" + (f" × {target:g}" if target != 1 else ""))
    ax.set_xlim(ship.x_aft, ship.x_fwd)
    ax.set_ylabel("block reaction [t]")
    ax.set_xlabel("x from AP [m]")
    ax.legend(fontsize=8, loc="upper center", ncol=3)
    ax.grid(True, alpha=0.3)

    # 3 — convergence over the analyses
    ax = axes[2]
    styles = {"start": ("s", "#222222"), "equidistribution": ("o", "#4c9a5f"),
              "slp": ("o", "#1f4e79"), "grid": ("*", "#e8a33d")}
    shown = set()
    for h in sp.history:
        marker, colour = styles.get(h.phase, ("o", "grey"))
        ax.plot(h.evaluations, h.max_utilisation, marker, color=colour, ms=7 if marker == "*" else 5,
                label=None if h.phase in shown else h.phase)
        shown.add(h.phase)
    ax.plot([h.evaluations for h in sp.history], [h.max_utilisation for h in sp.history],
            color="0.6", lw=0.8, zorder=0)
    ax.axhline(target, color=_OVER, ls="--", lw=1, label="target")
    ax.set_xlabel("number of analyses")
    ax.set_ylabel("max reaction / strength [–]")
    ax.legend(fontsize=8, loc="upper right", ncol=5)
    ax.grid(True, alpha=0.3)
    return True


