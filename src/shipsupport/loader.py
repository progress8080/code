"""Read a case from YAML (or an equivalent dict) and convert to SI.

Input conventions: lengths in m, weights/forces in tonnes (`_t`), moments in
t·m (`_tm`), distributed weight in t/m (`_t_per_m`), stresses and moduli as named
by the key suffix (`_MPa`, `_GPa`), gaps in mm. See examples/example_case.yaml.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from .model import (BendingEnvelope, Block, BlockType, Case, InputError, Layer, MeasuredProfile,
                    Section, Settings, ShearEnvelope, Ship, ShrinkageZone, SpacingDesign,
                    WeightCurve, WeightItem, WeldingDeformation, WeldJoint)
from .units import GPA, KN_PER_MM, MM, MPA, TONNE_FORCE

_MISSING = object()


def load_case(path: str | Path) -> Case:
    path = Path(path)
    with path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, Mapping):
        raise InputError(f"{path}: the top level must be a mapping")
    return case_from_dict(inline_data_files(dict(data), path.parent))


def read_survey_csv(path: str | Path) -> tuple[list[float], list[float]]:
    """Read a keel survey: CSV with a header row and the columns x_m (or x) and lift_mm
    (or lift); lines starting with '#' are comments."""
    import csv

    path = Path(path)
    try:
        lines = [ln for ln in path.read_text(encoding="utf-8-sig").splitlines()
                 if ln.strip() and not ln.lstrip().startswith("#")]
    except OSError as exc:
        raise InputError(f"Cannot read survey file {path}: {exc}") from None
    rows = list(csv.reader(lines))
    if len(rows) < 3:
        raise InputError(f"Survey file {path}: need a header and at least two points")
    header = [h.strip().lower() for h in rows[0]]
    try:
        ix = next(i for i, h in enumerate(header) if h in ("x", "x_m"))
        iz = next(i for i, h in enumerate(header) if h in ("lift_mm", "lift"))
    except StopIteration:
        raise InputError(f"Survey file {path}: header must contain 'x_m' and 'lift_mm'") from None
    xs, zs = [], []
    for n, row in enumerate(rows[1:], start=2):
        try:
            xs.append(float(row[ix]))
            zs.append(float(row[iz]))
        except (ValueError, IndexError):
            raise InputError(f"Survey file {path}, line {n}: expected numbers") from None
    return xs, zs


def inline_data_files(data: dict, base_dir: str | Path) -> dict:
    """Replace the survey file reference welding_deformation.measured.file (path relative
    to base_dir) by its contents, the x / lift_mm lists."""
    welding = data.get("welding_deformation")
    measured = (welding or {}).get("measured") or {}
    if measured.get("file"):
        path = Path(measured["file"])
        path = path if path.is_absolute() else Path(base_dir) / path
        xs, zs = read_survey_csv(path)
        measured = {k: v for k, v in measured.items() if k != "file"}
        measured.update(x=xs, lift_mm=zs, source=str(measured.get("source", path.name)))
        data["welding_deformation"] = {**welding, "measured": measured}
    return data


def _get(d: Mapping, key: str, ctx: str, default: Any = _MISSING) -> Any:
    if not isinstance(d, Mapping):
        raise InputError(f"{ctx}: expected a mapping")
    if d.get(key) is not None:
        return d[key]
    if default is _MISSING:
        raise InputError(f"{ctx}: missing required key '{key}'")
    return default


def _float(d: Mapping, key: str, ctx: str, default: Any = _MISSING) -> float | None:
    value = _get(d, key, ctx, default)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise InputError(f"{ctx}.{key}: expected a number, got {value!r}") from None


def _scaled(d: Mapping, key: str, ctx: str, scale: float, default: Any = None) -> float | None:
    value = _float(d, key, ctx, default)
    return None if value is None else value * scale


def _column(d: Mapping, key: str, ctx: str, scale: float = 1.0) -> np.ndarray:
    try:
        return np.asarray(_get(d, key, ctx), dtype=float) * scale
    except (TypeError, ValueError):
        raise InputError(f"{ctx}.{key}: expected a list of numbers") from None


def structure_positions(spec: Mapping) -> np.ndarray:
    """Structure positions from an explicit list `x: [..]`, or from `start` plus exactly
    two of `end`, `spacing`, `count`:

        start, end, spacing   every `spacing` from start up to end (inclusive if on grid)
        start, spacing, count `count` positions at `spacing`
        start, end, count     `count` positions evenly spread over [start, end]
    """
    ctx = "structure"
    if spec.get("x") is not None:
        try:
            return np.atleast_1d(np.asarray(spec["x"], dtype=float))
        except (TypeError, ValueError):
            raise InputError(f"{ctx}.x: expected a list of numbers") from None
    start = _float(spec, "start", ctx)
    end = _float(spec, "end", ctx, None)
    spacing = _float(spec, "spacing", ctx, None)
    count = _float(spec, "count", ctx, None)
    given = [k for k, v in (("end", end), ("spacing", spacing), ("count", count)) if v is not None]
    if len(given) != 2:
        raise InputError(f"{ctx}: give start plus exactly two of end, spacing, count "
                         f"(got {', '.join(given) or 'none'})")
    if count is not None:
        if count != int(count) or count < 1:
            raise InputError(f"{ctx}.count: must be a whole number >= 1 (got {count:g})")
        count = int(count)
    if spacing is not None and spacing <= 0:
        raise InputError(f"{ctx}.spacing: must be > 0")
    if end is not None and end < start:
        raise InputError(f"{ctx}: end must be >= start")

    if count is None:
        n = int(np.floor((end - start) / spacing + 1e-9))
        return start + spacing * np.arange(n + 1)
    if spacing is not None:
        return start + spacing * np.arange(count)
    if count == 1:
        if end != start:
            raise InputError(f"{ctx}: a single position needs end == start")
        return np.array([start])
    return np.linspace(start, end, count)


def _weight(d: Mapping, i: int) -> WeightItem:
    """A weight item with its LCG given as `lcg` or as the moment about AP `moment_tm`
    (LCG = moment / weight); `x` (or no extent) makes it a point load."""
    ctx = f"weights[{i}]"
    name = str(d.get("name", f"item {i + 1}"))
    weight_t = _float(d, "weight_t", ctx)
    if weight_t == 0:
        raise InputError(f"{ctx} ('{name}'): weight_t must not be zero")
    candidates = []
    for key in ("x", "lcg"):
        if d.get(key) is not None:
            candidates.append((key, _float(d, key, ctx)))
    if d.get("moment_tm") is not None:
        candidates.append(("moment_tm / weight_t", _float(d, "moment_tm", ctx) / weight_t))
    if not candidates:
        raise InputError(f"{ctx} ('{name}'): give lcg, moment_tm or (for a point load) x")
    (key0, lcg), rest = candidates[0], candidates[1:]
    for key, value in rest:
        if not math.isclose(value, lcg, rel_tol=1e-6, abs_tol=1e-6):
            raise InputError(f"{ctx} ('{name}'): {key} = {value:.6g} m does not agree with "
                             f"{key0} = {lcg:.6g} m; give one of them or consistent values")
    weight = weight_t * TONNE_FORCE
    has_a, has_b = d.get("xa") is not None, d.get("xb") is not None
    if d.get("x") is not None or not (has_a or has_b):
        if has_a or has_b:
            raise InputError(f"{ctx} ('{name}'): a point load (x) has no extent xa / xb")
        return WeightItem(name, weight, lcg, lcg, lcg)
    if has_a != has_b:
        raise InputError(f"{ctx} ('{name}'): give both xa and xb, or neither for a point load")
    return WeightItem(name, weight, lcg, _float(d, "xa", ctx), _float(d, "xb", ctx))


def _weight_curve(spec) -> WeightCurve:
    """weight_curve: {name, points: [{x, q_t_per_m}, ...]} or a plain list of points."""
    ctx = "weight_curve"
    if isinstance(spec, Mapping):
        name = str(spec.get("name") or "weight curve")
        points = _get(spec, "points", ctx)
    else:
        name, points = "weight curve", spec
    if not isinstance(points, (list, tuple)):
        raise InputError(f"{ctx}.points: expected a list of {{x, q_t_per_m}} points")
    xs, qs = [], []
    for k, p in enumerate(points):
        where = f"{ctx}.points[{k}]"
        xs.append(_float(p, "x", where))
        qs.append(_float(p, "q_t_per_m", where) * TONNE_FORCE)
    return WeightCurve(name, xs, qs)


def _section(d: Mapping, i: int) -> Section:
    ctx = f"sections[{i}]"
    if "EI_Nm2" in d:
        EI = _float(d, "EI_Nm2", ctx)
    else:
        EI = _float(d, "E_GPa", ctx, 206.0) * GPA * _float(d, "I_m4", ctx)
    shear_area = _float(d, "shear_area_m2", ctx, None)
    kGA = _float(d, "kGA_N", ctx, None)
    if kGA is None and shear_area is not None:
        kGA = _float(d, "G_GPa", ctx, 79.3) * GPA * shear_area
    return Section(
        xa=_float(d, "xa", ctx), xb=_float(d, "xb", ctx), EI=EI, kGA=kGA,
        z_deck=_float(d, "z_deck_m3", ctx, None), z_keel=_float(d, "z_keel_m3", ctx, None),
        shear_area=shear_area,
    )


def _block_type(name: str, d: Mapping) -> BlockType:
    ctx = f"blocks.types.{name}"
    layers = tuple(
        Layer(_float(layer, "thickness", f"{ctx}.layers[{j}]"),
              _float(layer, "E_MPa", f"{ctx}.layers[{j}]") * MPA)
        for j, layer in enumerate(d.get("layers") or []))
    return BlockType(
        name=name,
        contact_width=_float(d, "contact_width", ctx),
        contact_length=_float(d, "contact_length", ctx),
        layers=layers,
        stiffness_override=_scaled(d, "stiffness_kN_per_mm", ctx, KN_PER_MM),
        capacity=_scaled(d, "capacity_t", ctx, TONNE_FORCE),
        allowable_pressure=_scaled(d, "allowable_pressure_MPa", ctx, MPA),
        crush_stress=_scaled(d, "crush_stress_MPa", ctx, MPA),
        post_crush_ratio=_float(d, "post_crush_ratio", ctx, 0.05),
    )


def _broadcast(value, n: int, what: str) -> list[float | None]:
    """A scalar for all entries or a list of exactly n values (empty entries = None)."""
    if value is None:
        return [None] * n
    if isinstance(value, (list, tuple)) and len(value) != n:
        raise InputError(f"blocks.plan.{what}: {len(value)} values given, {n} expected")
    try:
        if isinstance(value, (list, tuple)):
            return [None if v is None or v == "" else float(v) for v in value]
        return [float(value)] * n
    except (TypeError, ValueError):
        raise InputError(f"blocks.plan.{what}: expected a number or a list of numbers") from None


def _blocks(d: Mapping) -> list[Block]:
    """The block plan: count, first position, and the spacing, strength and gap of every block.

    blocks.plan: {type, start_x, count, spacing: s | [s_1 … s_{N-1}],
                  strength_t: S | [S_1 … S_N]   (optional, default: the type's capacity),
                  gap_mm: g | [g_1 … g_N]        (optional initial gap, moves with the block)}

    The blocks are labelled B01, B02, … from aft to forward.
    """
    types = {str(name): _block_type(str(name), spec)
             for name, spec in _get(d, "types", "blocks").items()}
    ctx = "blocks.plan"
    if d.get("plan") is None and d.get("groups") is not None:
        raise InputError("blocks.groups is not supported; define the blocks with blocks.plan "
                         "{type, start_x, count, spacing, strength_t, gap_mm}")
    plan = _get(d, "plan", "blocks")
    type_name = str(_get(plan, "type", ctx))
    if type_name not in types:
        raise InputError(f"{ctx}: unknown block type '{type_name}'")
    spacing, strength, gap = plan.get("spacing"), plan.get("strength_t"), plan.get("gap_mm")
    count = plan.get("count")
    if count is None:
        if isinstance(strength, (list, tuple)):
            count = len(strength)
        elif isinstance(gap, (list, tuple)):
            count = len(gap)
        elif isinstance(spacing, (list, tuple)):
            count = len(spacing) + 1
        else:
            raise InputError(f"{ctx}: give count, or the spacing / strength of every block")
    try:
        whole = float(count) == int(float(count))
    except (TypeError, ValueError):
        whole = False
    if not whole or int(float(count)) < 2:
        raise InputError(f"{ctx}.count must be a whole number >= 2 (got {count!r})")
    n = int(float(count))
    gaps = _broadcast(_get(plan, "spacing", ctx), n - 1, "spacing")
    if any(g is None or g <= 0 for g in gaps):
        raise InputError(f"{ctx}.spacing: every spacing must be given and > 0")
    xs = _float(plan, "start_x", ctx) + np.concatenate([[0.0], np.cumsum(gaps)])
    strengths = _broadcast(strength, n, "strength_t")
    initial_gaps = _broadcast(gap, n, "gap_mm")
    return [Block(round(float(x), 6), types[type_name],
                  0.0 if g is None else g * MM, f"B{i + 1:02d}",
                  None if s is None else s * TONNE_FORCE)
            for i, (x, s, g) in enumerate(zip(xs, strengths, initial_gaps, strict=True))]


def _spacing(d: Mapping) -> SpacingDesign:
    ctx = "spacing_design"
    return SpacingDesign(
        method=str(d.get("method", "fem")),
        target_utilisation=_float(d, "target_utilisation", ctx, 1.0),
        min_spacing=_float(d, "min_spacing", ctx, 1.0),
        max_spacing=_float(d, "max_spacing", ctx, None),
        x_min=_float(d, "x_min", ctx, None),
        x_max=_float(d, "x_max", ctx, None),
        fix_end_blocks=bool(d.get("fix_end_blocks", False)),
        snap_to_structure=bool(d.get("snap_to_structure", True)),
        grid=_float(d, "grid", ctx, None),
        max_iterations=int(_float(d, "max_iterations", ctx, 20)),
    )


def _ship(data: Mapping) -> Ship:
    d = _get(data, "ship", "root")
    lpp = _float(d, "lpp", "ship")

    structure = data.get("structure")
    permissible = data.get("permissible") or {}
    bm = sf = None
    if permissible.get("bm"):
        rows = permissible["bm"]
        bm = BendingEnvelope(
            x=[_float(r, "x", "permissible.bm") for r in rows],
            sag=[_float(r, "sag_tm", "permissible.bm") * TONNE_FORCE for r in rows],
            hog=[_float(r, "hog_tm", "permissible.bm") * TONNE_FORCE for r in rows])
    if permissible.get("sf"):
        rows = permissible["sf"]
        sf = ShearEnvelope(
            x=[_float(r, "x", "permissible.sf") for r in rows],
            limit=[_float(r, "sf_t", "permissible.sf") * TONNE_FORCE for r in rows])

    weights = [_weight(w, i) for i, w in enumerate(data.get("weights") or [])]
    curve = _weight_curve(data["weight_curve"]) if data.get("weight_curve") else None
    if curve is not None:
        weights += curve.items()
    if not weights:
        raise InputError("root: give weight items (weights) and/or a weight curve (weight_curve)")

    return Ship(
        name=str(d.get("name", "Unnamed ship")),
        lpp=lpp,
        weights=weights,
        sections=[_section(s, i) for i, s in enumerate(_get(data, "sections", "root"))],
        x_aft=_float(d, "x_aft", "ship", 0.0),
        x_fwd=_float(d, "x_fwd", "ship", lpp),
        declared_weight=_scaled(d, "weight_t", "ship", TONNE_FORCE),
        declared_lcg=_float(d, "lcg", "ship", None),
        structure_positions=structure_positions(structure) if structure else None,
        permissible_bm=bm,
        permissible_sf=sf,
        weight_curve=curve,
    )


def _settings(d: Mapping) -> Settings:
    ctx = "settings"
    methods = d.get("methods", ["rigid", "fem"])
    return Settings(
        methods=tuple(str(m) for m in methods),
        max_element_length=_float(d, "max_element_length", ctx, 0.5),
        shear_deformation=bool(d.get("shear_deformation", True)),
        allowable_bending_stress=_scaled(d, "allowable_bending_MPa", ctx, MPA),
        allowable_shear_stress=_scaled(d, "allowable_shear_MPa", ctx, MPA),
        structure_tolerance=_float(d, "structure_tolerance", ctx, 0.15),
        weight_tolerance=_float(d, "weight_tolerance", ctx, 0.005),
        max_iterations=int(_float(d, "max_iterations", ctx, 200)),
    )


def _welding(d: Mapping) -> WeldingDeformation:
    """welding_deformation: joints (angle_mrad), zones (curvature or deck/bottom shrinkage),
    measured survey (lift_mm)."""
    # the angle may be omitted when it is identified from a survey
    joints = tuple(
        WeldJoint(_float(j, "x", f"welding_deformation.joints[{i}]"),
                  _float(j, "angle_mrad", f"welding_deformation.joints[{i}]", 0.0) * 1e-3)
        for i, j in enumerate(d.get("joints") or []))
    zones = []
    for i, z in enumerate(d.get("zones") or []):
        ctx = f"welding_deformation.zones[{i}]"
        if z.get("curvature_mm_per_m2") is not None:
            kappa = _float(z, "curvature_mm_per_m2", ctx) * 1e-3
        else:
            deck = _float(z, "deck_shrinkage_mm_per_m", ctx) * 1e-3
            bottom = _float(z, "bottom_shrinkage_mm_per_m", ctx, 0.0) * 1e-3
            depth = _float(z, "depth", ctx)
            if depth <= 0:
                raise InputError(f"{ctx}.depth: must be positive")
            kappa = (deck - bottom) / depth
        zones.append(ShrinkageZone(_float(z, "xa", ctx), _float(z, "xb", ctx), kappa))
    measured = None
    if d.get("measured"):
        m = d["measured"]
        ctx = "welding_deformation.measured"
        if m.get("file"):
            raise InputError(f"{ctx}.file: survey files are read by load_case(); for a dict "
                             f"input call inline_data_files(data, base_dir) first")
        measured = MeasuredProfile(
            x=_column(m, "x", ctx), lift=_column(m, "lift_mm", ctx, MM),
            interpolation=str(m.get("interpolation", "pchip")),
            identify_joints=bool(m.get("identify_joint_angles", False)),
            regularization=_float(m, "regularization", ctx, 0.0))
    return WeldingDeformation(joints=joints, zones=tuple(zones), measured=measured)


def case_from_dict(data: Mapping) -> Case:
    spacing = data.get("spacing_design")
    return Case(
        ship=_ship(data),
        blocks=_blocks(_get(data, "blocks", "root")),
        settings=_settings(data.get("settings") or {}),
        deformation=_welding(data["welding_deformation"]) if data.get("welding_deformation") else None,
        spacing=_spacing(spacing) if spacing and spacing.get("enabled", True) else None,
    )
