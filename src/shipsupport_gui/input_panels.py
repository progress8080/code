"""Input editor: tabs that edit a case dict with the same schema as the YAML file."""

from __future__ import annotations

import numpy as np
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (QCheckBox, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QScrollArea,
                             QSplitter, QTabWidget, QVBoxLayout, QWidget)

from shipsupport.loader import structure_positions
from shipsupport.model import InputError

from .plan_panel import BlockPlanEditor
from .widgets import Column, Field, FormEditor, RecordTable

# Sample welding deformation (placeholder values, same as examples/welding_case.yaml):
# a hull built from six 20 m construction blocks, kinks at the five erection joints and
# deck/bottom shrinkage in every block.
SAMPLE_WELDING = {
    "joints": [{"x": 20.0, "angle_mrad": 0.9}, {"x": 40.0, "angle_mrad": 0.5},
               {"x": 60.0, "angle_mrad": 0.4}, {"x": 80.0, "angle_mrad": 0.5},
               {"x": 100.0, "angle_mrad": 0.6}],
    "zones": [{"xa": float(a), "xb": float(a + 20), "deck_shrinkage_mm_per_m": 0.24,
               "bottom_shrinkage_mm_per_m": 0.20, "depth": 10.0} for a in range(0, 120, 20)],
}

# Sample weight distribution (placeholder values, same as examples/weight_curve_case.yaml) for
# a hull from x = 0 to 120 m: the lightship as a piecewise-linear curve with jumps (engine
# room, bulkheads, anchors) and three deadweight items given by their moments about AP.
SAMPLE_WEIGHT_CURVE = {
    "name": "lightship curve",
    "points": [(0.0, 8.0), (4.0, 16.0), (4.0, 30.0), (14.0, 36.0), (16.0, 36.0), (16.0, 64.0),
               (28.0, 64.0), (28.0, 46.0), (38.0, 42.0), (38.0, 34.0), (60.0, 36.0),
               (80.0, 35.0), (100.0, 30.0), (104.0, 28.0), (104.0, 36.0), (112.0, 33.0),
               (118.0, 16.0), (120.0, 0.0)],
    "items": [{"name": "residual fuel", "weight_t": 200.0, "moment_tm": 14000.0,
               "xa": 60.0, "xb": 80.0},
              {"name": "fresh water", "weight_t": 60.0, "moment_tm": 6000.0,
               "xa": 94.0, "xb": 106.0},
              {"name": "stores", "weight_t": 40.0, "moment_tm": 1200.0}],
}

DEFAULT_SETTINGS = {
    "methods": ["rigid", "fem"],
    "shear_deformation": True,
}


def parse_layers(text: str, ctx: str) -> list[dict]:
    """'0.9@30000; 0.2@1000' -> [{'thickness': 0.9, 'E_MPa': 30000}, ...]"""
    layers = []
    for raw in text.split(";"):
        part = raw.strip()
        if not part:
            continue
        try:
            thickness, modulus = (float(v) for v in part.split("@"))
        except ValueError:
            raise ValueError(f"{ctx}: layer '{part}' must be written thickness@E_MPa, "
                             f"e.g. 0.9@30000") from None
        layers.append({"thickness": thickness, "E_MPa": modulus})
    return layers


def format_layers(layers) -> str:
    return "; ".join(f"{float(layer['thickness']):.12g}@{float(layer['E_MPa']):.12g}" for layer in layers or [])


def curve_totals(points: list[dict]) -> tuple[float, float]:
    """Weight [t] and moment about AP [t·m] of a piecewise-linear weight curve
    (points with x and q_t_per_m), exact for every trapezoidal segment."""
    weight = moment = 0.0
    for a, b in zip(points[:-1], points[1:]):
        x1, x2, q1, q2 = a["x"], b["x"], a["q_t_per_m"], b["q_t_per_m"]
        length = x2 - x1
        weight += 0.5 * (q1 + q2) * length
        moment += length * (q1 * (2 * x1 + x2) + q2 * (x1 + 2 * x2)) / 6.0
    return weight, moment


def _scroll(widget: QWidget) -> QScrollArea:
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setWidget(widget)
    return area


def _column_widget(*widgets, stretch_last=False) -> QWidget:
    w = QWidget()
    layout = QVBoxLayout(w)
    for item in widgets:
        layout.addWidget(item)
    if not stretch_last:
        layout.addStretch()
    return w


class CaseEditor(QTabWidget):
    changed = pyqtSignal()
    sample_weight_curve_loaded = pyqtSignal()     # the sample replaced the weights: run it

    def __init__(self, parent=None):
        super().__init__(parent)
        # spacing_design of the loaded case: not edited in the GUI, written back unchanged
        self.spacing_design: dict | None = None

        # Ship ---------------------------------------------------------------
        self.ship_form = FormEditor("Ship", [
            Field("name", "Name", "str"),
            Field("lpp", "Lpp [m]"),
            Field("x_aft", "Hull aft end x [m]", default=0.0,
                  tooltip="x of the aft end of the hull (negative if aft of AP)"),
            Field("x_fwd", "Hull forward end x [m]", default="Lpp",
                  tooltip="x of the forward end of the hull; empty = Lpp"),
            Field("weight_t", "Declared weight [t]", tooltip="From the docking plan; cross-check only"),
            Field("lcg", "Declared LCG [m]", tooltip="From the docking plan; cross-check only"),
        ])
        note = QLabel("x is measured from the aft perpendicular (AP), positive forward.")
        note.setObjectName("hint")
        self.addTab(_scroll(_column_widget(self.ship_form, note)), "Ship")

        # Weights ------------------------------------------------------------
        self.weights = RecordTable("Weight items", [
            Column("name", "Name", "str", width=170),
            Column("weight_t", "Weight [t]"),
            Column("lcg", "LCG [m]", tooltip="LCG from AP; or give the moment"),
            Column("moment_tm", "Moment [t·m]",
                   tooltip="longitudinal moment about AP; LCG = moment / weight"),
            Column("xa", "xa [m]", tooltip="aft end of extent; leave xa and xb empty for a point load"),
            Column("xb", "xb [m]", tooltip="forward end of extent"),
        ], hint="Each item is spread as a trapezoid over [xa, xb] preserving weight and LCG "
                "(LCG must be in the middle third). Give the LCG or the moment about AP (both "
                "must agree if both are given). Leave xa/xb empty for a point load. "
                "Ctrl+V pastes from Excel.")
        self.curve_name = QLineEdit()
        self.curve_name.setPlaceholderText("weight curve")
        self.curve = RecordTable("Weight curve (piecewise linear, zero outside)", [
            Column("x", "x [m]"),
            Column("q_t_per_m", "q [t/m]", tooltip="weight per metre, >= 0"),
        ], hint="Points from aft to forward; the same x twice is a jump. Every segment is an "
                "exact trapezoid. Use it for the lightship distribution, alone or with items.")
        self.curve.add_button("Load sample curve", self.ask_load_sample_weight_curve)
        curve_page = QWidget()
        curve_layout = QVBoxLayout(curve_page)
        curve_layout.setContentsMargins(0, 0, 0, 0)
        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Curve name"))
        name_row.addWidget(self.curve_name)
        curve_layout.addLayout(name_row)
        curve_layout.addWidget(self.curve)
        weights_split = QSplitter(Qt.Vertical)
        weights_split.addWidget(self.weights)
        weights_split.addWidget(curve_page)
        weights_split.setSizes([400, 300])
        self.weight_summary = QLabel()
        self.weight_summary.setTextFormat(Qt.RichText)
        self.weights.changed.connect(self._update_weight_summary)
        self.curve.changed.connect(self._update_weight_summary)
        self.addTab(_column_widget(weights_split, self.weight_summary, stretch_last=True), "Weights")

        # Hull sections --------------------------------------------------------
        self.sections = RecordTable("Hull girder sections", [
            Column("xa", "xa [m]"),
            Column("xb", "xb [m]"),
            Column("I_m4", "I [m⁴]"),
            Column("E_GPa", "E [GPa]", tooltip="default 206"),
            Column("shear_area_m2", "As [m²]", tooltip="effective shear area; enables Timoshenko and τ"),
            Column("G_GPa", "G [GPa]", tooltip="default 79.3"),
            Column("z_deck_m3", "Z deck [m³]"),
            Column("z_keel_m3", "Z keel [m³]"),
            Column("EI_Nm2", "EI [N·m²]", tooltip="optional: overrides E·I"),
            Column("kGA_N", "κGA [N]", tooltip="optional: overrides G·As for shear stiffness"),
        ], hint="Sections must cover the whole hull extent.")
        self.addTab(_column_widget(self.sections, stretch_last=True), "Hull sections")

        # Blocks: types and the block plan ---------------------------------------
        self.block_types = RecordTable("Block types", [
            Column("name", "Name", "str", width=90),
            Column("contact_width", "Width [m]"),
            Column("contact_length", "Length [m]"),
            Column("layers", "Layers t@E_MPa; …", "str", width=220,
                   tooltip="Layers in series, e.g. 0.9@30000; 0.2@1000; 0.1@400"),
            Column("stiffness_kN_per_mm", "k override [kN/mm]"),
            Column("capacity_t", "Capacity [t]",
                   tooltip="default strength of a block of this type"),
            Column("allowable_pressure_MPa", "Allow. p [MPa]"),
            Column("crush_stress_MPa", "Crush σ [MPa]"),
            Column("post_crush_ratio", "Post-crush k ratio", tooltip="default 0.05"),
        ])
        self.block_types.changed.connect(self._update_type_names)
        self.plan_editor = BlockPlanEditor()
        blocks_split = QSplitter(Qt.Vertical)
        blocks_split.addWidget(self.block_types)
        blocks_split.addWidget(self.plan_editor)
        blocks_split.setSizes([140, 760])
        blocks_split.setStretchFactor(1, 1)
        self.blocks_tab = blocks_split
        self.addTab(self.blocks_tab, "Blocks")

        # Structure & limits ----------------------------------------------------
        self.structure_box = QGroupBox("Transverse structure under blocks (floors, bulkheads)")
        self.structure_box.setCheckable(True)
        self.structure_form = FormEditor("Structure", [
            Field("start", "Start x [m]"),
            Field("end", "End x [m]"),
            Field("spacing", "Spacing [m]"),
            Field("count", "Count", "int"),
            Field("x", "Explicit x list [m]", "floats"),
        ])
        QVBoxLayout(self.structure_box).addWidget(self.structure_form)
        self.bm_table = RecordTable("Permissible bending moment", [
            Column("x", "x [m]"), Column("sag_tm", "Sag [t·m]"), Column("hog_tm", "Hog [t·m]")])
        self.sf_table = RecordTable("Permissible shear force", [
            Column("x", "x [m]"), Column("sf_t", "SF [t]")])
        limits = QSplitter(Qt.Vertical)
        limits.addWidget(self.bm_table)
        limits.addWidget(self.sf_table)
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.addWidget(self.structure_box)
        layout.addWidget(limits, 1)
        self.addTab(page, "Structure && limits")

        # Welding deformation ----------------------------------------------------
        self.welding_box = QGroupBox("Welding deformation of the assembled hull")
        self.welding_box.setCheckable(True)
        self.weld_joints = RecordTable("Erection joints (angular distortion)", [
            Column("x", "x [m]"),
            Column("angle_mrad", "Kink angle [mrad]",
                   tooltip="Positive = concave-up kink that lifts the hull ends; may be left "
                           "empty when the angles are identified from the survey"),
        ])
        self.weld_zones = RecordTable("Construction blocks (shrinkage curvature)", [
            Column("xa", "xa [m]"), Column("xb", "xb [m]"),
            Column("curvature_mm_per_m2", "Curvature [mm/m²]",
                   tooltip="Direct curvature; overrides the shrinkage columns"),
            Column("deck_shrinkage_mm_per_m", "Deck shrink. [mm/m]"),
            Column("bottom_shrinkage_mm_per_m", "Bottom shrink. [mm/m]"),
            Column("depth", "Depth [m]"),
        ], hint="Curvature = (deck shrinkage − bottom shrinkage) / depth, or give it directly.")
        self.weld_survey = RecordTable("Surveyed keel line (optional)", [
            Column("x", "x [m]"), Column("lift_mm", "Lift [mm]", tooltip="Upward, any datum"),
        ], hint="Laser-tracker / total-station points: import a CSV (columns x_m, lift_mm) "
                "or paste from Excel with Ctrl+V.")
        self.weld_survey.add_button("Import CSV…", self._import_survey)
        self.survey_form = FormEditor("Survey", [
            Field("interpolation", "Interpolation", "str", default="pchip",
                  tooltip="pchip | linear | spline"),
            Field("identify_joint_angles", "Identify joint angles from survey", "bool",
                  default=False),
            Field("regularization", "Regularisation", default=0.0),
        ])
        for table in (self.weld_joints, self.weld_zones, self.weld_survey):
            table.setMinimumHeight(170)
        weld_layout = QVBoxLayout(self.welding_box)
        for w in (self.weld_joints, self.weld_zones, self.weld_survey, self.survey_form):
            weld_layout.addWidget(w)
        weld_note = QLabel("The welding gaps are computed at the block positions; the spacing "
                           "design recomputes them at every candidate layout.")
        weld_note.setWordWrap(True)
        weld_note.setObjectName("hint")
        self.welding_tab = _scroll(_column_widget(self.welding_box, weld_note))
        self.addTab(self.welding_tab, "Welding deformation")

        # Settings ---------------------------------------------------------------
        methods = QGroupBox("Methods")
        self.method_rigid = QCheckBox("Rigid hull (Method A)")
        self.method_fem = QCheckBox("Elastic hull FEM (Method B)")
        row = QHBoxLayout(methods)
        row.addWidget(self.method_rigid)
        row.addWidget(self.method_fem)
        row.addStretch()
        self.settings_form = FormEditor("Settings", [
            Field("max_element_length", "Max element length [m]", default=0.5),
            Field("shear_deformation", "Shear deformation (Timoshenko)", "bool", default=True),
            Field("allowable_bending_MPa", "Allowable bending stress [MPa]"),
            Field("allowable_shear_MPa", "Allowable shear stress [MPa]"),
            Field("structure_tolerance", "Block-to-structure tolerance [m]", default=0.15),
            Field("weight_tolerance", "Weight check tolerance [-]", default=0.005),
            Field("max_iterations", "Max contact iterations", "int", default=200),
        ])
        self.addTab(_scroll(_column_widget(methods, self.settings_form)), "Settings")

        for w in (self.ship_form, self.weights, self.curve, self.sections, self.block_types,
                  self.plan_editor, self.structure_form, self.bm_table, self.sf_table,
                  self.settings_form, self.weld_joints, self.weld_zones, self.weld_survey,
                  self.survey_form):
            w.changed.connect(self.changed)
        self.curve_name.textChanged.connect(lambda _t: self.changed.emit())
        for box in (self.structure_box, self.welding_box):
            box.toggled.connect(lambda _on: self.changed.emit())
        for cb in (self.method_rigid, self.method_fem):
            cb.stateChanged.connect(lambda _s: self.changed.emit())

    # ------------------------------------------------------------------ #

    def load_sample_welding(self, x_aft: float | None = None, x_fwd: float | None = None):
        """Fill the 'Welding deformation' tab with the sample welding data.

        The sample was made for a hull from x = 0 to 120 m; for another hull the joint
        and zone positions are mapped onto its extent [x_aft, x_fwd] (six equal
        construction blocks).
        """
        joints = [dict(j) for j in SAMPLE_WELDING["joints"]]
        zones = [dict(z) for z in SAMPLE_WELDING["zones"]]
        if x_aft is not None and x_fwd is not None and x_fwd > x_aft:
            scale = (x_fwd - x_aft) / 120.0
            for j in joints:
                j["x"] = round(x_aft + j["x"] * scale, 3)
            for z in zones:
                z["xa"], z["xb"] = round(x_aft + z["xa"] * scale, 3), round(x_aft + z["xb"] * scale, 3)
        self.weld_joints.set_rows(joints)
        self.weld_zones.set_rows(zones)
        self.weld_survey.set_rows([])
        self.survey_form.set_values({})
        self.welding_box.setChecked(True)
        self.changed.emit()

    def hull_extent(self) -> tuple[float, float]:
        """x_aft and x_fwd of the Ship tab (defaults 0 and Lpp; 0–120 m if not readable)."""
        try:
            ship = self.ship_form.values()
        except ValueError:
            ship = {}
        x_aft = float(ship.get("x_aft", 0.0))
        x_fwd = float(ship.get("x_fwd", ship.get("lpp", x_aft + 120.0)))
        return (x_aft, x_fwd) if x_fwd > x_aft else (0.0, 120.0)

    def load_sample_weight_curve(self, x_aft: float | None = None,
                                 x_fwd: float | None = None) -> tuple[float, float]:
        """Replace the weights by the sample weight distribution: the lightship as a weight
        curve and three deadweight items given by their moments.

        The sample was made for a hull from x = 0 to 120 m; for another hull the positions
        are mapped onto [x_aft, x_fwd] and the intensities scaled so that every weight is
        kept. The declared weight and LCG on the Ship tab are set to the new totals, which
        are returned (t, m).
        """
        if x_aft is None or x_fwd is None:
            x_aft, x_fwd = self.hull_extent()
        scale = (x_fwd - x_aft) / 120.0

        def at(x):
            return round(x_aft + x * scale, 4)

        points = [{"x": at(x), "q_t_per_m": round(q / scale, 6)}
                  for x, q in SAMPLE_WEIGHT_CURVE["points"]]
        items = []
        for item in SAMPLE_WEIGHT_CURVE["items"]:
            row = dict(item)
            lcg = at(row["moment_tm"] / row["weight_t"])
            row["moment_tm"] = round(row["weight_t"] * lcg, 6)
            if "xa" in row:
                row["xa"], row["xb"] = at(row["xa"]), at(row["xb"])
            items.append(row)
        self.weights.set_rows(items)
        self.curve.set_rows(points)
        self.curve_name.setText(SAMPLE_WEIGHT_CURVE["name"])
        w_curve, m_curve = curve_totals(points)
        total = w_curve + sum(r["weight_t"] for r in items)
        lcg = (m_curve + sum(r["moment_tm"] for r in items)) / total
        ship = self.ship_form.values()
        ship.update(weight_t=round(total, 3), lcg=round(lcg, 3))
        self.ship_form.set_values(ship)
        self.setCurrentIndex(1)
        self.changed.emit()
        return total, lcg

    def ask_load_sample_weight_curve(self, *_):
        """Ask before replacing existing weights, load the sample and request an analysis."""
        from PyQt5.QtWidgets import QMessageBox

        try:
            has_weights = bool(self.weights.rows() or self.curve.rows())
        except ValueError:
            has_weights = True
        if has_weights and QMessageBox.question(
                self, "Sample weight distribution",
                "Replace the current weight items and weight curve by the sample weight "
                "distribution (lightship curve + 3 deadweight items given by moments)?") \
                != QMessageBox.Yes:
            return
        self.load_sample_weight_curve()
        self.sample_weight_curve_loaded.emit()

    def _import_survey(self):
        from PyQt5.QtWidgets import QFileDialog, QMessageBox

        from shipsupport.loader import read_survey_csv

        path, _ = QFileDialog.getOpenFileName(self, "Import keel survey", "",
                                              "CSV files (*.csv *.txt);;All files (*)")
        if not path:
            return
        try:
            xs, zs = read_survey_csv(path)
        except InputError as exc:
            QMessageBox.warning(self, "Survey import", str(exc))
            return
        self.weld_survey.set_rows([{"x": x, "lift_mm": z} for x, z in zip(xs, zs, strict=True)])
        self.welding_box.setChecked(True)

    def _update_type_names(self):
        try:
            names = [str(r["name"]) for r in self.block_types.rows() if r.get("name")]
        except ValueError:
            return
        self.plan_editor.set_type_names(names)

    def _update_weight_summary(self):
        """Live Σ weight and LCG of the weight items and the weight curve together."""
        try:
            total = moment = 0.0
            items = self.weights.rows()
            for r in items:
                w = r.get("weight_t", 0.0)
                total += w
                moment += r["moment_tm"] if "moment_tm" in r else w * r.get("lcg", 0.0)
            points = [p for p in self.curve.rows() if "x" in p and "q_t_per_m" in p]
            w_curve, m_curve = curve_totals(points)
            total, moment = total + w_curve, moment + m_curve
            lcg = f"{moment / total:.3f} m" if total else "—"
            curve = (f"; curve {w_curve:.1f} t, LCG {m_curve / w_curve:.3f} m "
                     f"({len(points)} points)" if w_curve else "")
            self.weight_summary.setText(f"<b>Σ weight = {total:.1f} t, LCG = {lcg}</b> "
                                        f"({len(items)} items{curve})")
        except ValueError as exc:
            self.weight_summary.setText(f"<span style='color:red'>{exc}</span>")

    def structure_spacing(self) -> float | None:
        """The floor spacing of the structure section, if one is defined (smallest distance
        between two structure positions)."""
        if not self.structure_box.isChecked():
            return None
        try:
            values = self.structure_form.values()
            if not values:
                return None
            x = structure_positions(values)
        except (ValueError, InputError):
            return None
        steps = np.diff(np.unique(x))
        return round(float(steps.min()), 9) if steps.size else None

    def load_dict(self, data: dict):
        """Fill all editors from a case dict (YAML schema)."""
        data = data or {}
        self.ship_form.set_values(data.get("ship") or {})

        weights = []
        for item in data.get("weights") or []:
            row = dict(item)
            if "x" in row:                       # point load: shown with LCG and empty extent
                row["lcg"] = row.pop("x")
            weights.append(row)
        self.weights.set_rows(weights)
        curve = data.get("weight_curve") or {}
        if isinstance(curve, list):
            curve = {"points": curve}
        self.curve_name.blockSignals(True)
        self.curve_name.setText(str(curve.get("name") or ""))
        self.curve_name.blockSignals(False)
        self.curve.set_rows(list(curve.get("points") or []))
        self.sections.set_rows(data.get("sections") or [])

        blocks = data.get("blocks") or {}
        types = []
        for name, spec in (blocks.get("types") or {}).items():
            row = dict(spec or {})
            row["name"] = name
            row["layers"] = format_layers(row.get("layers"))
            types.append(row)
        self.block_types.set_rows(types)
        self.plan_editor.set_plan(blocks.get("plan"), next(iter(blocks.get("types") or {}), ""))
        self.spacing_design = dict(data["spacing_design"]) if data.get("spacing_design") else None

        structure = data.get("structure")
        self.structure_box.setChecked(bool(structure))
        self.structure_form.set_values(structure or {})

        permissible = data.get("permissible") or {}
        self.bm_table.set_rows(permissible.get("bm") or [])
        self.sf_table.set_rows(permissible.get("sf") or [])

        settings = {**DEFAULT_SETTINGS, **(data.get("settings") or {})}
        methods = settings.get("methods") or []
        self.method_rigid.setChecked("rigid" in methods)
        self.method_fem.setChecked("fem" in methods)
        self.settings_form.set_values(settings)

        welding = data.get("welding_deformation") or {}
        self.welding_box.setChecked(bool(welding))
        self.weld_joints.set_rows(welding.get("joints") or [])
        self.weld_zones.set_rows(welding.get("zones") or [])
        survey = welding.get("measured") or {}
        if survey.get("file"):
            # survey files are read by MainWindow.open_path (inline_data_files); a raw
            # reference here cannot be shown in the table
            raise ValueError("welding_deformation.measured.file must be read with "
                             "inline_data_files() before loading the case into the editor")
        xs, lifts = list(survey.get("x") or []), list(survey.get("lift_mm") or [])
        if len(xs) != len(lifts):
            raise ValueError("welding_deformation.measured: x and lift_mm need the same length")
        self.weld_survey.set_rows([{"x": x, "lift_mm": z} for x, z in zip(xs, lifts, strict=True)])
        self.survey_form.set_values(survey)
        self.setCurrentIndex(0)

    def to_dict(self) -> dict:
        """Collect a case dict (YAML schema). Raises ValueError on unreadable cells."""
        data: dict = {"ship": self.ship_form.values()}

        weights = []
        for i, w in enumerate(self.weights.rows()):
            has_a, has_b = "xa" in w, "xb" in w
            if has_a != has_b:
                raise ValueError(f"Weight items, row {i + 1}: give both xa and xb, "
                                 f"or neither for a point load")
            if "lcg" not in w and "moment_tm" not in w:
                raise ValueError(f"Weight items, row {i + 1}: give the LCG or the moment")
            if not has_a and "lcg" in w and "moment_tm" not in w:
                w["x"] = w.pop("lcg")           # point load
            weights.append(w)
        if weights:
            data["weights"] = weights
        points = self.curve.rows()
        for i, p in enumerate(points):
            if "x" not in p or "q_t_per_m" not in p:
                raise ValueError(f"Weight curve, row {i + 1}: give x and q")
        if points:
            curve = {"points": points}
            if self.curve_name.text().strip():
                curve = {"name": self.curve_name.text().strip(), **curve}
            data["weight_curve"] = curve
        data["sections"] = self.sections.rows()

        types = {}
        for i, t in enumerate(self.block_types.rows()):
            name = t.pop("name", None)
            if not name:
                raise ValueError(f"Block types, row {i + 1}: every type needs a name")
            if name in types:
                raise ValueError(f"Block types: duplicate name '{name}'")
            if "layers" in t:
                t["layers"] = parse_layers(t["layers"], f"Block type '{name}'")
            types[name] = t
        plan = self.plan_editor.plan()
        if not plan["type"]:
            plan["type"] = next(iter(types), "")
        data["blocks"] = {"types": types, "plan": plan}

        if self.structure_box.isChecked():
            structure = self.structure_form.values()
            if structure:
                data["structure"] = structure

        permissible = {}
        if self.bm_table.rows():
            permissible["bm"] = self.bm_table.rows()
        if self.sf_table.rows():
            permissible["sf"] = self.sf_table.rows()
        if permissible:
            data["permissible"] = permissible

        settings = self.settings_form.values()
        settings["methods"] = [m for m, cb in (("rigid", self.method_rigid), ("fem", self.method_fem))
                               if cb.isChecked()]
        if not settings["methods"]:
            raise ValueError("Settings: select at least one method")
        data["settings"] = settings

        if self.welding_box.isChecked():
            welding = {}
            if self.weld_joints.rows():
                welding["joints"] = self.weld_joints.rows()
            if self.weld_zones.rows():
                welding["zones"] = self.weld_zones.rows()
            survey = self.weld_survey.rows()
            if survey:
                for i, row in enumerate(survey):
                    if "x" not in row or "lift_mm" not in row:
                        raise ValueError(f"Surveyed keel line, row {i + 1}: give x and lift")
                welding["measured"] = {"x": [r["x"] for r in survey],
                                       "lift_mm": [r["lift_mm"] for r in survey],
                                       **self.survey_form.values()}
            if welding:
                data["welding_deformation"] = welding

        if self.spacing_design is not None:
            data["spacing_design"] = dict(self.spacing_design)
        return data

    def auto_spacing_design(self) -> dict:
        """spacing_design for the Auto-calculate spacing button: the loaded settings, or
        defaults with the floor spacing as the minimum spacing; always enabled."""
        if self.spacing_design is not None:
            return {**self.spacing_design, "enabled": True}
        design = {"enabled": True}
        floor = self.structure_spacing()
        if floor is not None:
            design["min_spacing"] = floor
        return design
