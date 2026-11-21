"""Main window: case editor on the left, results on the right."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml
from PyQt5.QtCore import QSettings, Qt, QUrl
from PyQt5.QtGui import QDesktopServices, QKeySequence
from PyQt5.QtWidgets import (QAction, QFileDialog, QLabel, QMainWindow, QMessageBox, QSplitter,
                             QStyle)

import shipsupport
from shipsupport import report
from shipsupport.analysis import AnalysisResult
from shipsupport.loader import inline_data_files
from shipsupport.units import MM, to_tonnes

from .input_panels import DEFAULT_SETTINGS, CaseEditor
from .results_panels import ResultsView
from .style import app_icon, badge_style, design_icon, run_icon
from .worker import AnalysisThread

PROJECT_ROOT = Path(shipsupport.__file__).resolve().parents[2]
EXAMPLE_CASE = PROJECT_ROOT / "examples" / "example_case.yaml"
EXAMPLES = [
    ("Example vessel (44 blocks of 200 t)", EXAMPLE_CASE),
    ("Block plan with block strengths + spacing design",
     PROJECT_ROOT / "examples" / "spacing_case.yaml"),
    ("Weight curve + items with moments + spacing design",
     PROJECT_ROOT / "examples" / "weight_curve_case.yaml"),
    ("Welded hull — welding model", PROJECT_ROOT / "examples" / "welding_case.yaml"),
    ("Welded hull — keel survey (CSV), identified joint angles",
     PROJECT_ROOT / "examples" / "welding_survey_case.yaml"),
]
THESIS_PDF = PROJECT_ROOT / "thesis" / "output" / "thesis.pdf"
FILE_FILTER = "Case files (*.yaml *.yml);;All files (*)"

NEW_CASE = {
    "ship": {"name": "New ship"},
    "sections": [],
    "blocks": {"types": {}},
    "settings": DEFAULT_SETTINGS,
}


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.path: Path | None = None
        self.modified = False
        self.result: AnalysisResult | None = None
        self.thread: AnalysisThread | None = None
        self.settings = QSettings("shipsupport", "shipsupport-gui")

        self.editor = CaseEditor()
        self.results = ResultsView()
        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.addWidget(self.editor)
        self.splitter.addWidget(self.results)
        self.splitter.setChildrenCollapsible(False)
        self.splitter.setStretchFactor(0, 9)
        self.splitter.setStretchFactor(1, 11)
        self.results.setMinimumWidth(450)
        self._sized = False
        self.setCentralWidget(self.splitter)

        self.state_label = QLabel()
        self.state_label.setObjectName("verdict")
        self.statusBar().addPermanentWidget(self.state_label)
        self.setWindowIcon(app_icon())
        self._create_actions()
        self.editor.changed.connect(self._on_input_changed)

        self.resize(1550, 950)
        self.new_case(confirm=False)

    # ------------------------------------------------------------------ UI

    def _action(self, text, slot, shortcut=None, icon=None, tip=""):
        action = QAction(text, self)
        action.triggered.connect(slot)
        if shortcut:
            action.setShortcut(QKeySequence(shortcut))
        if icon is not None:
            action.setIcon(self.style().standardIcon(icon))
        if tip:
            action.setStatusTip(tip)
        return action

    def _create_actions(self):
        # the lambda is needed: triggered(checked=False) must not arrive as new_case(confirm=False)
        self.new_action = self._action("&New", lambda: self.new_case(), QKeySequence.New,  # noqa: PLW0108
                                       QStyle.SP_FileIcon, "Start an empty case")
        self.open_action = self._action("&Open…", self.open_dialog, QKeySequence.Open,
                                        QStyle.SP_DialogOpenButton, "Open a YAML case file")
        self.save_action = self._action("&Save", self.save, QKeySequence.Save,
                                        QStyle.SP_DialogSaveButton, "Save the case")
        self.save_as_action = self._action("Save &as…", self.save_as, QKeySequence.SaveAs)
        self.export_action = self._action("&Export results…", self.export_results, "Ctrl+E",
                                          tip="Write summary, CSV tables and plots to a folder")
        self.export_action.setEnabled(False)
        exit_action = self._action("E&xit", self.close, QKeySequence.Quit)

        self.run_action = self._action("&Run analysis", self.run_analysis, "F5",
                                       QStyle.SP_MediaPlay, "Run the analysis (F5)")
        self.spacing_action = self._action(
            "Auto-calculate &spacing", self.auto_spacing, "Ctrl+F5",
            QStyle.SP_BrowserReload,
            "Design the spacing of the blocks so that every reaction fits its strength")
        theory_action = self._action("&Theory and methods (thesis PDF)", self.open_theory)
        theory_action.setEnabled(THESIS_PDF.exists())
        about_action = self._action("&About", self.about)

        file_menu = self.menuBar().addMenu("&File")
        for a in (self.new_action, self.open_action):
            file_menu.addAction(a)
        examples_menu = file_menu.addMenu("Open e&xample")
        for title, path in EXAMPLES:
            action = examples_menu.addAction(title)
            action.triggered.connect(lambda _checked=False, p=path: self.open_example(p))
            action.setEnabled(path.exists())
        self.sample_welding_action = self._action(
            "Add sample &welding data", self.add_sample_welding,
            tip="Fill the 'Welding deformation' tab with sample welding data and run")
        self.sample_curve_action = self._action(
            "Add sample weight &distribution curve", self.editor.ask_load_sample_weight_curve,
            tip="Replace the weights by a sample lightship curve and deadweight items, and run")
        self.editor.sample_weight_curve_loaded.connect(self._run_sample_weight_curve)
        for a in (self.sample_curve_action, self.sample_welding_action, None, self.save_action,
                  self.save_as_action, None,
                  self.export_action, None, exit_action):
            file_menu.addSeparator() if a is None else file_menu.addAction(a)
        self.results.sample_welding_requested.connect(self.add_sample_welding)
        self.editor.plan_editor.auto_requested.connect(self.auto_spacing)
        self.results.spacing_requested.connect(self.auto_spacing)
        analysis_menu = self.menuBar().addMenu("&Analysis")
        analysis_menu.addAction(self.run_action)
        analysis_menu.addAction(self.spacing_action)
        help_menu = self.menuBar().addMenu("&Help")
        help_menu.addAction(theory_action)
        help_menu.addAction(about_action)

        toolbar = self.addToolBar("Main")
        toolbar.setObjectName("main_toolbar")
        toolbar.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        toolbar.setMovable(False)
        for a in (self.new_action, self.open_action, self.save_action, None, self.run_action,
                  self.spacing_action):
            toolbar.addSeparator() if a is None else toolbar.addAction(a)
        self.run_action.setIcon(run_icon())
        self.spacing_action.setIcon(design_icon())
        toolbar.widgetForAction(self.run_action).setObjectName("runTool")       # main commands:
        toolbar.widgetForAction(self.spacing_action).setObjectName("designTool")  # bright colours

    def _update_title(self):
        name = self.path.name if self.path else "untitled"
        self.setWindowTitle(f"{name}{' *' if self.modified else ''} — Ship Support Force")

    def _set_state(self, text: str, colour: str = "black"):
        self.state_label.setText(text)
        self.state_label.setStyleSheet(badge_style(colour))
        self.state_label.setVisible(bool(text))

    def _on_input_changed(self):
        if not self.modified:
            self.modified = True
            self._update_title()
        if self.result is not None:
            self._set_state("Results out of date — press F5", "darkorange")

    # ------------------------------------------------------------------ files

    def _confirm_discard(self) -> bool:
        if not self.modified:
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes", "The case has unsaved changes. Save them first?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if answer == QMessageBox.Save:
            return self.save()
        return answer == QMessageBox.Discard

    def _load(self, data: dict, path: Path | None):
        self.editor.load_dict(data)
        self.path = path
        self.modified = False
        self.result = None
        self.export_action.setEnabled(False)
        self.results.clear()
        self._set_state("Not analysed")
        self._update_title()

    def new_case(self, confirm: bool = True):
        if confirm and not self._confirm_discard():
            return
        self._load(NEW_CASE, None)

    def open_path(self, path: Path) -> bool:
        try:
            with Path(path).open(encoding="utf-8") as f:
                data = yaml.safe_load(f)
            if not isinstance(data, dict):
                raise ValueError("the top level of the file must be a mapping")
            data = inline_data_files(data, Path(path).parent)   # survey CSV → table rows
            self._load(data, Path(path))
        except (OSError, ValueError, TypeError, AttributeError, yaml.YAMLError) as exc:
            # TypeError/AttributeError: sections of the file with an unexpected structure
            QMessageBox.critical(self, "Open failed", f"Could not read {path}:\n\n{exc}")
            return False
        self.settings.setValue("last_dir", str(Path(path).parent))
        self.statusBar().showMessage(f"Opened {path}", 5000)
        return True

    def open_dialog(self):
        if not self._confirm_discard():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open case", self._last_dir(), FILE_FILTER)
        if path:
            self.open_path(Path(path))

    def open_example(self, path: Path = EXAMPLE_CASE):
        """Open a shipped example as an untitled copy (Save asks for a new file name)."""
        if self._confirm_discard():
            if self.open_path(path):
                self.path = None   # don't overwrite the shipped example by accident
                self.modified = True
                self._update_title()

    def _run_sample_weight_curve(self):
        self.statusBar().showMessage("Sample weight distribution curve added (placeholder "
                                     "values) — running the analysis…", 5000)
        if self.thread is None or not self.thread.isRunning():
            self.run_analysis()

    def add_sample_welding(self, *_):
        """Add the sample welding deformation to the current case, then run."""
        if self.thread is not None and self.thread.isRunning():
            self.statusBar().showMessage("Wait for the running analysis to finish first.", 5000)
            return
        try:
            ship = self.editor.ship_form.values()
        except ValueError as exc:
            QMessageBox.warning(self, "Input error", str(exc))
            return
        lpp = ship.get("lpp")
        x_aft = ship.get("x_aft", 0.0)
        x_fwd = ship.get("x_fwd", lpp)
        self.editor.load_sample_welding(x_aft, x_fwd)
        self.editor.setCurrentIndex(self.editor.indexOf(self.editor.welding_tab))
        self.statusBar().showMessage("Sample welding data added (placeholder values) — running "
                                     "the analysis…", 5000)
        self.run_analysis()

    def _last_dir(self) -> str:
        return str(self.settings.value("last_dir", str(Path.home())))

    def save(self) -> bool:
        return self.save_as() if self.path is None else self._write(self.path)

    def save_as(self) -> bool:
        start = str(self.path) if self.path else str(Path(self._last_dir()) / "case.yaml")
        path, _ = QFileDialog.getSaveFileName(self, "Save case", start, FILE_FILTER)
        if not path:
            return False
        return self._write(Path(path))

    def _write(self, path: Path) -> bool:
        try:
            data = self.editor.to_dict()
        except ValueError as exc:
            QMessageBox.warning(self, "Input error", str(exc))
            return False
        try:
            with path.open("w", encoding="utf-8") as f:
                f.write("# Ship support force case — saved by shipsupport GUI\n")
                yaml.safe_dump(data, f, sort_keys=False, allow_unicode=True, width=100)
        except OSError as exc:
            QMessageBox.critical(self, "Save failed", f"Could not write {path}:\n\n{exc}")
            return False
        self.path = path
        self.modified = False
        self.settings.setValue("last_dir", str(path.parent))
        self._update_title()
        self.statusBar().showMessage(f"Saved {path}", 5000)
        return True

    def export_results(self):
        if self.result is None:
            return
        folder = QFileDialog.getExistingDirectory(self, "Export results to folder", self._last_dir())
        if not folder:
            return
        out = Path(folder)
        self.settings.setValue("last_dir", str(out))
        try:
            (out / "summary.txt").write_text(report.format_summary(self.result), encoding="utf-8")
            written = [out / "summary.txt", *report.write_csv(self.result, out),
                       report.plot_results(self.result, out / "results.png"),
                       report.plot_block_forces(self.result, out / "block_forces.png")]
            for extra in (report.plot_welding_data(self.result, out / "welding_data.png"),
                          report.plot_spacing(self.result, out / "spacing_design.png")):
                if extra:
                    written.append(extra)
        except (OSError, ValueError, RuntimeError) as exc:   # file system or plotting errors
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Export complete",
                                "Written:\n" + "\n".join(p.name for p in written) + f"\n\nto {out}")

    # ------------------------------------------------------------------ analysis

    def run_analysis(self, *_, spacing_design: bool = False):
        """Run the analysis (F5) in a background thread. The spacing design runs only with
        spacing_design=True (Auto-calculate spacing); otherwise a spacing_design section of
        the loaded case is kept for saving but not run."""
        if self.thread is not None and self.thread.isRunning():
            self.statusBar().showMessage("An analysis is already running; press F5 again when "
                                         "it has finished.", 5000)
            return
        try:
            data = self.editor.to_dict()
        except ValueError as exc:
            QMessageBox.warning(self, "Input error", str(exc))
            return
        data.pop("spacing_design", None)
        if spacing_design:
            data["spacing_design"] = self.editor.auto_spacing_design()
        for action in (self.run_action, self.spacing_action):
            action.setEnabled(False)
        self.editor.plan_editor.auto_button.setEnabled(False)
        self._set_state("Designing the block spacing…" if spacing_design else "Running…", "blue")
        if self.thread is not None:                 # release the previous, finished thread
            self.thread.deleteLater()
        self.thread = AnalysisThread(data, self)
        self.thread.succeeded.connect(self.show_result)
        self.thread.failed.connect(self.show_failure)
        self.thread.finished.connect(self._analysis_finished)
        self.thread.start()

    def _analysis_finished(self):
        for action in (self.run_action, self.spacing_action):
            action.setEnabled(True)
        self.editor.plan_editor.auto_button.setEnabled(True)

    def auto_spacing(self, *_):
        """Run the analysis with the spacing design and write the designed spacing into the
        block plan (Blocks tab)."""
        self.editor.setCurrentIndex(self.editor.indexOf(self.editor.blocks_tab))
        self.editor.plan_editor.set_status("<span style='color:blue'>Designing the spacing… "
                                           "(every candidate is analysed with the full support "
                                           "model)</span>")
        self.run_analysis(spacing_design=True)

    def _apply_spacing(self, result: AnalysisResult):
        """Write the designed positions into the block plan and report the outcome."""
        sp = result.spacing
        plan = self.editor.plan_editor
        blocks = sp.final_blocks
        plan.set_blocks(blocks[0].x, sp.spacing,
                        [None if b.strength is None else round(to_tonnes(b.strength), 6)
                         for b in blocks],
                        blocks[0].type.name, [round(b.gap / MM, 6) for b in blocks])
        u0, u1 = np.max(sp.initial_utilisation), np.max(sp.final_utilisation)
        colour = "green" if sp.fits else "red"
        plan.set_status(
            f"<span style='color:{colour}'><b>{'Fits' if sp.fits else 'Does not fit'}</b></span>"
            f" — max reaction / strength {u0:.3f} → <b>{u1:.3f}</b> (target {sp.target:g}); "
            f"spacing {sp.spacing.min():.2f} – {sp.spacing.max():.2f} m; {sp.evaluations} "
            f"analyses. The designed spacing is now in the block plan."
            + ("" if sp.fits else f"<br><span style='color:red'>{sp.message}</span>"))

    def show_result(self, result: AnalysisResult):
        self.result = result
        self.results.set_result(result)
        self.export_action.setEnabled(True)
        if result.spacing is not None:
            self._apply_spacing(result)
            sp = result.spacing
            self.editor.plan_editor.show_reactions(
                [b.x for b in sp.final_blocks], to_tonnes(sp.final_evaluation.reactions),
                sp.final_utilisation, sp.target)
        else:
            design_ = result.supports[result.design_method]
            table = result.strength[result.design_method].blocks
            self.editor.plan_editor.show_reactions(design_.block_x, to_tonnes(design_.reactions),
                                                   table.capacity_util)

        design = result.supports[result.design_method]
        i = int(np.argmax(design.reactions))
        verdicts = [s.passed for s in result.strength.values()]
        if any(v is False for v in verdicts):
            text, colour = "CHECKS FAILED", "red"
        elif any(v for v in verdicts):
            text, colour = "All checks passed", "green"
        else:
            text, colour = "No limits defined", "gray"
        extra = ""
        if result.spacing is not None:
            sp = result.spacing
            extra = (f"; designed spacing: max reaction/strength "
                     f"{np.max(sp.final_utilisation):.2f} — "
                     f"{'fits' if sp.fits else 'does not fit'}")
        self._set_state(f"{text} — max R {to_tonnes(design.reactions[i]):.1f} t at "
                        f"x = {design.block_x[i]:.2f} m ({result.design_method}){extra}; "
                        f"{len(result.warnings)} warning(s)", colour)

    def show_failure(self, message: str):
        self.result = None
        self.export_action.setEnabled(False)
        self.results.show_error(message)
        self._set_state("Analysis failed", "red")
        QMessageBox.critical(self, "Analysis failed", message)

    # ------------------------------------------------------------------ misc

    def open_theory(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(THESIS_PDF)))

    def about(self):
        QMessageBox.about(
            self, "About",
            f"<b>Ship Support Force</b> {shipsupport.__version__}<br><br>"
            "Dry-docking block reactions (rigid hull and elastic hull FEM, also on a welded "
            "hull), automatic block spacing design from the block strengths, and hull girder "
            "strength checks.<br><br>"
            "Results depend on the input data; verify against the docking plan and the "
            "governing standard before use.")

    def showEvent(self, event):
        super().showEvent(event)
        if not self._sized:   # split 45 / 55 once the real window width is known
            self._sized = True
            width = self.splitter.width()
            self.splitter.setSizes([int(width * 0.42), width - int(width * 0.42)])

    def closeEvent(self, event):
        if not self._confirm_discard():
            event.ignore()
            return
        if self.thread is not None and self.thread.isRunning():
            # the analysis cannot be interrupted; destroying a running QThread would crash
            self.statusBar().showMessage("Waiting for the running analysis to finish…")
            self.thread.wait()
        event.accept()
