"""Result views: summary text, block and check tables, embedded matplotlib plots."""

from __future__ import annotations

import os

os.environ.setdefault("QT_API", "pyqt5")

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QBrush, QColor
from PyQt5.QtWidgets import (QAbstractItemView, QHBoxLayout, QHeaderView, QLabel, QListWidget,
                             QListWidgetItem, QPlainTextEdit, QPushButton, QScrollArea,
                             QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from shipsupport.analysis import AnalysisResult
from shipsupport.report import (draw_block_forces, draw_results, draw_spacing, draw_welding_data,
                                format_summary)
from shipsupport.units import KN_PER_MM, MM, MPA, TONNE_FORCE, to_tonnes

from .widgets import monospace_font

FAIL = QColor(255, 205, 205)
WARN = QColor(255, 235, 190)
PASS = QColor(210, 245, 210)
MUTED = QColor(130, 130, 130)

_DISPLAY = {"N": (TONNE_FORCE, "t"), "N·m": (TONNE_FORCE, "t·m"), "Pa": (MPA, "MPa")}


def _cell(value, fmt="{:.2f}", background: QColor | None = None, foreground: QColor | None = None):
    if isinstance(value, str):
        item = QTableWidgetItem(value)
    elif value is None or (isinstance(value, float) and np.isnan(value)):
        item = QTableWidgetItem("—")
        item.setTextAlignment(Qt.AlignCenter)
    else:
        item = QTableWidgetItem(fmt.format(value))
        item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    if background is not None:
        item.setBackground(QBrush(background))
    if foreground is not None:
        item.setForeground(QBrush(foreground))
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item


def _util_colour(util) -> QColor | None:
    if util is None or np.isnan(util):
        return None
    return FAIL if util > 1.0 else WARN if util > 0.9 else None


def _read_only_table(headers: list[str]) -> QTableWidget:
    table = QTableWidget(0, len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setEditTriggers(QAbstractItemView.NoEditTriggers)
    table.setSelectionMode(QAbstractItemView.ContiguousSelection)
    table.verticalHeader().setVisible(False)
    table.verticalHeader().setDefaultSectionSize(24)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
    table.horizontalHeader().setStretchLastSection(True)
    return table


class PlotView(QWidget):
    """A matplotlib figure with its navigation toolbar.

    The canvas never shrinks below a readable size (MIN_WIDTH pixels wide, height from the
    figure's aspect ratio); in a narrow window the chart scrolls instead of being squeezed.
    """

    MIN_WIDTH = 760

    def __init__(self, figsize=(10, 11), parent=None, action_text: str | None = None):
        super().__init__(parent)
        self.figure = Figure(figsize=figsize, facecolor="white")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(self.MIN_WIDTH, int(self.MIN_WIDTH * figsize[1] / figsize[0]))
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.canvas)
        self.scroll = scroll
        self.message = QLabel()
        self.message.setObjectName("emptyState")
        self.message.setAlignment(Qt.AlignCenter)
        self.message.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(scroll, 1)
        layout.addWidget(self.message)
        self.action = None
        if action_text:                       # optional button shown with the empty-state message
            self.action = QPushButton(action_text)
            self.action.setObjectName("primary")
            row = QHBoxLayout()
            row.addStretch()
            row.addWidget(self.action)
            row.addStretch()
            layout.addLayout(row)
            self.action.setVisible(False)

    def show_message(self, text: str, with_action: bool = False):
        self.figure.clear()
        self.canvas.draw_idle()
        self.scroll.setVisible(not text)
        self.message.setText(text)
        self.message.setVisible(bool(text))
        if self.action is not None:
            self.action.setVisible(with_action)

    def show_plot(self):
        self.scroll.setVisible(True)
        self.message.setVisible(False)
        if self.action is not None:
            self.action.setVisible(False)
        self.canvas.draw_idle()


class ResultsView(QTabWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.summary.setFont(monospace_font())
        self.blocks = _read_only_table([])
        self.checks = _read_only_table(["Method", "Check", "Demand", "Limit", "Unit",
                                        "Utilisation", "x [m]", "Status"])
        self.plots = PlotView()
        self.force_plot = PlotView(figsize=(18, 12.5))
        self.weld_data_plot = PlotView(figsize=(13, 11.5),
                                       action_text="Add sample welding data and run")
        self.spacing_plot = PlotView(figsize=(13, 14), action_text="Auto-calculate spacing")
        # emitted when the user presses the empty-state buttons (handled by the main window)
        self.sample_welding_requested = self.weld_data_plot.action.clicked
        self.spacing_requested = self.spacing_plot.action.clicked
        self.warnings = QListWidget()

        self.addTab(self.summary, "Summary")
        self.addTab(self.force_plot, "Block forces")
        self.addTab(self.plots, "Plots")
        self.addTab(self.blocks, "Block reactions")
        self.addTab(self.checks, "Checks")
        self.addTab(self.weld_data_plot, "Welding data")
        self.addTab(self.spacing_plot, "Spacing design")
        self.addTab(self.warnings, "Warnings")
        self.clear()

    def clear(self, message: str = "Run the analysis (F5) to see results."):
        self.summary.setPlainText(message)
        self.blocks.setRowCount(0)
        self.checks.setRowCount(0)
        self.warnings.clear()
        self.setTabText(self.indexOf(self.warnings), "Warnings")
        self.plots.show_message(message)
        self.force_plot.show_message(message)
        self.weld_data_plot.show_message(message)
        self.spacing_plot.show_message(message)

    def show_error(self, text: str):
        self.clear()
        self.summary.setPlainText(f"ANALYSIS FAILED\n\n{text}")
        self.setCurrentIndex(0)

    def set_result(self, result: AnalysisResult):
        self.summary.setPlainText(format_summary(result))
        self._fill_blocks(result)
        self._fill_checks(result)

        draw_results(self.plots.figure, result)
        self.plots.show_plot()
        draw_block_forces(self.force_plot.figure, result)
        self.force_plot.show_plot()
        if draw_welding_data(self.weld_data_plot.figure, result):
            self.weld_data_plot.show_plot()
        else:
            self.weld_data_plot.show_message(
                "This case has no welding deformation data. Enter it on the input tab "
                "'Welding deformation' (joint kinks, shrinkage zones or a keel survey), or add "
                "the sample data of a hull built from six welded construction blocks:",
                with_action=True)
        if draw_spacing(self.spacing_plot.figure, result):
            self.spacing_plot.show_plot()
        else:
            self.spacing_plot.show_message(
                "No spacing design in this analysis. Define the number, spacing and strength "
                "of every block in the block plan (input tab 'Blocks') and press "
                "Auto-calculate spacing to find the spacing at which every reaction fits the "
                "strength of its block:", with_action=True)

        self.warnings.clear()
        for w in result.warnings:
            QListWidgetItem(w, self.warnings).setBackground(QBrush(WARN))
        if not result.warnings:
            QListWidgetItem("No warnings.", self.warnings)
        self.setTabText(self.indexOf(self.warnings),
                        f"Warnings ({len(result.warnings)})" if result.warnings else "Warnings")

    def _fill_blocks(self, result: AnalysisResult):
        design = result.design_method
        table = result.strength[design].blocks
        springs = result.springs
        has_weld = result.weld_gap is not None
        headers = (["Block", "x [m]", "Gap [mm]"]
                   + (["Weld gap [mm]"] if has_weld else [])
                   + ["k [kN/mm]"] + [f"R {m} [t]" for m in result.block_supports]
                   + ["p [MPa]", "p util", "Capacity util", "δ [mm]", "State", "Off structure [m]"])
        self.blocks.clear()
        self.blocks.setColumnCount(len(headers))
        self.blocks.setHorizontalHeaderLabels(headers)
        self.blocks.setRowCount(table.x.size)
        delta = result.supports[design].delta
        for i in range(table.x.size):
            state = table.state[i]
            muted = MUTED if state == "open" else None
            cells = [_cell(table.label[i]), _cell(table.x[i], "{:.3f}"),
                     _cell(springs.gap[i] / MM, "{:.2f}")]
            if has_weld:
                cells.append(_cell(result.weld_gap[i] / MM, "{:.2f}"))
            cells.append(_cell(springs.k[i] / KN_PER_MM, "{:.1f}"))
            cells += [_cell(to_tonnes(r.reactions[i]), foreground=muted)
                      for r in result.block_supports.values()]
            cells += [
                _cell(table.pressure[i] / MPA, "{:.3f}"),
                _cell(table.pressure_util[i], "{:.3f}", _util_colour(table.pressure_util[i])),
                _cell(table.capacity_util[i], "{:.3f}", _util_colour(table.capacity_util[i])),
                _cell(delta[i] / MM, "{:.3f}"),
                _cell(state, background=FAIL if state == "crushed" else None, foreground=muted),
                _cell(table.structure_offset[i] if table.off_structure[i] else "",
                      "{:.3f}", WARN if table.off_structure[i] else None),
            ]
            for c, item in enumerate(cells):
                self.blocks.setItem(i, c, item)

    def _fill_checks(self, result: AnalysisResult):
        rows = [(m, c) for m, s in result.strength.items() for c in s.checks]
        self.checks.setRowCount(len(rows))
        for i, (method, c) in enumerate(rows):
            scale, unit = _DISPLAY[c.unit]
            status = {True: "OK", False: "FAIL", None: "no limit"}[c.passed]
            colour = {True: PASS, False: FAIL, None: None}[c.passed]
            cells = [_cell(method), _cell(c.name), _cell(c.demand / scale),
                     _cell(None if c.capacity is None else c.capacity / scale), _cell(unit),
                     _cell(c.utilisation, "{:.3f}", _util_colour(c.utilisation)),
                     _cell(c.x, "{:.3f}"), _cell(status, background=colour)]
            for col, item in enumerate(cells):
                self.checks.setItem(i, col, item)
