"""Block plan editor: count, spacing, strength and gap of every block, and the automatic
spacing design button (part of the 'Blocks' input tab)."""

from __future__ import annotations

import numpy as np
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (QAbstractItemView, QComboBox, QGridLayout, QGroupBox, QHeaderView,
                             QLabel, QLineEdit, QPushButton, QSpinBox, QTableWidgetItem,
                             QVBoxLayout, QWidget)

from .widgets import PasteTable, format_value

COL_BLOCK, COL_X, COL_SPACING, COL_STRENGTH, COL_GAP, COL_R, COL_U = range(7)
HEADERS = ["Block", "x [m]", "Spacing [m]", "Strength [t]", "Gap [mm]", "R [t]", "R / S"]
HEADER_TIPS = ["block label, B01 = aftmost", "position from AP (from the first x and the spacings)",
               "distance to the next block (editable)",
               "force the block can carry (editable; empty = block type capacity)",
               "initial gap between block top and keel (editable)",
               "reaction of the block in the last analysis",
               "reaction / strength in the last analysis (red > 1)"]
EDITABLE = (COL_SPACING, COL_STRENGTH, COL_GAP)
READ_ONLY_BG = QColor("#f3f5f8")


def _float(text: str, what: str) -> float:
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"Block plan: {what} '{text}' is not a valid number") from None


def _compact(values: list):
    """One value if all entries are equal, else the list (as in the YAML schema)."""
    return values[0] if values and all(v == values[0] for v in values) else values


class BlockPlanEditor(QWidget):
    """Edits blocks.plan: the number of blocks, the first position, and the spacing,
    strength and initial gap of every block. The table is also the main view of the
    result: it shows the reaction and reaction / strength of every block, and after an
    automatic spacing design the designed spacing."""

    changed = pyqtSignal()
    auto_requested = pyqtSignal()          # "Auto-calculate spacing"

    def __init__(self, parent=None):
        super().__init__(parent)
        self._updating = False

        self.plan_box = QGroupBox("Block plan: number, spacing and strength of every block")
        self.count = QSpinBox()
        self.count.setRange(2, 1000)
        self.count.setValue(44)
        self.count.setToolTip("Number of blocks; rows are added (copying the last row) or "
                              "removed at the forward end")
        self.start_x = QLineEdit("6.4")
        self.start_x.setToolTip("x of the first (aftmost) block, from AP")
        self.block_type = QComboBox()
        self.block_type.setEditable(True)
        self.block_type.setToolTip("Block type of every block (types table above)")
        self.all_spacing = QLineEdit("2.4")
        self.all_strength = QLineEdit("200")
        set_spacing = QPushButton("Set all spacings")
        set_strength = QPushButton("Set all strengths")
        set_strength.setToolTip("An empty strength means the block type's capacity")
        grid = QGridLayout()
        for row, (label_a, widget_a, label_b, widget_b) in enumerate((
                ("Number of blocks", self.count, "First block x [m]", self.start_x),
                ("Block type", self.block_type, "", None),
                ("Spacing [m]", self.all_spacing, "", set_spacing),
                ("Strength [t]", self.all_strength, "", set_strength))):
            grid.addWidget(QLabel(label_a), row, 0)
            grid.addWidget(widget_a, row, 1)
            if label_b:
                grid.addWidget(QLabel(label_b), row, 2)
                grid.addWidget(widget_b, row, 3)
            elif widget_b is not None:
                grid.addWidget(widget_b, row, 2, 1, 2)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)

        self.table = PasteTable(0, len(HEADERS))
        for c, (text, tip) in enumerate(zip(HEADERS, HEADER_TIPS, strict=True)):
            item = QTableWidgetItem(text)
            item.setToolTip(tip)
            self.table.setHorizontalHeaderItem(c, item)
        self.table.setSelectionMode(QAbstractItemView.ContiguousSelection)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(24)
        header = self.table.horizontalHeader()
        header.setMinimumSectionSize(56)
        header.setSectionResizeMode(QHeaderView.Stretch)
        for c in (COL_BLOCK, COL_X):
            header.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        self.table.setMinimumHeight(300)
        self.table.setToolTip("Edit the spacing to the next block, the strength (empty = type "
                              "capacity) and the initial gap of every block; paste columns from "
                              "a spreadsheet with Ctrl+V")
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setTextFormat(Qt.RichText)

        self.auto_button = QPushButton("Auto-calculate spacing")
        self.auto_button.setObjectName("primary")
        self.auto_button.setToolTip("Find the spacing of the blocks (number and strengths kept) "
                                    "at which no reaction exceeds its block strength, and "
                                    "write it into the block plan (Ctrl+F5)")
        self.status = QLabel("Auto-calculate spacing designs the spacing of these blocks so "
                             "that every reaction fits its block strength; the result is "
                             "written into the table.")
        self.status.setObjectName("hint")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.RichText)

        plan_layout = QVBoxLayout(self.plan_box)
        plan_layout.addLayout(grid)
        plan_layout.addWidget(self.table, 1)
        plan_layout.addWidget(self.info)
        plan_layout.addWidget(self.auto_button)
        plan_layout.addWidget(self.status)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.plan_box)

        self.count.valueChanged.connect(self._resize)
        set_spacing.clicked.connect(lambda: self._set_all(COL_SPACING, self.all_spacing))
        set_strength.clicked.connect(lambda: self._set_all(COL_STRENGTH, self.all_strength))
        self.auto_button.clicked.connect(self.auto_requested)
        self.table.itemChanged.connect(self._edited)
        self.start_x.textChanged.connect(self._edited)
        self.block_type.currentTextChanged.connect(lambda _t: self.changed.emit())
        self._resize(self.count.value())

    # ------------------------------------------------------------------ table

    def _item(self, row: int, col: int, text: str = "") -> QTableWidgetItem:
        item = self.table.item(row, col)
        if item is None:
            item = QTableWidgetItem(text)
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            if col not in EDITABLE:
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setBackground(READ_ONLY_BG)
            self.table.setItem(row, col, item)
        else:
            item.setText(text)
        return item

    def _resize(self, count: int):
        rows = self.table.rowCount()
        if count == rows:
            return
        self._updating = True
        self.table.blockSignals(True)
        last_spacing = self._text(rows - 2, COL_SPACING) if rows >= 2 else self.all_spacing.text()
        last_strength = self._text(rows - 1, COL_STRENGTH) if rows >= 1 else ""
        self.table.setRowCount(count)
        for r in range(rows, count):
            for c in range(len(HEADERS)):
                self._item(r, c)
            self._item(r, COL_STRENGTH, last_strength)
        for r in range(max(rows - 1, 0), count - 1):          # every gap needs a spacing
            if not self._text(r, COL_SPACING):
                self._item(r, COL_SPACING, last_spacing)
        self.table.blockSignals(False)
        self._updating = False
        self._edited()

    def _text(self, row: int, col: int) -> str:
        item = self.table.item(row, col) if row >= 0 else None
        return item.text().strip() if item else ""

    def _set_all(self, col: int, source: QLineEdit):
        text = source.text().strip()
        if text or col != COL_STRENGTH:              # an empty strength = type capacity
            _float(text, HEADERS[col])
        self.table.blockSignals(True)
        last = self.table.rowCount() - (1 if col == COL_SPACING else 0)
        for r in range(last):
            self._item(r, col, text)
        self.table.blockSignals(False)
        self._edited()

    def _edited(self, *_):
        if self._updating:
            return
        self._refresh()
        self.changed.emit()

    def _refresh(self):
        """Recompute labels, x positions, the last row and the summary line."""
        self._updating = True
        self.table.blockSignals(True)
        n = self.table.rowCount()
        try:
            x = _float(self.start_x.text().strip() or "0", "first block x")
        except ValueError:
            x = None
        total_strength, missing, positions = 0.0, 0, []
        for r in range(n):
            self._item(r, COL_BLOCK, f"B{r + 1:02d}")
            self._item(r, COL_X, "" if x is None else f"{x:.3f}")
            positions.append(x)
            if r == n - 1:
                spacing_item = self._item(r, COL_SPACING, "")
                spacing_item.setFlags(spacing_item.flags() & ~Qt.ItemIsEditable)
                spacing_item.setBackground(READ_ONLY_BG)
            elif x is not None:
                try:
                    x = x + _float(self._text(r, COL_SPACING), "spacing")
                except ValueError:
                    x = None
            try:
                total_strength += _float(self._text(r, COL_STRENGTH), "strength")
            except ValueError:
                missing += 1
        self.table.blockSignals(False)
        self._updating = False
        if positions and positions[0] is not None and positions[-1] is not None:
            gaps = np.diff(positions)
            strength = (f"total strength {total_strength:.0f} t" if not missing else
                        f"{missing} block(s) with the type capacity as strength")
            self.info.setText(
                f"{n} blocks, block line {positions[0]:.2f} – {positions[-1]:.2f} m "
                f"({positions[-1] - positions[0]:.2f} m), spacing {gaps.min():.2f} – "
                f"{gaps.max():.2f} m, {strength}")
        else:
            self.info.setText("<span style='color:red'>Enter a number for the first block x "
                              "and every spacing.</span>")

    # ------------------------------------------------------------------ data

    def plan(self) -> dict:
        """blocks.plan as in the YAML schema (raises ValueError on unreadable cells)."""
        n = self.table.rowCount()

        def column(col, what, rows):
            return [None if not self._text(r, col) else
                    _float(self._text(r, col), f"B{r + 1:02d} {what}") for r in range(rows)]

        spacing = column(COL_SPACING, "spacing", n - 1)
        missing = [f"B{r + 1:02d}" for r, s in enumerate(spacing) if s is None]
        if missing:
            raise ValueError(f"Block plan: give the spacing to the next block for "
                             f"{', '.join(missing)}")
        plan = {"type": self.block_type.currentText().strip(),
                "start_x": _float(self.start_x.text().strip(), "first block x"), "count": n,
                "spacing": _compact(spacing)}
        strength = column(COL_STRENGTH, "strength", n)
        if any(s is not None for s in strength):
            plan["strength_t"] = _compact(strength)
        gap = column(COL_GAP, "gap", n)
        if any(g for g in gap):
            plan["gap_mm"] = _compact(gap)
        return plan

    def set_type_names(self, names: list[str]):
        """Offer the block type names; the current text is kept."""
        current = self.block_type.currentText()
        self.block_type.blockSignals(True)
        self.block_type.clear()
        self.block_type.addItems(names)
        self.block_type.setEditText(current or (names[0] if names else ""))
        self.block_type.blockSignals(False)

    def set_plan(self, plan: dict | None, default_type: str = ""):
        plan = plan or {}
        spacing, strength, gap = plan.get("spacing"), plan.get("strength_t"), plan.get("gap_mm")
        count = plan.get("count") or next(
            (len(v) + (1 if v is spacing else 0) for v in (strength, gap, spacing)
             if isinstance(v, list)), 44)
        count = int(count)
        if spacing is None:
            spacing = float(self.all_spacing.text() or 2.4)
        gaps = spacing if isinstance(spacing, list) else [spacing] * (count - 1)
        strengths = strength if isinstance(strength, list) else [strength] * count
        initial = gap if isinstance(gap, list) else [gap] * count
        self.set_blocks(float(plan.get("start_x", 6.4)), gaps, strengths,
                        str(plan.get("type") or default_type), initial)

    def set_blocks(self, start_x: float, gaps, strengths, type_name: str | None = None,
                   initial_gaps=None):
        """Fill the table: first block x, spacing to the next block, strength [t] (None =
        type capacity) and initial gap [mm] of every block."""
        initial_gaps = initial_gaps if initial_gaps is not None else [None] * len(strengths)
        self._updating = True
        self.count.blockSignals(True)
        self.count.setValue(len(strengths))
        self.count.blockSignals(False)
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.setRowCount(len(strengths))
        for r, (s, g) in enumerate(zip(strengths, initial_gaps, strict=True)):
            for c in range(len(HEADERS)):
                self._item(r, c)
            if r < len(gaps):
                self._item(r, COL_SPACING, format_value(round(float(gaps[r]), 6), "float"))
            self._item(r, COL_STRENGTH, "" if s is None else format_value(float(s), "float"))
            self._item(r, COL_GAP, "" if not g else format_value(float(g), "float"))
        self.table.blockSignals(False)
        self.start_x.blockSignals(True)
        self.start_x.setText(format_value(round(float(start_x), 6), "float"))
        self.start_x.blockSignals(False)
        if type_name is not None:
            self.block_type.blockSignals(True)
            self.block_type.setEditText(type_name)
            self.block_type.blockSignals(False)
        self._updating = False
        self._refresh()
        self.changed.emit()

    def x_positions(self) -> np.ndarray | None:
        try:
            return np.array([_float(self._text(r, COL_X), "x") for r in range(self.table.rowCount())])
        except ValueError:
            return None

    # ------------------------------------------------------------------ results

    def show_reactions(self, x, reactions_t, utilisation, target: float = 1.0) -> bool:
        """Show the reaction and R/S of every block if the result belongs to this plan."""
        own = self.x_positions()
        self.table.blockSignals(True)
        matches = own is not None and x is not None and len(own) == len(x) and \
            np.allclose(own, x, atol=1e-3)
        for r in range(self.table.rowCount()):
            if matches:
                u = utilisation[r]
                self._item(r, COL_R, f"{reactions_t[r]:.1f}")
                item = self._item(r, COL_U, "" if not np.isfinite(u) else f"{u:.3f}")
                item.setBackground(READ_ONLY_BG if not np.isfinite(u) else
                                   QColor("#f4c7c3") if u > target else
                                   QColor("#fce8b2") if u > 0.9 * target else QColor("#d9ead3"))
            else:
                self._item(r, COL_R, "")
                self._item(r, COL_U, "").setBackground(READ_ONLY_BG)
        self.table.blockSignals(False)
        return matches

    def set_status(self, html: str):
        self.status.setText(html)
