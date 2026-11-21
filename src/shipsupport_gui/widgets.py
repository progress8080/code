"""Reusable editor widgets: a record table with Excel copy/paste and a simple form."""

from __future__ import annotations

from dataclasses import dataclass

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QKeySequence
from PyQt5.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QFormLayout, QHBoxLayout,
                             QHeaderView, QLabel, QLineEdit, QPushButton, QTableWidget,
                             QTableWidgetItem, QVBoxLayout, QWidget)

# Value kinds: "float", "int", "str", "floats" (comma/space separated list)


def format_value(value, kind: str) -> str:
    if value is None:
        return ""
    if kind == "floats":
        values = value if isinstance(value, (list, tuple)) else [value]
        return ", ".join(f"{float(v):.12g}" for v in values)
    if kind == "float":
        return f"{float(value):.12g}"
    if kind == "int":
        return str(int(value))
    return str(value)


def parse_value(text: str, kind: str, ctx: str):
    try:
        if kind == "float":
            return float(text)
        if kind == "int":
            value = float(text)
            if not value.is_integer():          # e.g. a block count of 2.5 is an input error
                raise ValueError(f"{ctx}: '{text}' must be a whole number")
            return int(value)
        if kind == "floats":
            parts = text.replace(";", ",").replace(",", " ").split()
            return [float(p) for p in parts]
    except ValueError as exc:
        if "whole number" in str(exc):
            raise
        raise ValueError(f"{ctx}: '{text}' is not a valid number") from None
    return text


@dataclass(frozen=True)
class Column:
    key: str
    header: str
    kind: str = "float"
    tooltip: str = ""
    width: int = 0


class PasteTable(QTableWidget):
    """QTableWidget with Ctrl+C / Ctrl+V of tab-separated blocks (Excel compatible) and Delete."""

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.Copy):
            self.copy_selection()
        elif event.matches(QKeySequence.Paste):
            self.paste()
        elif event.key() in (Qt.Key_Delete, Qt.Key_Backspace) and self.state() != self.EditingState:
            items = self.selectedItems()
            self.blockSignals(True)                    # one change notification, not one per cell
            for item in items:
                item.setText("")
            self.blockSignals(False)
            if items:
                self.itemChanged.emit(items[0])
        else:
            super().keyPressEvent(event)

    def copy_selection(self):
        ranges = self.selectedRanges()
        if not ranges:
            return
        r = ranges[0]
        lines = []
        for row in range(r.topRow(), r.bottomRow() + 1):
            cells = []
            for col in range(r.leftColumn(), r.rightColumn() + 1):
                item = self.item(row, col)
                cells.append(item.text() if item else "")
            lines.append("\t".join(cells))
        QApplication.clipboard().setText("\n".join(lines))

    def paste(self):
        text = QApplication.clipboard().text()
        if not text:
            return
        rows = [line.split("\t") for line in text.rstrip("\r\n").splitlines()]
        start_row = max(self.currentRow(), 0)
        start_col = max(self.currentColumn(), 0)
        self.blockSignals(True)
        for i, cells in enumerate(rows):
            row = start_row + i
            if row >= self.rowCount():
                self.insertRow(row)
            for j, cell in enumerate(cells):
                col = start_col + j
                if col >= self.columnCount():
                    continue
                if self.item(row, col) is None:         # rows inserted by Qt have no items yet
                    self.setItem(row, col, QTableWidgetItem(""))
                self.item(row, col).setText(cell.strip())
        self.blockSignals(False)
        self.itemChanged.emit(self.item(start_row, start_col))


class RecordTable(QWidget):
    """Editable list of records (dicts). Empty cells are omitted from the records."""

    changed = pyqtSignal()

    def __init__(self, title: str, columns: list[Column], hint: str = "", parent=None):
        super().__init__(parent)
        self.title = title
        self.columns = columns

        header = QHBoxLayout()
        label = QLabel(title)
        label.setObjectName("sectionTitle")
        header.addWidget(label)
        header.addStretch()
        self._header = header
        for text, slot in (("Add", self.add_row), ("Duplicate", self.duplicate_row),
                           ("Remove", self.remove_rows)):
            self.add_button(text, slot)

        self.table = PasteTable(0, len(columns))
        self.table.setSelectionMode(QAbstractItemView.ContiguousSelection)
        self.table.verticalHeader().setDefaultSectionSize(24)
        self.table.setAlternatingRowColors(True)
        for c, col in enumerate(columns):
            item = QTableWidgetItem(col.header)
            if col.tooltip:
                item.setToolTip(col.tooltip)
            self.table.setHorizontalHeaderItem(c, item)
            if col.width:
                self.table.setColumnWidth(c, col.width)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(lambda _item: self.changed.emit())

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(header)
        if hint:
            hint_label = QLabel(hint)
            hint_label.setWordWrap(True)
            hint_label.setObjectName("hint")
            layout.addWidget(hint_label)
        layout.addWidget(self.table)

    def add_button(self, text: str, slot) -> QPushButton:
        button = QPushButton(text)
        button.clicked.connect(slot)
        self._header.addWidget(button)
        return button

    def _insert(self, row: int, record: dict | None = None):
        self.table.insertRow(row)
        for c, col in enumerate(self.columns):
            text = format_value(record.get(col.key), col.kind) if record else ""
            item = QTableWidgetItem(text)
            if col.kind in ("float", "int"):
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(row, c, item)

    def set_rows(self, records: list[dict]):
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        for record in records:
            self._insert(self.table.rowCount(), record)
        self.table.blockSignals(False)
        self.changed.emit()

    def rows(self) -> list[dict]:
        out = []
        for r in range(self.table.rowCount()):
            record = {}
            for c, col in enumerate(self.columns):
                item = self.table.item(r, c)
                text = item.text().strip() if item else ""
                if text:
                    record[col.key] = parse_value(text, col.kind,
                                                  f"{self.title}, row {r + 1}, column '{col.header}'")
            if record:
                out.append(record)
        return out

    def add_row(self):
        self.table.blockSignals(True)
        row = self.table.currentRow() + 1 if self.table.currentRow() >= 0 else self.table.rowCount()
        self._insert(row)
        self.table.blockSignals(False)
        self.table.setCurrentCell(row, 0)
        self.changed.emit()

    def duplicate_row(self):
        row = self.table.currentRow()
        if row < 0:
            return
        texts = [self.table.item(row, c).text() if self.table.item(row, c) else ""
                 for c in range(len(self.columns))]
        self.table.blockSignals(True)
        self._insert(row + 1)
        for c, text in enumerate(texts):
            self.table.item(row + 1, c).setText(text)
        self.table.blockSignals(False)
        self.changed.emit()

    def remove_rows(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        if not rows and self.table.currentRow() >= 0:
            rows = [self.table.currentRow()]
        for row in rows:
            self.table.removeRow(row)
        self.changed.emit()


@dataclass(frozen=True)
class Field:
    key: str
    label: str
    kind: str = "float"       # float, int, str, floats, bool
    default: object = None    # shown as placeholder (text allowed, e.g. "Lpp"); bool default state
    tooltip: str = ""


class FormEditor(QWidget):
    """Label/value form. Empty fields are omitted from values() so loader defaults apply."""

    changed = pyqtSignal()

    def __init__(self, title: str, fields: list[Field], parent=None):
        super().__init__(parent)
        self.title = title
        self.fields = fields
        self.editors: dict[str, QWidget] = {}
        self.setMaximumWidth(680)              # forms stay readable in a wide window
        layout = QFormLayout(self)
        layout.setVerticalSpacing(7)
        layout.setHorizontalSpacing(12)
        layout.setLabelAlignment(Qt.AlignRight)
        for f in fields:
            if f.kind == "bool":
                editor = QCheckBox()
                editor.setChecked(bool(f.default))
                editor.stateChanged.connect(lambda _s: self.changed.emit())
            else:
                editor = QLineEdit()
                if isinstance(f.default, str) and f.kind != "str":
                    editor.setPlaceholderText(f"default {f.default}")   # e.g. "default Lpp"
                elif f.default is not None:
                    editor.setPlaceholderText(f"default {format_value(f.default, f.kind)}")
                if f.kind in ("float", "int"):
                    editor.setAlignment(Qt.AlignRight)
                editor.textChanged.connect(lambda _t: self.changed.emit())
            if f.tooltip:
                editor.setToolTip(f.tooltip)
            self.editors[f.key] = editor
            layout.addRow(f.label, editor)

    def values(self) -> dict:
        out = {}
        for f in self.fields:
            editor = self.editors[f.key]
            if f.kind == "bool":
                out[f.key] = editor.isChecked()
                continue
            text = editor.text().strip()
            if text:
                out[f.key] = parse_value(text, f.kind, f"{self.title}: {f.label}")
        return out

    def set_values(self, values: dict):
        for f in self.fields:
            editor = self.editors[f.key]
            editor.blockSignals(True)
            if f.kind == "bool":
                editor.setChecked(bool(values.get(f.key, f.default)))
            else:
                editor.setText(format_value(values.get(f.key), f.kind))
            editor.blockSignals(False)
        self.changed.emit()


def monospace_font(size: int = 9) -> QFont:
    font = QFont("Consolas")
    font.setStyleHint(QFont.Monospace)
    font.setPointSize(size)
    return font
