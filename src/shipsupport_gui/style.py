"""Visual style of the GUI: one palette, one stylesheet and the application icon.

Colours are defined once here and used by the stylesheet, the status badge and the
icon, so that the whole program has a consistent look (light theme, blue accent).
"""

from __future__ import annotations

from PyQt5.QtCore import QPointF, QRectF, Qt
from PyQt5.QtGui import QBrush, QColor, QFont, QIcon, QPainter, QPalette, QPen, QPixmap, QPolygonF
from PyQt5.QtWidgets import QApplication

ACCENT = "#1f4e79"          # navy blue: titles, primary buttons, selected tabs
ACCENT_HOVER = "#2a6496"
ACCENT_LIGHT = "#e3ecf5"
WINDOW = "#f4f6f9"
SURFACE = "#ffffff"
BORDER = "#cfd6df"
TEXT = "#1d2733"
MUTED = "#6b7785"
RUN = "#16a34a"             # bright green: Run analysis
RUN_HOVER = "#15803d"
DESIGN = "#1e88e5"          # bright blue: Auto-calculate spacing
DESIGN_HOVER = "#1565c0"

# status badge colours: (background, text)
STATE_COLOURS = {
    "green": ("#d9f0dd", "#1e6b30"),
    "red": ("#fbdcdc", "#a11d1d"),
    "blue": ("#dde9f6", ACCENT),
    "gray": ("#e6e9ed", "#4a5560"),
    "black": ("#e6e9ed", TEXT),
}

STYLESHEET = f"""
QMainWindow, QDialog {{ background: {WINDOW}; }}
QWidget {{ color: {TEXT}; }}
QToolTip {{ background: {SURFACE}; color: {TEXT}; border: 1px solid {BORDER}; padding: 4px; }}

QMenuBar {{ background: {SURFACE}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item {{ padding: 4px 10px; background: transparent; }}
QMenuBar::item:selected {{ background: {ACCENT_LIGHT}; border-radius: 4px; }}
QMenu {{ background: {SURFACE}; border: 1px solid {BORDER}; padding: 4px; }}
QMenu::item {{ padding: 5px 24px 5px 22px; border-radius: 4px; }}
QMenu::item:selected {{ background: {ACCENT_LIGHT}; color: {TEXT}; }}
QMenu::separator {{ height: 1px; background: {BORDER}; margin: 4px 8px; }}

QToolBar#main_toolbar {{ background: {SURFACE}; border: none; border-bottom: 1px solid {BORDER};
                        padding: 4px 6px; spacing: 6px; }}
QToolBar#main_toolbar QToolButton {{ padding: 5px 10px; border-radius: 5px; }}
QToolBar#main_toolbar QToolButton:hover {{ background: {ACCENT_LIGHT}; }}
QToolBar#main_toolbar QToolButton#runTool, QToolBar#main_toolbar QToolButton#designTool {{
            color: white; font-weight: 600; padding: 6px 14px; border-radius: 6px;
            margin: 0 3px; }}
QToolBar#main_toolbar QToolButton#runTool {{ background: {RUN}; }}
QToolBar#main_toolbar QToolButton#runTool:hover {{ background: {RUN_HOVER}; }}
QToolBar#main_toolbar QToolButton#designTool {{ background: {DESIGN}; }}
QToolBar#main_toolbar QToolButton#designTool:hover {{ background: {DESIGN_HOVER}; }}
QToolBar#main_toolbar QToolButton#runTool:disabled,
QToolBar#main_toolbar QToolButton#designTool:disabled {{ background: #b8c2cc; }}

QSplitter::handle {{ background: {WINDOW}; }}
QSplitter::handle:horizontal {{ width: 6px; }}
QSplitter::handle:vertical {{ height: 6px; }}

QTabWidget::pane {{ border: 1px solid {BORDER}; border-radius: 6px; background: {SURFACE};
                   top: -1px; }}
QTabBar::tab {{ background: transparent; color: {MUTED}; padding: 7px 14px; margin-right: 2px;
               border: 1px solid transparent; border-bottom: none;
               border-top-left-radius: 6px; border-top-right-radius: 6px; }}
QTabBar::tab:hover {{ color: {TEXT}; background: {ACCENT_LIGHT}; }}
QTabBar::tab:selected {{ color: {ACCENT}; background: {SURFACE}; border-color: {BORDER};
                        border-top: 3px solid {ACCENT}; }}

QGroupBox {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 6px;
            margin-top: 14px; padding: 12px 8px 8px 8px; font-weight: 600; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; color: {ACCENT}; }}
QGroupBox::indicator {{ width: 14px; height: 14px; }}

QLineEdit {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 4px;
             padding: 4px 6px; min-height: 18px; }}
QLineEdit:focus {{ border: 1px solid {ACCENT}; }}
QSpinBox, QComboBox {{ min-height: 24px; }}
QLineEdit:disabled {{ background: {WINDOW}; color: {MUTED}; }}

QPushButton {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 5px;
              padding: 5px 12px; }}
QPushButton:hover {{ background: {ACCENT_LIGHT}; border-color: #a9bccf; }}
QPushButton:pressed {{ background: #cddcec; }}
QPushButton:disabled {{ color: {MUTED}; background: {WINDOW}; }}
QPushButton#primary {{ background: {DESIGN}; color: white; border: none; font-weight: 600;
                      border-radius: 6px; padding: 8px 18px; font-size: 10.5pt; }}
QPushButton#primary:hover {{ background: {DESIGN_HOVER}; }}
QPushButton#primary:disabled {{ background: #b8c2cc; }}

QTableView, QTableWidget {{ background: {SURFACE}; alternate-background-color: #f7f9fb;
            gridline-color: #e3e7ec; border: 1px solid {BORDER}; border-radius: 4px;
            selection-background-color: #cfe0f1; selection-color: {TEXT}; }}
QHeaderView::section {{ background: #eef2f6; color: {TEXT}; font-weight: 600; padding: 5px 6px;
                       border: none; border-right: 1px solid #dde3ea;
                       border-bottom: 1px solid {BORDER}; }}
QTableCornerButton::section {{ background: #eef2f6; border: none; }}

QPlainTextEdit, QListWidget {{ background: {SURFACE}; border: 1px solid {BORDER};
                              border-radius: 4px; }}
QScrollArea {{ border: none; background: transparent; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}

QLabel#hint {{ color: {MUTED}; }}
QLabel#sectionTitle {{ color: {ACCENT}; font-weight: 600; font-size: 10pt; }}
QLabel#emptyState {{ color: {MUTED}; font-size: 11pt; padding: 20px; }}
QLabel#verdict {{ border-radius: 9px; padding: 2px 10px; font-weight: 600; }}
QStatusBar {{ background: {SURFACE}; border-top: 1px solid {BORDER}; }}
"""


def apply_style(app: QApplication) -> None:
    """Fusion style with the program's palette, font and stylesheet."""
    app.setStyle("Fusion")
    palette = QPalette()
    for role, colour in ((QPalette.Window, WINDOW), (QPalette.Base, SURFACE),
                         (QPalette.AlternateBase, "#f7f9fb"), (QPalette.Button, SURFACE),
                         (QPalette.Text, TEXT), (QPalette.WindowText, TEXT),
                         (QPalette.ButtonText, TEXT), (QPalette.Highlight, "#cfe0f1"),
                         (QPalette.HighlightedText, TEXT), (QPalette.ToolTipBase, SURFACE),
                         (QPalette.ToolTipText, TEXT), (QPalette.PlaceholderText, "#9aa4af")):
        palette.setColor(role, QColor(colour))
    app.setPalette(palette)
    font = QFont("Segoe UI", 9)
    font.setStyleHint(QFont.SansSerif)
    app.setFont(font)
    app.setStyleSheet(STYLESHEET)
    app.setWindowIcon(app_icon())


def app_icon() -> QIcon:
    """Program icon drawn in code: a hull resting on three keel blocks."""
    icon = QIcon()
    for size in (16, 24, 32, 48, 64, 128, 256):
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.transparent)
        p = QPainter(pixmap)
        p.setRenderHint(QPainter.Antialiasing)
        s = size / 64.0
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(QColor(ACCENT)))
        p.drawRoundedRect(QRectF(0, 0, size, size), 12 * s, 12 * s)
        hull = QPolygonF([QPointF(8 * s, 16 * s), QPointF(56 * s, 16 * s), QPointF(50 * s, 36 * s),
                          QPointF(14 * s, 36 * s)])
        p.setBrush(QBrush(QColor("#e3ecf5")))
        p.drawPolygon(hull)
        p.setBrush(QBrush(QColor("#e8a33d")))
        for x in (17, 29, 41):
            p.drawRect(QRectF(x * s, 39 * s, 7 * s, 9 * s))
        p.setPen(QPen(QColor("#e3ecf5"), max(1.0, 2.5 * s)))
        p.drawLine(QPointF(10 * s, 50 * s), QPointF(54 * s, 50 * s))
        p.end()
        icon.addPixmap(pixmap)
    return icon


def _white_icon(draw) -> QIcon:
    """A white icon for the coloured toolbar buttons, drawn on a 32-unit canvas."""
    icon = QIcon()
    for size in (16, 24, 32, 48):
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.transparent)
        p = QPainter(pixmap)
        p.setRenderHint(QPainter.Antialiasing)
        p.scale(size / 32.0, size / 32.0)
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(QColor("white")))
        draw(p)
        p.end()
        icon.addPixmap(pixmap)
    return icon


def run_icon() -> QIcon:
    """Play triangle: run the analysis."""
    return _white_icon(lambda p: p.drawPolygon(QPolygonF(
        [QPointF(9, 6), QPointF(26, 16), QPointF(9, 26)])))


def design_icon() -> QIcon:
    """Three blocks with a double arrow between them: design the block spacing."""
    def draw(p):
        for x in (3, 13.5, 24):
            p.drawRoundedRect(QRectF(x, 18, 5, 9), 1, 1)
        p.setPen(QPen(QColor("white"), 2.2, Qt.SolidLine, Qt.RoundCap))
        p.drawLine(QPointF(7, 10), QPointF(25, 10))
        p.setPen(Qt.NoPen)
        p.drawPolygon(QPolygonF([QPointF(3, 10), QPointF(9, 6), QPointF(9, 14)]))
        p.drawPolygon(QPolygonF([QPointF(29, 10), QPointF(23, 6), QPointF(23, 14)]))
    return _white_icon(draw)


def badge_style(colour: str) -> str:
    """Stylesheet of the verdict badge in the status bar for a state colour name."""
    background, text = STATE_COLOURS.get(colour, STATE_COLOURS["black"])
    return f"background: {background}; color: {text};"
