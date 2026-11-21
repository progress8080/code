"""Run the analysis off the GUI thread."""

from __future__ import annotations

import traceback

from PyQt5.QtCore import QThread, pyqtSignal

from shipsupport.analysis import AnalysisResult, run_case
from shipsupport.loader import case_from_dict
from shipsupport.model import InputError
from shipsupport.solvers import ContactError


def run_analysis(data: dict) -> AnalysisResult:
    """Same pipeline as the CLI: YAML-schema dict -> Case -> run_case."""
    return run_case(case_from_dict(data))


class AnalysisThread(QThread):
    succeeded = pyqtSignal(object)   # AnalysisResult
    failed = pyqtSignal(str)

    def __init__(self, data: dict, parent=None):
        super().__init__(parent)
        self.data = data

    def run(self):
        try:
            self.succeeded.emit(run_analysis(self.data))
        except (InputError, ContactError) as exc:
            self.failed.emit(str(exc))
        except Exception:  # unexpected: show the traceback rather than crash the GUI
            self.failed.emit("Unexpected error:\n\n" + traceback.format_exc())
