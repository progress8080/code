"""Dry-docking block reactions, block spacing design and hull girder strength analysis."""

from .analysis import AnalysisResult, run_case
from .loader import case_from_dict, load_case
from .model import (BendingEnvelope, Block, BlockType, Case, InputError, Layer, Section, Settings,
                    ShearEnvelope, Ship, SpacingDesign, WeightCurve, WeightItem)

__version__ = "0.1.0"

__all__ = [
    "AnalysisResult", "run_case", "load_case", "case_from_dict",
    "BendingEnvelope", "Block", "BlockType", "Case", "InputError", "Layer", "Section",
    "Settings", "ShearEnvelope", "Ship", "SpacingDesign", "WeightCurve", "WeightItem",
]
