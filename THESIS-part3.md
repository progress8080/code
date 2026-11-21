> **Part 3 of 3** of the thesis *Implementation of a Simulation System for Predicting Keel-Block Reaction Forces and Optimizing Support Spacing During Ship Construction* (Appendix B (source code)). [Part 1](THESIS-part1.md) – the front matter and Chapters 1–3 · [Part 2](THESIS-part2.md) – Chapters 4–7, the references and Appendix A.

# Appendix B – Source Code

This appendix lists the folder structure of the software developed in this thesis and the complete source code of the calculation package, the graphical user interface and the automated tests. The listing is generated from the files on disk when the thesis is built, so it always corresponds to the version of the software described in Chapters 4 and 5. The thesis sources themselves and generated files are omitted.

In total the listed files contain 7075 lines.

## B.1 Folder Structure

```text
ship-support-force/
├── examples/                                             # example cases (placeholder data)
│   ├── example_case.yaml
│   ├── spacing_case.yaml
│   ├── weight_curve_case.yaml
│   ├── welding_case.yaml
│   ├── welding_survey.csv
│   └── welding_survey_case.yaml
├── src/                                                  # source code
│   ├── shipsupport/                                      # calculation package
│   │   ├── solvers/                                      # rigid-hull and FEM solvers, contact iteration
│   │   │   ├── __init__.py
│   │   │   ├── contact.py
│   │   │   ├── fem.py
│   │   │   ├── result.py
│   │   │   └── rigid.py
│   │   ├── __init__.py
│   │   ├── __main__.py
│   │   ├── analysis.py
│   │   ├── cli.py
│   │   ├── deformation.py
│   │   ├── loader.py
│   │   ├── model.py
│   │   ├── report.py
│   │   ├── spacing.py
│   │   ├── strength.py
│   │   ├── units.py
│   │   └── weights.py
│   └── shipsupport_gui/                                  # PyQt5 desktop application
│       ├── __init__.py
│       ├── __main__.py
│       ├── app.py
│       ├── input_panels.py
│       ├── main_window.py
│       ├── plan_panel.py
│       ├── results_panels.py
│       ├── style.py
│       ├── widgets.py
│       └── worker.py
├── tests/                                                # automated test suite
│   ├── helpers.py
│   ├── test_deformation.py
│   ├── test_example_case.py
│   ├── test_fem.py
│   ├── test_gui.py
│   ├── test_loader.py
│   ├── test_rigid.py
│   ├── test_spacing.py
│   ├── test_strength.py
│   └── test_weights.py
├── .gitignore                                            # files not tracked by version control
├── pyproject.toml                                        # package metadata, dependencies, entry points
└── README.md                                             # installation and usage
```

## B.2 Project configuration

### B.2.1 pyproject.toml

34 lines

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "shipsupport"
version = "0.1.0"
description = "Dry-docking block reactions, block spacing design and hull girder strength analysis"
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.24",
    "scipy>=1.10",
    "pyyaml>=6.0",
]

[project.optional-dependencies]
plots = ["matplotlib>=3.7"]
gui = ["matplotlib>=3.7", "PyQt5>=5.15"]
dev = ["pytest>=7", "matplotlib>=3.7", "PyQt5>=5.15"]
thesis = ["matplotlib>=3.7", "PyQt5>=5.15", "pypdfium2>=4"]

[project.scripts]
shipsupport = "shipsupport.cli:main"

[project.gui-scripts]
shipsupport-gui = "shipsupport_gui.app:main"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

## B.3 Calculation package

### B.3.1 src/shipsupport/\_\_init\_\_.py

14 lines

```python
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
```

### B.3.2 src/shipsupport/\_\_main\_\_.py

3 lines

```python
from .cli import main

raise SystemExit(main())
```

### B.3.3 src/shipsupport/units.py

14 lines

```python
"""Unit conversions. Every calculation inside the package is in SI (N, m, Pa)."""

G = 9.80665                      # m/s², standard gravity

TONNE_FORCE = 1000.0 * G         # N per tonne (mass) under standard gravity
MM = 1e-3                        # m per mm
MPA = 1e6                        # Pa per MPa
GPA = 1e9                        # Pa per GPa
KN_PER_MM = 1e3 / MM             # N/m per kN/mm


def to_tonnes(force_n):
    """Force in N -> tonnes (also N·m -> t·m, N/m -> t/m)."""
    return force_n / TONNE_FORCE
```

### B.3.4 src/shipsupport/model.py

461 lines

```python
"""Input data model.

All quantities are SI (N, m, Pa). The longitudinal coordinate x is measured
from the aft perpendicular (AP), positive forward. Vertical displacement is
positive downward, into the blocks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


class InputError(ValueError):
    """Invalid or inconsistent input data."""


def _require_positive(value, what):
    if value is None or not math.isfinite(value) or value <= 0:
        raise InputError(f"{what} must be a positive finite number (got {value!r})")


def _optional_positive(value, what):
    if value is not None:
        _require_positive(value, what)


def _as_array(values, what):
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        raise InputError(f"{what} must be a non-empty list of finite numbers")
    return arr


# --------------------------------------------------------------------------- #
# Weights and hull girder
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class WeightItem:
    """A weight spread as a trapezoid over [xa, xb], or a point load if xa == xb.

    The trapezoid preserves both the weight and its LCG, which is only possible
    with non-negative end intensities when the LCG lies in the middle third.
    """

    name: str
    weight: float   # N (negative allowed for deductions)
    lcg: float      # m from AP
    xa: float       # m
    xb: float       # m

    def __post_init__(self):
        if not math.isfinite(self.weight) or self.weight == 0:
            raise InputError(f"Weight item '{self.name}': weight must be non-zero and finite")
        if self.xb < self.xa:
            raise InputError(f"Weight item '{self.name}': xb ({self.xb}) < xa ({self.xa})")
        if self.is_point:
            if abs(self.lcg - self.xa) > 1e-9:
                raise InputError(f"Weight item '{self.name}': a point load needs lcg == xa == xb")
            return
        third = self.length / 3.0
        xbar = self.lcg - self.xa
        tol = 1e-9 * self.length
        if not (third - tol <= xbar <= 2.0 * third + tol):
            raise InputError(
                f"Weight item '{self.name}': LCG {self.lcg:.3f} m must lie in the middle third of "
                f"its extent [{self.xa:.3f}, {self.xb:.3f}] m for a non-negative trapezoid; "
                f"widen the extent or split the item")

    @property
    def length(self) -> float:
        return self.xb - self.xa

    @property
    def is_point(self) -> bool:
        return self.length <= 1e-9


@dataclass(frozen=True, eq=False)
class WeightCurve:
    """Piecewise-linear weight distribution q(x) ≥ 0 (N/m) through the points (x_k, q_k).

    The load is zero outside the first and last point. x must not decrease; the same x
    twice is a jump of q. Every segment is an exact trapezoid, so the curve is converted
    into one WeightItem per loaded segment (see items()) and the load model is unchanged.
    """

    name: str
    x: np.ndarray        # m from AP
    q: np.ndarray        # N/m

    def __post_init__(self):
        object.__setattr__(self, "x", np.asarray(self.x, dtype=float).ravel())
        object.__setattr__(self, "q", np.asarray(self.q, dtype=float).ravel())
        if self.x.size != self.q.size or self.x.size < 2:
            raise InputError(f"Weight curve '{self.name}': give at least two points (x, q)")
        for k, (x, q) in enumerate(zip(self.x, self.q, strict=True)):
            where = f"Weight curve '{self.name}', point {k + 1} (x = {x:g} m)"
            if not (math.isfinite(x) and math.isfinite(q)):
                raise InputError(f"{where}: x and q must be finite numbers")
            if q < 0:
                raise InputError(f"{where}: q = {q:g} N/m must not be negative")
            if k and x < self.x[k - 1]:
                raise InputError(f"{where}: x must not be smaller than at point {k} "
                                 f"(x = {self.x[k - 1]:g} m); give the points from aft to forward")
        if not self.items():
            raise InputError(f"Weight curve '{self.name}' carries no weight")

    def items(self) -> list[WeightItem]:
        """One exact trapezoidal item per segment with load: weight (q_k + q_k+1) L / 2 and
        LCG x_k + L (q_k + 2 q_k+1) / (3 (q_k + q_k+1)). Jumps and unloaded segments are skipped."""
        out = []
        n = self.x.size - 1
        for k in range(n):
            xa, xb, qa, qb = self.x[k], self.x[k + 1], self.q[k], self.q[k + 1]
            length = xb - xa
            if length <= 1e-9 or qa + qb <= 0:
                continue
            lcg = xa + length * (qa + 2.0 * qb) / (3.0 * (qa + qb))
            out.append(WeightItem(f"{self.name} {k + 1}/{n}", 0.5 * (qa + qb) * length,
                                  float(lcg), float(xa), float(xb)))
        return out


@dataclass(frozen=True)
class Section:
    """Hull girder properties, constant over [xa, xb]."""

    xa: float
    xb: float
    EI: float                          # N·m²
    kGA: float | None = None           # N, effective shear stiffness G·As (None = no shear deformation)
    z_deck: float | None = None        # m³
    z_keel: float | None = None        # m³
    shear_area: float | None = None    # m², effective shear area for τ = V / As

    def __post_init__(self):
        if self.xb <= self.xa:
            raise InputError(f"Section [{self.xa}, {self.xb}]: xb must be greater than xa")
        _require_positive(self.EI, f"Section [{self.xa}, {self.xb}] EI")
        _optional_positive(self.kGA, f"Section [{self.xa}, {self.xb}] kGA")
        _optional_positive(self.z_deck, f"Section [{self.xa}, {self.xb}] z_deck")
        _optional_positive(self.z_keel, f"Section [{self.xa}, {self.xb}] z_keel")
        _optional_positive(self.shear_area, f"Section [{self.xa}, {self.xb}] shear_area")


# --------------------------------------------------------------------------- #
# Blocks
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Layer:
    thickness: float   # m
    E: float           # Pa, modulus in the loading direction (perpendicular to grain for timber)

    def __post_init__(self):
        _require_positive(self.thickness, "Layer thickness")
        _require_positive(self.E, "Layer E")


@dataclass(frozen=True)
class BlockType:
    """A block build-up. Layers act as springs in series over the contact area.

    If crush_stress is given the block is bilinear: stiffness k up to the crush
    force, then post_crush_ratio * k.
    """

    name: str
    contact_width: float                     # m, transverse
    contact_length: float                    # m, longitudinal
    layers: tuple[Layer, ...] = ()
    stiffness_override: float | None = None  # N/m, replaces the layer calculation
    capacity: float | None = None            # N, rated block load
    allowable_pressure: float | None = None  # Pa, allowable bearing pressure
    crush_stress: float | None = None        # Pa
    post_crush_ratio: float = 0.05

    def __post_init__(self):
        _require_positive(self.contact_width, f"Block type '{self.name}' contact_width")
        _require_positive(self.contact_length, f"Block type '{self.name}' contact_length")
        if not self.layers and self.stiffness_override is None:
            raise InputError(f"Block type '{self.name}': give layers or a stiffness override")
        _optional_positive(self.stiffness_override, f"Block type '{self.name}' stiffness")
        _optional_positive(self.capacity, f"Block type '{self.name}' capacity")
        _optional_positive(self.allowable_pressure, f"Block type '{self.name}' allowable_pressure")
        _optional_positive(self.crush_stress, f"Block type '{self.name}' crush_stress")
        if not 0.0 <= self.post_crush_ratio <= 1.0:
            raise InputError(f"Block type '{self.name}': post_crush_ratio must be in [0, 1]")

    @property
    def area(self) -> float:
        return self.contact_width * self.contact_length

    @property
    def stiffness(self) -> float:
        """Initial vertical stiffness, N/m."""
        if self.stiffness_override is not None:
            return self.stiffness_override
        compliance = sum(layer.thickness / layer.E for layer in self.layers)
        return self.area / compliance

    @property
    def crush_force(self) -> float | None:
        return None if self.crush_stress is None else self.crush_stress * self.area


@dataclass(frozen=True)
class Block:
    x: float                 # m from AP
    type: BlockType
    gap: float = 0.0         # m, positive = block top below the keel reference line
    label: str = ""
    strength: float | None = None   # N, force this block can carry; None = its type's capacity

    def __post_init__(self):
        _optional_positive(self.strength, f"Block at x = {self.x} strength")

    @property
    def effective_strength(self) -> float | None:
        """Force the block can carry (N): its own strength, else its type's capacity."""
        return self.strength if self.strength is not None else self.type.capacity


# --------------------------------------------------------------------------- #
# Permissible hull girder limits
# --------------------------------------------------------------------------- #

@dataclass(frozen=True, eq=False)
class BendingEnvelope:
    """Permissible still-water bending moment; sag and hog are positive magnitudes (N·m)."""

    x: np.ndarray
    sag: np.ndarray
    hog: np.ndarray

    def __post_init__(self):
        for name in ("x", "sag", "hog"):
            object.__setattr__(self, name, _as_array(getattr(self, name), f"Bending envelope {name}"))
        if not (self.x.size == self.sag.size == self.hog.size):
            raise InputError("Bending envelope: x, sag and hog must have equal length")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Bending envelope: x must be strictly increasing")
        if np.any(self.sag <= 0) or np.any(self.hog <= 0):
            raise InputError("Bending envelope: limits must be positive magnitudes")

    def at(self, x):
        return np.interp(x, self.x, self.sag), np.interp(x, self.x, self.hog)


@dataclass(frozen=True, eq=False)
class ShearEnvelope:
    """Permissible still-water shear force magnitude (N)."""

    x: np.ndarray
    limit: np.ndarray

    def __post_init__(self):
        for name in ("x", "limit"):
            object.__setattr__(self, name, _as_array(getattr(self, name), f"Shear envelope {name}"))
        if self.x.size != self.limit.size:
            raise InputError("Shear envelope: x and limit must have equal length")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Shear envelope: x must be strictly increasing")
        if np.any(self.limit <= 0):
            raise InputError("Shear envelope: limits must be positive")

    def at(self, x):
        return np.interp(x, self.x, self.limit)


# --------------------------------------------------------------------------- #
# Welding deformation of the assembled hull
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class WeldJoint:
    """Angular distortion (kink) at an erection joint between two hull blocks."""

    x: float        # m from AP
    angle: float    # rad; positive = concave-up kink, which lifts both hull ends

    def __post_init__(self):
        if not (math.isfinite(self.x) and math.isfinite(self.angle)):
            raise InputError("Weld joint x and angle must be finite")


@dataclass(frozen=True)
class ShrinkageZone:
    """Uniform welding-shrinkage curvature over a hull block [xa, xb]."""

    xa: float
    xb: float
    curvature: float   # 1/m; positive = concave up (deck shortens more than bottom), lifts the ends

    def __post_init__(self):
        if self.xb <= self.xa:
            raise InputError(f"Shrinkage zone [{self.xa}, {self.xb}]: xb must be greater than xa")
        if not math.isfinite(self.curvature):
            raise InputError("Shrinkage zone curvature must be finite")


PROFILE_INTERPOLATIONS = ("pchip", "linear", "spline")


@dataclass(frozen=True, eq=False)
class MeasuredProfile:
    """Surveyed keel line after welding: upward lift (m) at stations x (m)."""

    x: np.ndarray
    lift: np.ndarray
    interpolation: str = "pchip"
    identify_joints: bool = False   # fit joint angles to the survey instead of interpolating it
    regularization: float = 0.0     # Tikhonov weight for the joint-angle fit

    def __post_init__(self):
        object.__setattr__(self, "x", _as_array(self.x, "Measured profile x"))
        object.__setattr__(self, "lift", _as_array(self.lift, "Measured profile lift"))
        if self.x.size != self.lift.size or self.x.size < 2:
            raise InputError("Measured profile: x and lift need the same length (at least 2 points)")
        if np.any(np.diff(self.x) <= 0):
            raise InputError("Measured profile: x must be strictly increasing")
        if self.interpolation not in PROFILE_INTERPOLATIONS:
            raise InputError(f"Measured profile: interpolation must be one of {PROFILE_INTERPOLATIONS}")
        if self.regularization < 0:
            raise InputError("Measured profile: regularization must be >= 0")


@dataclass(eq=False)
class WeldingDeformation:
    """Initial (stress-free) keel shape of the welded hull.

    The shape is the sum of joint kinks, zone curvatures and an interpolated survey.
    Only its deviation from a straight line matters: a straight-line part is a
    rigid-body tilt that the support solvers absorb.
    """

    joints: tuple[WeldJoint, ...] = ()
    zones: tuple[ShrinkageZone, ...] = ()
    measured: MeasuredProfile | None = None

    def __post_init__(self):
        self.joints = tuple(self.joints)
        self.zones = tuple(self.zones)
        if not self.joints and not self.zones and self.measured is None:
            raise InputError("Welding deformation needs joints, zones or a measured profile")
        if self.measured is not None and self.measured.identify_joints and not self.joints:
            raise InputError("Identifying joint angles from the survey needs the joint positions")


# --------------------------------------------------------------------------- #
# Ship, settings, case
# --------------------------------------------------------------------------- #

@dataclass(eq=False)
class Ship:
    name: str
    lpp: float
    weights: list[WeightItem]
    sections: list[Section]
    x_aft: float = 0.0                     # m, aft end of the hull
    x_fwd: float | None = None             # m, forward end of the hull (default Lpp)
    declared_weight: float | None = None   # N, from the docking plan, for cross-checking
    declared_lcg: float | None = None      # m
    structure_positions: np.ndarray | None = None   # m, floors / bulkheads / web frames
    permissible_bm: BendingEnvelope | None = None
    permissible_sf: ShearEnvelope | None = None
    # the weight curve of the input, if any; its segments are already among the weights
    weight_curve: WeightCurve | None = None

    def __post_init__(self):
        _require_positive(self.lpp, "Ship Lpp")
        if self.x_fwd is None:
            self.x_fwd = self.lpp
        if self.x_fwd <= self.x_aft:
            raise InputError("Ship x_fwd must be greater than x_aft")
        if not self.weights:
            raise InputError("Ship needs at least one weight item")
        if not self.sections:
            raise InputError("Ship needs at least one hull section")
        if self.total_weight <= 0:
            raise InputError("Ship total weight must be positive")
        if self.structure_positions is not None:
            self.structure_positions = np.sort(_as_array(self.structure_positions, "Structure positions"))

    @property
    def total_weight(self) -> float:
        return float(sum(w.weight for w in self.weights))

    @property
    def lcg(self) -> float:
        return float(sum(w.weight * w.lcg for w in self.weights)) / self.total_weight


@dataclass
class Settings:
    methods: tuple[str, ...] = ("rigid", "fem")
    max_element_length: float = 0.5               # m
    shear_deformation: bool = True                 # Timoshenko elements where kGA is given
    allowable_bending_stress: float | None = None  # Pa
    allowable_shear_stress: float | None = None    # Pa
    structure_tolerance: float = 0.15              # m, allowed block offset from structure
    weight_tolerance: float = 0.005                # relative, declared vs computed weight
    max_iterations: int = 200

    def __post_init__(self):
        _require_positive(self.max_element_length, "max_element_length")
        _optional_positive(self.allowable_bending_stress, "allowable_bending_stress")
        _optional_positive(self.allowable_shear_stress, "allowable_shear_stress")
        if self.structure_tolerance < 0:
            raise InputError("structure_tolerance must be >= 0")
        if self.weight_tolerance < 0:
            raise InputError("weight_tolerance must be >= 0")
        if self.max_iterations < 1:
            raise InputError("max_iterations must be >= 1")


@dataclass
class SpacingDesign:
    """Spacing design: the positions of a given number of blocks with given strengths
    (see spacing.py). The block count and the order of the blocks are kept."""

    method: str = "fem"                      # support model: fem | rigid
    target_utilisation: float = 1.0          # required max of reaction / strength
    min_spacing: float = 1.0                 # m
    max_spacing: float | None = None         # m, None = 2 × mean spacing of the given line
    x_min: float | None = None               # m, aft limit of the block line (None = first block)
    x_max: float | None = None               # m, forward limit (None = last block)
    fix_end_blocks: bool = False             # keep the first and last block where they are
    snap_to_structure: bool = True           # final positions at structure positions
    grid: float | None = None                # m, rounding step when not snapping (None = none)
    max_iterations: int = 20                 # sequential-LP iterations

    def __post_init__(self):
        if self.method not in ("rigid", "fem"):
            raise InputError("Spacing design method must be 'rigid' or 'fem'")
        _require_positive(self.target_utilisation, "Spacing design target_utilisation")
        _require_positive(self.min_spacing, "Spacing design min_spacing")
        if self.max_spacing is not None and self.max_spacing < self.min_spacing:
            raise InputError("Spacing design max_spacing must not be smaller than min_spacing")
        if self.x_min is not None and self.x_max is not None and self.x_max <= self.x_min:
            raise InputError("Spacing design x_max must be greater than x_min")
        _optional_positive(self.grid, "Spacing design grid")
        if self.max_iterations < 0:
            raise InputError("Spacing design max_iterations must be >= 0")


@dataclass(eq=False)
class Case:
    ship: Ship
    blocks: list[Block]
    settings: Settings = field(default_factory=Settings)
    deformation: WeldingDeformation | None = None
    spacing: SpacingDesign | None = None

    def __post_init__(self):
        if len(self.blocks) < 2:
            raise InputError("At least two blocks are required")
```

### B.3.5 src/shipsupport/loader.py

420 lines

```python
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
```

### B.3.6 src/shipsupport/weights.py

118 lines

```python
"""Longitudinal weight distribution.

Weight items are spread as trapezoids onto a mesh whose nodes include every
item boundary, so each element carries an exactly linear load and all
integrals are exact.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from .model import InputError, WeightItem

NODE_TOL = 1e-6  # m, points closer than this are the same node


def trapezoid_ends(item: WeightItem) -> tuple[float, float]:
    """End intensities (N/m) of the trapezoid that preserves the item's weight and LCG.

    From W = L(wa + wb)/2 and W·x̄ = L²(wa + 2wb)/6.
    """
    L = item.length
    xbar = item.lcg - item.xa
    wa = 2.0 * item.weight / L * (2.0 - 3.0 * xbar / L)
    wb = 2.0 * item.weight / L * (3.0 * xbar / L - 1.0)
    return wa, wb


def build_nodes(key_points: Iterable[float], x_aft: float, x_fwd: float,
                max_element_length: float) -> np.ndarray:
    """Sorted mesh nodes containing every key point, no element longer than max_element_length."""
    if max_element_length <= 0:
        raise InputError("max_element_length must be positive")
    pts = np.asarray(list(key_points) + [x_aft, x_fwd], dtype=float)
    outside = (pts < x_aft - NODE_TOL) | (pts > x_fwd + NODE_TOL)
    if outside.any():
        raise InputError(f"x = {pts[outside][0]:.3f} m lies outside the hull extent "
                         f"[{x_aft:.3f}, {x_fwd:.3f}] m")
    pts = np.sort(np.clip(pts, x_aft, x_fwd))
    pts = pts[np.concatenate([[True], np.diff(pts) > NODE_TOL])]
    pts[0], pts[-1] = x_aft, x_fwd

    nodes = [pts[:1]]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(1, math.ceil((b - a) / max_element_length - 1e-9))
        nodes.append(a + (b - a) * np.arange(1, n + 1) / n)
    return np.concatenate(nodes)


def node_index(nodes: np.ndarray, x) -> np.ndarray:
    """Index of the node at each x; raises if x is not a mesh node."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    idx = np.clip(np.searchsorted(nodes, x), 1, len(nodes) - 1)
    idx = np.where(np.abs(x - nodes[idx - 1]) <= np.abs(nodes[idx] - x), idx - 1, idx)
    off = np.abs(nodes[idx] - x) > NODE_TOL
    if off.any():
        raise InputError(f"x = {x[off][0]:.4f} m is not a mesh node")
    return idx


@dataclass(frozen=True, eq=False)
class LoadModel:
    """Piecewise-linear distributed load plus nodal point loads (downward positive)."""

    nodes: np.ndarray
    q1: np.ndarray      # N/m at the aft end of each element
    q2: np.ndarray      # N/m at the forward end of each element
    point: np.ndarray   # N at each node

    @property
    def lengths(self) -> np.ndarray:
        return np.diff(self.nodes)

    @property
    def total(self) -> float:
        return float(np.sum(self.lengths * (self.q1 + self.q2) / 2.0) + self.point.sum())

    @property
    def moment(self) -> float:
        """First moment of the load about x = 0, N·m."""
        x1, x2, L = self.nodes[:-1], self.nodes[1:], self.lengths
        distributed = np.sum(L * (self.q1 * (2.0 * x1 + x2) + self.q2 * (x1 + 2.0 * x2)) / 6.0)
        return float(distributed + np.dot(self.point, self.nodes))

    @property
    def lcg(self) -> float:
        return self.moment / self.total

    def intensity_curve(self) -> tuple[np.ndarray, np.ndarray]:
        """(x, w) pairs per element, suitable for plotting jumps in w(x)."""
        x = np.column_stack([self.nodes[:-1], self.nodes[1:]]).ravel()
        w = np.column_stack([self.q1, self.q2]).ravel()
        return x, w


def build_load_model(items: Sequence[WeightItem], nodes: np.ndarray) -> LoadModel:
    x1, x2 = nodes[:-1], nodes[1:]
    xm = 0.5 * (x1 + x2)
    q1 = np.zeros(xm.size)
    q2 = np.zeros(xm.size)
    point = np.zeros(nodes.size)

    for item in items:
        if item.is_point:
            point[node_index(nodes, item.lcg)[0]] += item.weight
            continue
        node_index(nodes, [item.xa, item.xb])  # elements must not straddle item ends
        wa, wb = trapezoid_ends(item)
        slope = (wb - wa) / item.length
        inside = (xm > item.xa) & (xm < item.xb)
        q1[inside] += wa + slope * (x1[inside] - item.xa)
        q2[inside] += wa + slope * (x2[inside] - item.xa)

    return LoadModel(nodes, q1, q2, point)
```

### B.3.7 src/shipsupport/deformation.py

169 lines

```python
"""Welding deformation of the assembled hull and the resulting block gaps.

A hull built from prefabricated blocks is not straight after welding. Two
mechanisms are modelled, both as an upward keel lift z(x):

* angular distortion φ at an erection joint x_j (a kink: the slope jumps by φ),
      z_j(x) = φ · |x − x_j| / 2
* a uniform curvature κ over a hull block [xa, xb] from differential shrinkage
  of deck and bottom (z'' = κ inside the block, straight outside),
      z(x) = κ · ψ(x),  ψ = 0 (x < xa), (x − xa)²/2 (xa ≤ x ≤ xb),
                            (xb − xa)(x − (xa + xb)/2) (x > xb)

Positive φ and κ make the keel concave up, i.e. they lift the hull ends.
Alternatively a surveyed keel line is interpolated, or used to identify the
joint angles by regularised linear least squares.

Only the deviation from a straight line matters (a straight line is a rigid
tilt absorbed by the support solvers), so the gap at block i is the lift
relative to the least-squares line through the block positions, shifted so
that the lowest block point has zero gap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from .model import MeasuredProfile, WeldingDeformation, WeldJoint


def kink_shape(x: np.ndarray, xj: float) -> np.ndarray:
    """Lift per radian of a joint kink at xj (symmetric form)."""
    return 0.5 * np.abs(np.asarray(x, dtype=float) - xj)


def zone_shape(x: np.ndarray, xa: float, xb: float) -> np.ndarray:
    """Lift per unit curvature of a zone [xa, xb]: z'' = 1 inside, 0 outside."""
    x = np.asarray(x, dtype=float)
    inside = (x - xa) ** 2 / 2.0
    after = (xb - xa) * (x - 0.5 * (xa + xb))
    return np.where(x < xa, 0.0, np.where(x <= xb, inside, after))


def interpolate_profile(profile: MeasuredProfile, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if profile.interpolation == "linear":
        # linear extrapolation with the end slopes
        y = np.interp(x, profile.x, profile.lift)
        lo, hi = x < profile.x[0], x > profile.x[-1]
        s0 = (profile.lift[1] - profile.lift[0]) / (profile.x[1] - profile.x[0])
        s1 = (profile.lift[-1] - profile.lift[-2]) / (profile.x[-1] - profile.x[-2])
        y[lo] = profile.lift[0] + s0 * (x[lo] - profile.x[0])
        y[hi] = profile.lift[-1] + s1 * (x[hi] - profile.x[-1])
        return y
    if profile.interpolation == "spline" and profile.x.size >= 3:
        return CubicSpline(profile.x, profile.lift, bc_type="natural", extrapolate=True)(x)
    return PchipInterpolator(profile.x, profile.lift, extrapolate=True)(x)


@dataclass(frozen=True)
class Identification:
    """Result of fitting joint angles to a surveyed keel line."""

    joints: tuple[WeldJoint, ...]
    rms_residual: float          # m
    max_residual: float          # m


def identify_joint_angles(profile: MeasuredProfile, joint_x, regularization: float = 0.0
                          ) -> Identification:
    """Least-squares fit  z ≈ c0 + c1·x + Σ φ_j |x − x_j|/2  with Tikhonov weight on φ.

    The straight-line terms c0, c1 absorb the unknown survey datum and tilt.
    """
    joint_x = np.asarray(joint_x, dtype=float)
    x, z = profile.x, profile.lift
    scale = max(np.ptp(x), 1.0)
    A = np.column_stack([np.ones_like(x), x / scale] +
                        [kink_shape(x, xj) for xj in joint_x])
    rhs = z.copy()
    if regularization > 0:
        reg = np.zeros((joint_x.size, A.shape[1]))
        reg[:, 2:] = np.sqrt(regularization) * np.eye(joint_x.size) * scale
        A = np.vstack([A, reg])
        rhs = np.concatenate([rhs, np.zeros(joint_x.size)])
    coef, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    fitted = A[:x.size] @ coef
    residual = z - fitted
    joints = tuple(WeldJoint(float(xj), float(phi)) for xj, phi in zip(joint_x, coef[2:]))
    return Identification(joints, float(np.sqrt(np.mean(residual**2))),
                          float(np.max(np.abs(residual))))


def effective_model(deformation: WeldingDeformation) -> tuple[WeldingDeformation, Identification | None]:
    """The deformation model actually used for the gaps.

    * survey with identify_joints: the joint angles are fitted to the survey after
      the known zone-shrinkage part has been subtracted (the survey contains it);
    * survey without identification: the survey is the total keel shape, so any
      joints and zones given as well are not added on top of it;
    * no survey: joints and zones as given.
    """
    m = deformation.measured
    if m is None:
        return deformation, None
    if not m.identify_joints:
        return WeldingDeformation(measured=m), None
    zone_part = sum((z.curvature * zone_shape(m.x, z.xa, z.xb) for z in deformation.zones),
                    np.zeros_like(m.x))
    residual = MeasuredProfile(m.x, m.lift - zone_part, m.interpolation)
    ident = identify_joint_angles(residual, [j.x for j in deformation.joints], m.regularization)
    return WeldingDeformation(joints=ident.joints, zones=deformation.zones), ident


def keel_lift(deformation: WeldingDeformation, x) -> np.ndarray:
    """Upward keel lift z(x) of the welded hull (defined up to a straight line)."""
    x = np.asarray(x, dtype=float)
    z = np.zeros_like(x)
    for j in deformation.joints:
        z += j.angle * kink_shape(x, j.x)
    for zone in deformation.zones:
        z += zone.curvature * zone_shape(x, zone.xa, zone.xb)
    if deformation.measured is not None:
        z += interpolate_profile(deformation.measured, x)
    return z


def gaps_from_lift(lift_at_blocks: np.ndarray, block_x: np.ndarray) -> np.ndarray:
    """Block gaps (m, ≥ 0): lift relative to the least-squares line through the
    block positions, shifted so the lowest point has zero gap."""
    lift = np.asarray(lift_at_blocks, dtype=float)
    x = np.asarray(block_x, dtype=float)
    if np.ptp(x) > 0:
        slope, intercept = np.polyfit(x, lift, 1)
        lift = lift - (intercept + slope * x)
    return lift - lift.min()


def lift_components(deformation: WeldingDeformation, x_curve, block_x) -> dict[str, np.ndarray]:
    """Keel lift split into its parts (joint kinks, zone shrinkage, survey), each with
    the least-squares line through the block positions removed. Their sum plus a
    constant equals detrended_lift."""
    block_x = np.asarray(block_x, dtype=float)
    x_curve = np.asarray(x_curve, dtype=float)
    parts = {}
    if deformation.joints:
        parts["joint kinks"] = WeldingDeformation(joints=deformation.joints)
    if deformation.zones:
        parts["zone shrinkage"] = WeldingDeformation(zones=deformation.zones)
    if deformation.measured is not None:
        parts["survey (interpolated)"] = WeldingDeformation(measured=deformation.measured)
    out = {}
    for name, part in parts.items():
        zb = keel_lift(part, block_x)
        slope, intercept = np.polyfit(block_x, zb, 1) if np.ptp(block_x) > 0 else (0.0, 0.0)
        out[name] = keel_lift(part, x_curve) - (intercept + slope * x_curve)
    return out


def detrended_lift(deformation: WeldingDeformation, x_curve, block_x) -> np.ndarray:
    """The same normalisation as gaps_from_lift, applied to a curve for plotting."""
    block_x = np.asarray(block_x, dtype=float)
    zb = keel_lift(deformation, block_x)
    slope, intercept = np.polyfit(block_x, zb, 1) if np.ptp(block_x) > 0 else (0.0, 0.0)
    shift = (zb - (intercept + slope * block_x)).min()
    x_curve = np.asarray(x_curve, dtype=float)
    return keel_lift(deformation, x_curve) - (intercept + slope * x_curve) - shift
```

### B.3.8 src/shipsupport/solvers/\_\_init\_\_.py

10 lines

```python
from .contact import (CRUSHED, ELASTIC, OPEN, STATE_NAMES, BlockSprings, ContactError,
                      UnstableSupportError)
from .fem import solve_fem
from .result import SupportResult
from .rigid import solve_rigid

__all__ = [
    "BlockSprings", "ContactError", "UnstableSupportError", "SupportResult",
    "solve_fem", "solve_rigid", "OPEN", "ELASTIC", "CRUSHED", "STATE_NAMES",
]
```

### B.3.9 src/shipsupport/solvers/contact.py

169 lines

```python
"""Compression-only (optionally crushing) block springs and the contact-state
iteration shared by the rigid and FEM solvers.

Each block force is piecewise linear in its compression δ = v − gap:

    open     δ ≤ 0          R = 0
    elastic  0 < δ ≤ δy     R = k·δ
    crushed  δ > δy         R = k₂·δ + (k − k₂)·δy

For a fixed set of states the problem is linear, so we iterate on the states
until they are consistent with the solution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from ..model import Block

OPEN, ELASTIC, CRUSHED = 0, 1, 2
STATE_NAMES = {OPEN: "open", ELASTIC: "elastic", CRUSHED: "crushed"}


class ContactError(RuntimeError):
    """The contact iteration failed."""


class UnstableSupportError(ContactError):
    """Too few blocks in contact to support the ship."""


@dataclass(frozen=True, eq=False)
class BlockSprings:
    x: np.ndarray
    k: np.ndarray
    gap: np.ndarray
    crush_deflection: np.ndarray   # δy, inf where the block does not crush
    k_crushed: np.ndarray

    @classmethod
    def from_blocks(cls, blocks: Sequence[Block]) -> "BlockSprings":
        x = np.array([b.x for b in blocks], dtype=float)
        k = np.array([b.type.stiffness for b in blocks], dtype=float)
        gap = np.array([b.gap for b in blocks], dtype=float)
        fy = np.array([np.inf if b.type.crush_force is None else b.type.crush_force for b in blocks])
        k2 = np.array([b.type.post_crush_ratio for b in blocks]) * k
        return cls(x, k, gap, fy / k, k2)

    def __len__(self) -> int:
        return self.x.size

    def states(self, delta: np.ndarray) -> np.ndarray:
        s = np.full(delta.shape, ELASTIC)
        s[delta > self.crush_deflection] = CRUSHED
        s[delta <= 0.0] = OPEN
        return s

    def linearize(self, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(kt, c) such that R = kt·δ + c for the given states."""
        dy = np.where(np.isfinite(self.crush_deflection), self.crush_deflection, 0.0)
        kt = np.select([states == ELASTIC, states == CRUSHED], [self.k, self.k_crushed], 0.0)
        c = np.where(states == CRUSHED, (self.k - self.k_crushed) * dy, 0.0)
        return kt, c

    def force(self, delta: np.ndarray) -> np.ndarray:
        kt, c = self.linearize(self.states(delta))
        return kt * delta + c

    def violation(self, delta: np.ndarray, states: np.ndarray) -> np.ndarray:
        """How far each block's δ lies outside the range of its assumed state."""
        dy = self.crush_deflection
        v = np.zeros(delta.shape)
        v = np.where(states == OPEN, np.maximum(delta, 0.0), v)
        v = np.where(states == ELASTIC, np.maximum(-delta, 0.0) + np.maximum(delta - dy, 0.0), v)
        v = np.where(states == CRUSHED, np.maximum(dy - delta, 0.0), v)
        return v


# linear_solve(kt, c, gap) -> (δ at each block, extra data)
LinearSolve = Callable[[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, object]]


def is_supported(x: np.ndarray, kt: np.ndarray) -> bool:
    stiff = kt > 0.0
    return np.count_nonzero(stiff) >= 2 and np.ptp(x[stiff]) >= 1e-9


def check_support(x: np.ndarray, kt: np.ndarray) -> None:
    if not is_supported(x, kt):
        stiff = kt > 0.0
        raise UnstableSupportError(
            f"Only {np.count_nonzero(stiff)} block position(s) carry load with stiffness; at least "
            f"two separated blocks are needed. Check block positions, gaps and the LCG.")


def _iterate(springs: BlockSprings, linear_solve: LinearSolve, gap: np.ndarray,
             states: np.ndarray, max_iterations: int):
    """State iteration for one gap vector, starting from the given states.

    All states are updated at once; if a set of states repeats (cycling), or a
    full update would leave the hull without support, only the worst violator
    is changed per iteration.
    """
    seen: set[bytes] = set()
    one_at_a_time = False
    for iteration in range(1, max_iterations + 1):
        kt, c = springs.linearize(states)
        check_support(springs.x, kt)
        delta, extra = linear_solve(kt, c, gap)
        new_states = springs.states(delta)
        changed = new_states != states
        if not changed.any():
            return delta, states, iteration, extra

        seen.add(states.tobytes())
        unsupported = not is_supported(springs.x, springs.linearize(new_states)[0])
        if one_at_a_time or unsupported or new_states.tobytes() in seen:
            one_at_a_time = True
            worst = np.where(changed, springs.violation(delta, states), -np.inf)
            i = int(np.argmax(worst))
            new_states = states.copy()
            new_states[i] = springs.states(delta)[i]
        states = new_states
    raise ContactError(f"Block contact iteration did not converge in {max_iterations} iterations")


def solve_contact(springs: BlockSprings, linear_solve: LinearSolve,
                  max_iterations: int = 200) -> tuple[np.ndarray, np.ndarray, int, object]:
    """Iterate block states to consistency; returns (δ, states, iterations, extra).

    The state iteration is tried directly. If it fails — large, uneven gaps such as
    a welding end lift can make it cycle or lose all support on the way — a gap
    continuation is used: the problem is first solved with all gaps equal to the
    smallest gap (a rigid shift, so every block starts in contact), and the gaps
    are then increased in steps to their true values, each step starting from the
    previous contact states. Failed steps are halved.
    """
    initial = np.full(len(springs), ELASTIC)
    try:
        return _iterate(springs, linear_solve, springs.gap, initial, max_iterations)
    except ContactError as direct_error:
        first_error = direct_error

    g0 = float(springs.gap.min())
    target = springs.gap
    try:
        delta, states, total, extra = _iterate(springs, linear_solve, np.full_like(target, g0),
                                               initial, max_iterations)
    except ContactError:
        raise first_error
    lam, step = 0.0, 0.25
    while lam < 1.0:
        nxt = min(1.0, lam + step)
        gap = g0 + nxt * (target - g0)
        try:
            delta, states_next, iterations, extra = _iterate(springs, linear_solve, gap, states,
                                                             max_iterations)
        except ContactError:
            step /= 2.0
            if step < 1.0 / 4096:
                raise first_error
            continue
        total += iterations
        lam, states = nxt, states_next
        step = min(2.0 * step, 0.5)
    return delta, states, total, extra
```

### B.3.10 src/shipsupport/solvers/result.py

42 lines

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .contact import CRUSHED, OPEN


@dataclass(eq=False)
class SupportResult:
    method: str
    block_x: np.ndarray
    block_nodes: np.ndarray
    reactions: np.ndarray     # N, upward on the hull
    delta: np.ndarray         # m, block compression (negative = gap still open)
    states: np.ndarray
    iterations: int
    nodes: np.ndarray
    v: np.ndarray             # m, hull displacement at nodes, downward positive
    theta: np.ndarray         # rad, dv/dx at nodes
    total_load: float         # N
    load_moment: float        # N·m about x = 0

    @property
    def force_residual(self) -> float:
        """(ΣR − W) / W."""
        return (self.reactions.sum() - self.total_load) / self.total_load

    @property
    def moment_residual(self) -> float:
        """(ΣR·x − W·LCG) / (W·span)."""
        span = max(np.ptp(self.nodes), 1e-12)
        return (np.dot(self.reactions, self.block_x) - self.load_moment) / (self.total_load * span)

    @property
    def n_open(self) -> int:
        return int(np.count_nonzero(self.states == OPEN))

    @property
    def n_crushed(self) -> int:
        return int(np.count_nonzero(self.states == CRUSHED))
```

### B.3.11 src/shipsupport/solvers/rigid.py

54 lines

```python
"""Method A: rigid hull on elastic blocks.

The hull displacement is linear, v(x) = a + b·(x − x0), so equilibrium alone
gives a 2×2 system for each set of block states:

    [ Σkt      Σkt·xc  ] [a]   [ W  + Σ(kt·g − c)    ]
    [ Σkt·xc   Σkt·xc² ] [b] = [ M0 + Σ(kt·g − c)·xc ]

with xc = x − x0 (x0 = block centroid, for conditioning) and M0 the load
moment about x0.
"""

from __future__ import annotations

import numpy as np

from ..weights import LoadModel, node_index
from .contact import BlockSprings, UnstableSupportError, solve_contact
from .result import SupportResult


def solve_rigid(load: LoadModel, springs: BlockSprings, max_iterations: int = 200) -> SupportResult:
    x0 = float(springs.x.mean())
    xc = springs.x - x0
    W = load.total
    M0 = load.moment - W * x0

    def linear_solve(kt, c, gap):
        r = kt * gap - c
        A = np.array([[kt.sum(), np.dot(kt, xc)],
                      [np.dot(kt, xc), np.dot(kt, xc * xc)]])
        try:
            a, b = np.linalg.solve(A, [W + r.sum(), M0 + np.dot(r, xc)])
        except np.linalg.LinAlgError as exc:
            raise UnstableSupportError("Rigid-hull system is singular") from exc
        return a + b * xc - gap, (a, b)

    delta, states, iterations, (a, b) = solve_contact(springs, linear_solve, max_iterations)
    kt, c = springs.linearize(states)
    nodes = load.nodes
    return SupportResult(
        method="rigid",
        block_x=springs.x,
        block_nodes=node_index(nodes, springs.x),
        reactions=kt * delta + c,
        delta=delta,
        states=states,
        iterations=iterations,
        nodes=nodes,
        v=a + b * (nodes - x0),
        theta=np.full(nodes.size, b),
        total_load=W,
        load_moment=load.moment,
    )
```

### B.3.12 src/shipsupport/solvers/fem.py

116 lines

```python
"""Method B: elastic hull girder (Euler–Bernoulli / Timoshenko beam elements)
on compression-only block springs.

DOFs per node: v (downward) and θ = dv/dx. The beam is free–free; it is held
only by the springs, so at least two separated blocks must be in contact.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..model import InputError, Section
from ..weights import LoadModel, node_index
from .contact import BlockSprings, UnstableSupportError, solve_contact
from .result import SupportResult


def element_properties(sections: Sequence[Section], nodes: np.ndarray,
                       shear_deformation: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """EI and kGA for each element (kGA = inf means no shear deformation)."""
    xm = 0.5 * (nodes[:-1] + nodes[1:])
    EI = np.full(xm.size, np.nan)
    kGA = np.full(xm.size, np.inf)
    for s in sections:
        m = np.isnan(EI) & (xm >= s.xa) & (xm <= s.xb)
        EI[m] = s.EI
        if shear_deformation and s.kGA is not None:
            kGA[m] = s.kGA
    missing = np.isnan(EI)
    if missing.any():
        raise InputError(f"No hull section covers x = {xm[missing][0]:.3f} m")
    return EI, kGA


def element_stiffness(L: np.ndarray, EI: np.ndarray, kGA: np.ndarray) -> np.ndarray:
    """Timoshenko beam element stiffness, shape (n, 4, 4); Φ = 0 gives Euler–Bernoulli."""
    L = np.asarray(L, dtype=float)
    phi = 12.0 * EI / (kGA * L**2)
    c = EI / (L**3 * (1.0 + phi))
    one = np.ones_like(L)
    k = np.array([
        [12 * one,  6 * L,                  -12 * one,  6 * L],
        [6 * L,     (4 + phi) * L**2,       -6 * L,     (2 - phi) * L**2],
        [-12 * one, -6 * L,                 12 * one,   -6 * L],
        [6 * L,     (2 - phi) * L**2,       -6 * L,     (4 + phi) * L**2],
    ])
    return np.moveaxis(k, -1, 0) * c[:, None, None]


def element_loads(L: np.ndarray, q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Consistent nodal loads for a linear load q1 → q2, shape (n, 4)."""
    return np.column_stack([
        L * (7 * q1 + 3 * q2) / 20,
        L**2 * (3 * q1 + 2 * q2) / 60,
        L * (3 * q1 + 7 * q2) / 20,
        -L**2 * (2 * q1 + 3 * q2) / 60,
    ])


def assemble(load: LoadModel, EI: np.ndarray, kGA: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    nodes = load.nodes
    ndof = 2 * nodes.size
    L = load.lengths
    e = np.arange(L.size)
    dofs = np.column_stack([2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3])

    ke = element_stiffness(L, EI, kGA)
    rows = np.repeat(dofs, 4, axis=1).ravel()
    cols = np.tile(dofs, (1, 4)).ravel()
    K = sp.coo_matrix((ke.ravel(), (rows, cols)), shape=(ndof, ndof)).tocsr()

    F = np.zeros(ndof)
    np.add.at(F, dofs, element_loads(L, load.q1, load.q2))
    F[0::2] += load.point
    return K, F


def solve_fem(load: LoadModel, sections: Sequence[Section], springs: BlockSprings,
              shear_deformation: bool = True, max_iterations: int = 200) -> SupportResult:
    nodes = load.nodes
    ndof = 2 * nodes.size
    EI, kGA = element_properties(sections, nodes, shear_deformation)
    K0, F0 = assemble(load, EI, kGA)
    block_nodes = node_index(nodes, springs.x)
    vdofs = 2 * block_nodes

    def linear_solve(kt, c, gap):
        kdiag = np.zeros(ndof)
        np.add.at(kdiag, vdofs, kt)
        F = F0.copy()
        np.add.at(F, vdofs, kt * gap - c)
        u = spla.spsolve((K0 + sp.diags(kdiag)).tocsc(), F)
        if not np.all(np.isfinite(u)):
            raise UnstableSupportError("FEM system is singular")
        return u[vdofs] - gap, u

    delta, states, iterations, u = solve_contact(springs, linear_solve, max_iterations)
    kt, c = springs.linearize(states)
    return SupportResult(
        method="fem",
        block_x=springs.x,
        block_nodes=block_nodes,
        reactions=kt * delta + c,
        delta=delta,
        states=states,
        iterations=iterations,
        nodes=nodes,
        v=u[0::2],
        theta=u[1::2],
        total_load=load.total,
        load_moment=load.moment,
    )
```

### B.3.13 src/shipsupport/strength.py

205 lines

```python
"""Hull girder shear force, bending moment, stresses and block checks.

Sign conventions:
    V(x) = Σ upward forces − Σ downward loads to the left (aft) of x
    M(x) = ∫ V dx, so sagging is positive and hogging negative
    stresses are positive in tension (sagging puts the deck in compression)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .model import Block, Section, Settings, Ship
from .solvers.contact import STATE_NAMES
from .solvers.result import SupportResult
from .weights import LoadModel


@dataclass(frozen=True)
class Check:
    name: str
    unit: str                  # "N", "N·m" or "Pa"
    demand: float              # governing value (signed)
    capacity: float | None     # limit at the governing location
    utilisation: float         # |demand| / capacity, nan if no limit is defined
    x: float                   # m, governing location

    @property
    def passed(self) -> bool | None:
        return None if np.isnan(self.utilisation) else bool(self.utilisation <= 1.0)


@dataclass(eq=False)
class BlockTable:
    x: np.ndarray
    label: list[str]
    reaction: np.ndarray           # N
    pressure: np.ndarray           # Pa
    pressure_util: np.ndarray      # nan if no allowable
    capacity_util: np.ndarray      # nan if no capacity
    state: list[str]
    structure_offset: np.ndarray   # m to nearest structure, nan if none defined
    off_structure: np.ndarray      # bool


@dataclass(eq=False)
class StrengthResult:
    method: str
    x: np.ndarray                  # sample points; nodes appear twice to show jumps in V
    V: np.ndarray                  # N
    M: np.ndarray                  # N·m
    sigma_deck: np.ndarray         # Pa
    sigma_keel: np.ndarray         # Pa
    tau: np.ndarray                # Pa
    bm_util: np.ndarray
    sf_util: np.ndarray
    closure_shear: float           # V at the forward end / max |V|
    closure_moment: float          # M at the forward end / max |M|
    blocks: BlockTable
    checks: list[Check]

    @property
    def passed(self) -> bool | None:
        results = [c.passed for c in self.checks if c.passed is not None]
        return None if not results else all(results)


def shear_moment_curves(load: LoadModel, block_nodes: np.ndarray, reactions: np.ndarray,
                        samples_per_element: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact V(x) and M(x) for linear element loads, nodal point loads and block reactions."""
    nodes = load.nodes
    nodal = -load.point.astype(float)
    np.add.at(nodal, block_nodes, reactions)
    frac = np.arange(1, samples_per_element) / samples_per_element

    xs, Vs, Ms = [], [], []
    V = M = 0.0
    for j in range(nodes.size):
        xs.append([nodes[j], nodes[j]])
        Vs.append([V, V + nodal[j]])
        Ms.append([M, M])
        V += nodal[j]
        if j == nodes.size - 1:
            break
        L = nodes[j + 1] - nodes[j]
        q1 = load.q1[j]
        dq = load.q2[j] - q1
        s = L * frac
        xs.append(nodes[j] + s)
        Vs.append(V - (q1 * s + dq * s**2 / (2 * L)))
        Ms.append(M + V * s - (q1 * s**2 / 2 + dq * s**3 / (6 * L)))
        M += V * L - (q1 * L**2 / 2 + dq * L**2 / 6)
        V -= q1 * L + dq * L / 2
    return np.concatenate(xs), np.concatenate(Vs), np.concatenate(Ms)


def section_property(sections: Sequence[Section], x: np.ndarray, attr: str) -> np.ndarray:
    out = np.full(x.shape, np.nan)
    for s in sections:
        value = getattr(s, attr)
        if value is None:
            continue
        m = np.isnan(out) & (x >= s.xa - 1e-9) & (x <= s.xb + 1e-9)
        out[m] = value
    return out


def governing(name: str, unit: str, x, demand, capacity) -> Check | None:
    """The location with the highest utilisation, or the largest |demand| if no limits exist."""
    x = np.asarray(x, dtype=float)
    demand = np.asarray(demand, dtype=float)
    capacity = np.broadcast_to(np.asarray(capacity, dtype=float), demand.shape)
    finite = np.isfinite(demand)
    if not finite.any():
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        util = np.abs(demand) / capacity
    valid = finite & np.isfinite(util)
    if valid.any():
        i = np.flatnonzero(valid)[np.argmax(util[valid])]
        return Check(name, unit, float(demand[i]), float(capacity[i]), float(util[i]), float(x[i]))
    i = np.flatnonzero(finite)[np.argmax(np.abs(demand[finite]))]
    return Check(name, unit, float(demand[i]), None, float("nan"), float(x[i]))


def block_table(support: SupportResult, blocks: Sequence[Block], ship: Ship,
                settings: Settings) -> BlockTable:
    R = support.reactions
    area = np.array([b.type.area for b in blocks])
    allowable = np.array([np.nan if b.type.allowable_pressure is None else b.type.allowable_pressure
                          for b in blocks])
    capacity = np.array([np.nan if b.effective_strength is None else b.effective_strength for b in blocks])
    pressure = R / area
    if ship.structure_positions is not None:
        offset = np.min(np.abs(support.block_x[:, None] - ship.structure_positions[None, :]), axis=1)
        off = offset > settings.structure_tolerance + 1e-9
    else:
        offset = np.full(R.size, np.nan)
        off = np.zeros(R.size, dtype=bool)
    return BlockTable(
        x=support.block_x,
        label=[b.label for b in blocks],
        reaction=R,
        pressure=pressure,
        pressure_util=pressure / allowable,
        capacity_util=R / capacity,
        state=[STATE_NAMES[int(s)] for s in support.states],
        structure_offset=offset,
        off_structure=off,
    )


def compute_strength(load: LoadModel, support: SupportResult, blocks: Sequence[Block], ship: Ship,
                     settings: Settings, samples_per_element: int = 8) -> StrengthResult:
    x, V, M = shear_moment_curves(load, support.block_nodes, support.reactions, samples_per_element)

    z_deck = section_property(ship.sections, x, "z_deck")
    z_keel = section_property(ship.sections, x, "z_keel")
    shear_area = section_property(ship.sections, x, "shear_area")
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_deck = -M / z_deck
        sigma_keel = M / z_keel
        tau = V / shear_area

    if ship.permissible_bm is not None:
        sag, hog = ship.permissible_bm.at(x)
        bm_cap = np.where(M >= 0.0, sag, hog)
    else:
        bm_cap = np.full(x.shape, np.nan)
    sf_cap = ship.permissible_sf.at(x) if ship.permissible_sf is not None else np.full(x.shape, np.nan)
    bm_util = np.abs(M) / bm_cap
    sf_util = np.abs(V) / sf_cap

    table = block_table(support, blocks, ship, settings)
    allow_b = np.nan if settings.allowable_bending_stress is None else settings.allowable_bending_stress
    allow_s = np.nan if settings.allowable_shear_stress is None else settings.allowable_shear_stress
    allowable_p = np.array([np.nan if b.type.allowable_pressure is None else b.type.allowable_pressure
                            for b in blocks])
    capacity = np.array([np.nan if b.effective_strength is None else b.effective_strength for b in blocks])

    checks = [
        governing("Hull bending moment", "N·m", x, M, bm_cap),
        governing("Hull shear force", "N", x, V, sf_cap),
        governing("Deck bending stress", "Pa", x, sigma_deck, allow_b),
        governing("Keel bending stress", "Pa", x, sigma_keel, allow_b),
        governing("Hull shear stress", "Pa", x, tau, allow_s),
        governing("Block bearing pressure", "Pa", table.x, table.pressure, allowable_p),
        governing("Block load vs capacity", "N", table.x, table.reaction, capacity),
    ]

    v_max = max(np.max(np.abs(V)), 1e-12)
    m_max = max(np.max(np.abs(M)), 1e-12)
    return StrengthResult(
        method=support.method,
        x=x, V=V, M=M,
        sigma_deck=sigma_deck, sigma_keel=sigma_keel, tau=tau,
        bm_util=bm_util, sf_util=sf_util,
        closure_shear=float(V[-1] / v_max),
        closure_moment=float(M[-1] / m_max),
        blocks=table,
        checks=[c for c in checks if c is not None],
    )
```

### B.3.14 src/shipsupport/spacing.py

509 lines

```python
"""Spacing design: the positions of a given number of blocks with given strengths.

Problem
    N blocks in a fixed order, block i with strength S_i. Find the positions
    x_1 < … < x_N within the allowed extent [x_lo, x_hi], with
    s_min ≤ x_{i+1} − x_i ≤ s_max, that minimise the largest utilisation

        U(x) = max_i R_i(x) / S_i,

    where R_i(x) are the block reactions of the support model. The layout fits when
    U ≤ u_target. The block count and the order of the blocks (and so which block has
    which strength) are not changed; only the spacing is designed.

Method (three phases, every candidate is analysed with the full support model)
    1. Strength-weighted equidistribution, a fixed-point iteration: from the reactions
       of the current layout, a load density along the block line is estimated and the
       blocks are placed so that the load between two neighbours is proportional to
       their mean strength. Fast, global, gives a good start.
    2. Sequential linear programming (SLP) on the minimax problem: the sensitivities
       ∂U_i/∂x_k are computed by finite differences and the linearised problem
           min t  s.t.  U_i + Σ_k J_ik Δ_k ≤ t,  spacing and extent limits,  |Δ_k| ≤ δ
       is solved as a linear programme (HiGHS). A step is accepted only if the analysed
       U decreases; the trust radius δ grows after success and halves after failure.
    3. Rounding to admissible positions (structure positions or a grid) and a 1-opt
       local search that moves the governing block or a neighbour by one grid step.

The analysis is injected as ``evaluate(blocks)`` (dependency injection; in the program it
is analysis.spacing_evaluator, which also recomputes the welding gaps at the candidate
block positions), so the design can be tested with a substitute analysis.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import linprog

from .model import Block, InputError, SpacingDesign

_TOL = 1e-9   # m, positions closer than this are the same position


# --------------------------------------------------------------------------- helpers

def utilisation(reactions: np.ndarray, blocks: Sequence[Block]) -> np.ndarray:
    """R_i / strength_i for every block; NaN where a block has no strength."""
    strength = np.array([np.nan if b.effective_strength is None else b.effective_strength
                         for b in blocks], dtype=float)
    with np.errstate(invalid="ignore"):
        return np.asarray(reactions, dtype=float) / strength


def relabel(blocks: Sequence[Block]) -> list[Block]:
    """Sort by position and label B01, B02, … from aft to forward."""
    ordered = sorted(blocks, key=lambda b: b.x)
    return [dataclasses.replace(b, label=f"B{i + 1:02d}") for i, b in enumerate(ordered)]


def _admissible(y: np.ndarray, s_min: float, s_max: float,
                grid: np.ndarray | None = None) -> np.ndarray | None:
    """Enforce the minimum and maximum spacing (and grid / structure positions) with fixed
    end blocks: one sweep from aft to forward, each block kept in the window that its
    predecessor and the room needed by the remaining blocks leave for it. None if no
    admissible layout is found."""
    n = y.size
    lo, hi = y[0], y[-1]
    if hi - lo < (n - 1) * s_min - _TOL or hi - lo > (n - 1) * s_max + _TOL:
        return None
    y = y.copy()
    for i in range(1, n - 1):
        rest = n - 1 - i                              # gaps still to come
        lower = max(y[i - 1] + s_min, hi - rest * s_max)
        upper = min(y[i - 1] + s_max, hi - rest * s_min)
        if upper < lower - _TOL:
            return None
        if grid is None:
            y[i] = min(max(y[i], lower), upper)
            continue
        ok = grid[(grid >= lower - _TOL) & (grid <= upper + _TOL)]
        if ok.size == 0:
            return None
        y[i] = ok[np.argmin(np.abs(ok - y[i]))]
    gap = hi - y[-2] if n > 2 else hi - lo
    if gap < s_min - _TOL or gap > s_max + _TOL:
        return None
    return y


@dataclass(frozen=True)
class SpacingStep:
    phase: str                 # start | equidistribution | slp | grid
    iteration: int
    max_utilisation: float
    min_utilisation: float
    evaluations: int           # analyses performed so far


@dataclass(eq=False)
class SpacingResult:
    status: str                # fits | does_not_fit
    message: str
    initial_blocks: list[Block]
    final_blocks: list[Block]
    initial_evaluation: object
    final_evaluation: object
    history: list[SpacingStep] = field(default_factory=list)
    evaluations: int = 0
    target: float = 1.0
    suggested_count: int | None = None     # blocks probably needed when it does not fit
    total_strength: float = 0.0            # N, Σ S_i
    total_load: float = 0.0                # N, Σ R_i

    @property
    def fits(self) -> bool:
        return self.status == "fits"

    @property
    def initial_utilisation(self) -> np.ndarray:
        return utilisation(self.initial_evaluation.reactions, self.initial_blocks)

    @property
    def final_utilisation(self) -> np.ndarray:
        return utilisation(self.final_evaluation.reactions, self.final_blocks)

    @property
    def spacing(self) -> np.ndarray:
        """Designed distance from each block to the next (N − 1 values)."""
        return np.diff([b.x for b in self.final_blocks])

    @property
    def initial_spacing(self) -> np.ndarray:
        return np.diff([b.x for b in self.initial_blocks])


Evaluate = Callable[[Sequence[Block]], object]   # returns an object with .reactions


class _Counter:
    """Wraps evaluate(): counts the analyses and returns (evaluation, U_i)."""

    def __init__(self, evaluate: Evaluate, blocks: Sequence[Block], strength: np.ndarray):
        self.evaluate, self.blocks, self.strength = evaluate, list(blocks), strength
        self.count = 0

    def blocks_at(self, x: np.ndarray) -> list[Block]:
        return [dataclasses.replace(b, x=round(float(v), 6)) for b, v in zip(self.blocks, x, strict=True)]

    def __call__(self, x: np.ndarray):
        self.count += 1
        evaluation = self.evaluate(self.blocks_at(x))
        return evaluation, np.maximum(np.asarray(evaluation.reactions, dtype=float), 0.0) / self.strength


def _peak(u: np.ndarray) -> float:
    return float(np.max(u))


def _score(u: np.ndarray) -> tuple[float, float]:
    """Largest utilisation first, then a smooth measure of the high ones (tie-break)."""
    return round(float(np.max(u)), 9), float(np.sum(u ** 8))


# --------------------------------------------------------------------------- phase 1

def equidistribution_step(x: np.ndarray, reactions: np.ndarray, strength: np.ndarray) -> np.ndarray:
    """New positions with the ends fixed: the load between two neighbouring blocks is made
    proportional to their mean strength.

    The reactions spread uniformly over the tributary cells give a load density; its
    cumulative Q(x) is inverted at the cumulative interval weights (S_i + S_{i+1}) / 2.
    """
    edges = np.concatenate([[x[0]], 0.5 * (x[1:] + x[:-1]), [x[-1]]])
    load = np.maximum(reactions, 0.0)
    load = np.maximum(load, 0.05 * max(float(np.mean(load)), 1e-12))     # lifted-off blocks
    Q = np.concatenate([[0.0], np.cumsum(load)])
    w = 0.5 * (strength[:-1] + strength[1:])
    F = np.concatenate([[0.0], np.cumsum(w)]) / np.sum(w)
    y = np.interp(F * Q[-1], Q, edges)
    y[0], y[-1] = x[0], x[-1]
    return y


# --------------------------------------------------------------------------- phase 2

def _sensitivities(run: _Counter, x: np.ndarray, u: np.ndarray, free: np.ndarray,
                   lo: float, hi: float, s_min: float, s_max: float, h: float) -> np.ndarray:
    """J[i, k] = ∂U_i / ∂x_k by one-sided finite differences (admissible direction)."""
    n = x.size
    J = np.zeros((n, n))
    for k in np.flatnonzero(free):
        for step in (h, -h):
            y = x.copy()
            y[k] += step
            gaps = np.diff(y)
            if y[0] < lo - _TOL or y[-1] > hi + _TOL or np.any(gaps < s_min - _TOL) or \
                    np.any(gaps > s_max + _TOL):
                continue
            _, uk = run(y)
            J[:, k] = (uk - u) / step
            break
    return J


def _slp_step(x: np.ndarray, u: np.ndarray, J: np.ndarray, free: np.ndarray, lo: float,
              hi: float, s_min: float, s_max: float, radius: float, smoothing: float):
    """Solve the linearised minimax problem

        min  t + μ Σ_i |x'_{i+2} − 2 x'_{i+1} + x'_i|          (x' = x + Δ)
        s.t. U_i + J_i Δ ≤ t,   s_min ≤ x'_{i+1} − x'_i ≤ s_max,   x'_1 ≥ lo,  x'_N ≤ hi,
             |Δ_k| ≤ radius

    as a linear programme. The small second term (μ = ``smoothing``) prefers a spacing
    that changes gradually among layouts with practically the same largest utilisation.
    Returns (Δx, predicted max U) or None.
    """
    n = x.size
    idx = np.flatnonzero(free)
    m = idx.size
    if m == 0:
        return None
    k2 = max(n - 2, 0)
    # variables: Δ+ (m), Δ− (m), a (n − 2, |second difference|), t;   Δ = Δ+ − Δ−
    nv = 2 * m + k2 + 1
    c = np.zeros(nv)
    c[:2 * m] = 1e-6                                       # no needless moves
    c[2 * m:2 * m + k2] = smoothing
    c[-1] = 1.0
    col = {k: j for j, k in enumerate(idx)}

    def delta(k, sign=1.0):
        r = np.zeros(nv)
        if k in col:
            r[col[k]] += sign
            r[m + col[k]] -= sign
        return r

    rows, rhs = [], []
    Jf = J[:, idx]
    for i in range(n):                                     # U_i + J_i Δ − t ≤ 0
        r = np.zeros(nv)
        r[:m], r[m:2 * m], r[-1] = Jf[i], -Jf[i], -1.0
        rows.append(r)
        rhs.append(-u[i])
    for i in range(n - 1):
        gap = x[i + 1] - x[i]
        d = delta(i + 1) - delta(i)                         # change of gap i
        if d.any():
            rows.append(-d)
            rhs.append(gap - s_min)
            if np.isfinite(s_max):
                rows.append(d)
                rhs.append(s_max - gap)
    for i in range(k2):                                    # ± second difference ≤ a_i
        d2 = x[i + 2] - 2 * x[i + 1] + x[i]
        r = delta(i + 2) - 2 * delta(i + 1) + delta(i)
        a = np.zeros(nv)
        a[2 * m + i] = 1.0
        rows.append(r - a)
        rhs.append(-d2)
        rows.append(-r - a)
        rhs.append(d2)
    ends = (0, n - 1)                                      # only the end blocks touch the extent
    up = [max(min(radius, hi - x[k]) if k in ends else radius, 0.0) for k in idx]
    down = [max(min(radius, x[k] - lo) if k in ends else radius, 0.0) for k in idx]
    bounds = [(0.0, b) for b in up] + [(0.0, b) for b in down] + [(0.0, None)] * k2 + [(None, None)]
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), bounds=bounds, method="highs")
    if res.status != 0:
        return None
    step = np.zeros(n)
    step[idx] = res.x[:m] - res.x[m:2 * m]
    return step, float(res.x[-1])


# --------------------------------------------------------------------------- phase 3

def _grid_neighbour(grid: np.ndarray, x: float, direction: int) -> float | None:
    """The next admissible grid position from x in the given direction."""
    if direction > 0:
        ahead = grid[grid > x + 1e-6]
        return float(ahead[0]) if ahead.size else None
    behind = grid[grid < x - 1e-6]
    return float(behind[-1]) if behind.size else None


def _local_moves(x: np.ndarray, u: np.ndarray, grid: np.ndarray, free: np.ndarray, lo: float,
                 hi: float, s_min: float, s_max: float, hot_tolerance: float = 0.02) -> list:
    """Neighbouring grid layouts for the local search: every governing block (within 2 % of
    the largest utilisation, at most three) and its neighbours are moved by one or two grid
    positions, alone or together with the next block; only admissible layouts are returned."""
    n = x.size
    peak = float(np.max(u))
    hot = [int(i) for i in np.argsort(-u)[:3] if u[i] >= (1.0 - hot_tolerance) * peak]
    movable = sorted({k for i in hot for k in (i - 1, i, i + 1) if 0 <= k < n and free[k]})
    trials, seen = [], set()
    for k in movable:
        for direction in (-1, 1):
            for steps in (1, 2):
                for group in ((k,), (k, k + 1)):
                    if group[-1] >= n or not all(free[g] for g in group):
                        continue
                    z = x.copy()
                    for g in group:
                        v = z[g]
                        for _ in range(steps):
                            v = _grid_neighbour(grid, v, direction)
                            if v is None:
                                break
                        if v is None:
                            break
                        z[g] = v
                    else:
                        gaps = np.diff(z)
                        key = tuple(np.round(z, 6))
                        if key in seen or z[0] < lo - _TOL or z[-1] > hi + _TOL or \
                                np.any(gaps < s_min - _TOL) or np.any(gaps > s_max + _TOL):
                            continue
                        seen.add(key)
                        trials.append(z)
    return trials


def _round_to_grid(x: np.ndarray, grid: np.ndarray, fixed_ends: bool, lo: float, hi: float,
                   s_min: float, s_max: float) -> np.ndarray | None:
    """Snap the free end blocks to the nearest grid position within the extent, then the
    interior blocks in one admissible sweep (see _admissible)."""
    y = x.copy()
    if not fixed_ends:
        inside = grid[(grid >= lo - _TOL) & (grid <= hi + _TOL)]
        if inside.size == 0:
            return None
        for k in (0, -1):
            y[k] = inside[np.argmin(np.abs(inside - y[k]))]
    return _admissible(y, s_min, s_max, grid)


# --------------------------------------------------------------------------- driver

def design_spacing(blocks: Sequence[Block], evaluate: Evaluate, design: SpacingDesign,
                   x_bounds: tuple[float, float], structure: np.ndarray | None = None,
                   equidistribution_iterations: int = 15, fd_step: float = 0.05,
                   local_search_evaluations: int | None = None) -> SpacingResult:
    """Design the positions of the given blocks (count and order kept).

    ``x_bounds`` is the allowed extent of the block line (used unless the design sets
    x_min / x_max); ``structure`` gives admissible final positions if
    ``design.snap_to_structure``; otherwise ``design.grid`` (if set) is used.
    """
    current = relabel(blocks)
    n = len(current)
    strength = np.array([np.nan if b.effective_strength is None else b.effective_strength
                         for b in current], dtype=float)
    if np.any(~np.isfinite(strength)) or np.any(strength <= 0):
        missing = ", ".join(b.label for b, s in zip(current, strength, strict=True)
                            if not np.isfinite(s) or s <= 0)
        raise InputError(f"Spacing design needs a strength for every block; missing: {missing}")
    target = design.target_utilisation
    x0 = np.array([b.x for b in current])
    lo = design.x_min if design.x_min is not None else x_bounds[0]
    hi = design.x_max if design.x_max is not None else x_bounds[1]
    if design.fix_end_blocks:
        lo, hi = x0[0], x0[-1]
    s_min = design.min_spacing
    # default maximum spacing: twice the mean spacing of the given line (an unlimited gap
    # would let the minimax optimum leave long unsupported spans amidships)
    s_max = design.max_spacing if design.max_spacing is not None else \
        max(2.0 * (x0[-1] - x0[0]) / max(n - 1, 1), s_min)
    if (n - 1) * s_min > hi - lo + 1e-6:
        raise InputError(f"Spacing design: {n} blocks at a minimum spacing of {s_min:g} m need "
                         f"{(n - 1) * s_min:.2f} m, but the block line is only {hi - lo:.2f} m "
                         f"long ({lo:g}–{hi:g} m)")
    if np.isfinite(s_max) and (n - 1) * s_max < hi - lo - 1e-6:
        if design.fix_end_blocks:
            raise InputError(f"Spacing design: {n} blocks at a maximum spacing of {s_max:g} m "
                             f"cannot span the fixed block line {lo:g}–{hi:g} m")
        hi = lo + (n - 1) * s_max      # free ends: shorten the line to what the blocks can span
    run = _Counter(evaluate, current, strength)
    history: list[SpacingStep] = []

    def record(phase, it, u):
        history.append(SpacingStep(phase, it, _peak(u), float(np.min(u)), run.count))

    first_eval, u_first = run(x0)
    record("start", 0, u_first)

    # starting layout: ends at the extent limits, interior mapped proportionally
    x = x0.copy()
    if not design.fix_end_blocks:
        x = lo + (x0 - x0[0]) * (hi - lo) / max(x0[-1] - x0[0], _TOL)
        x[0], x[-1] = lo, hi
    x = _admissible(x, s_min, s_max)
    if x is None:
        raise InputError("Spacing design: no admissible start layout within the spacing limits")
    evaluation, u = run(x)
    best = (x, evaluation, u)

    # phase 1 — strength-weighted equidistribution (fixed point, under-relaxed)
    for it in range(1, equidistribution_iterations + 1):
        y = equidistribution_step(x, evaluation.reactions, strength)
        y = _admissible(x + 0.5 * (y - x), s_min, s_max)
        if y is None or np.max(np.abs(y - x)) < 1e-3:
            break
        x = y
        evaluation, u = run(x)
        record("equidistribution", it, u)
        if _score(u) < _score(best[2]):
            best = (x, evaluation, u)
    x, evaluation, u = best

    # phase 2 — sequential linear programming on the minimax problem
    free = np.ones(n, dtype=bool)
    if design.fix_end_blocks:
        free[[0, -1]] = False
    mean_gap = (hi - lo) / max(n - 1, 1)
    radius, radius_max, radius_min = 0.25 * mean_gap, mean_gap, 0.005
    smoothing = 0.002 / mean_gap        # utilisation traded for 1 m of spacing variation
    J = None
    for it in range(1, design.max_iterations + 1):
        if J is None:
            J = _sensitivities(run, x, u, free, lo, hi, s_min, s_max, fd_step)
        step = _slp_step(x, u, J, free, lo, hi, s_min, s_max, radius, smoothing)
        if step is None:
            break
        delta, predicted = step
        if _peak(u) - predicted < 1e-5:
            break                                          # linearised optimum reached
        y = x + delta
        y[0], y[-1] = max(y[0], lo), min(y[-1], hi)
        if np.any(np.diff(y) < s_min - 1e-6):
            y = _admissible(y, s_min, s_max)
            if y is None:
                radius *= 0.5
                continue
        ev, uy = run(y)
        if _score(uy) < _score(u):
            x, evaluation, u = y, ev, uy
            record("slp", it, u)
            radius = min(2.0 * radius, radius_max)
            J = None
        else:
            radius *= 0.5
            if radius < radius_min:
                break

    # phase 3 — admissible positions (structure / grid) and 1-opt local search
    grid = None
    if design.snap_to_structure and structure is not None and np.size(structure):
        grid = np.sort(np.asarray(structure, dtype=float))
    elif design.grid:
        grid = np.arange(np.ceil((lo - _TOL) / design.grid), np.floor((hi + _TOL) / design.grid) + 1) \
            * design.grid
    if grid is not None:
        y = _round_to_grid(x, grid, design.fix_end_blocks, lo, hi, s_min, s_max)
        if y is not None:
            x, (evaluation, u) = y, run(y)
            budget = local_search_evaluations if local_search_evaluations is not None else 25 * n
            start = run.count
            while run.count - start < budget:
                trials = _local_moves(x, u, grid, free, lo, hi, s_min, s_max)
                best_trial = None
                for z in trials:
                    if run.count - start >= budget:
                        break
                    ev, uz = run(z)
                    if best_trial is None or _score(uz) < _score(best_trial[2]):
                        best_trial = (z, ev, uz)
                if best_trial is None or not _score(best_trial[2]) < _score(u):
                    break                                  # local optimum on the grid
                x, evaluation, u = best_trial
            record("grid", 0, u)

    final = run.blocks_at(x)
    peak = _peak(u)
    total_strength = float(np.sum(strength))
    total_load = float(np.sum(np.maximum(evaluation.reactions, 0.0)))
    if peak <= target * (1 + 1e-9):
        status = "fits"
        message = (f"all {n} blocks within strength: max reaction/strength {peak:.3f} "
                   f"(target {target:g}) after {run.count} analyses")
        suggested = None
    else:
        status = "does_not_fit"
        suggested = max(n + 1, int(np.ceil(n * peak / target)))
        if total_strength * target < total_load:
            reason = (f"the total strength Σ S = {total_strength / 9806.65:.1f} t × {target:g} "
                      f"is less than the load {total_load / 9806.65:.1f} t")
        else:
            reason = f"the best spacing gives max reaction/strength {peak:.3f}"
        worst = int(np.argmax(u))
        if worst in (0, n - 1) and not design.fix_end_blocks:
            side = "aft" if worst == 0 else "forward"
            limit = "x_min" if worst == 0 else "x_max"
            reason += (f", governed by the {side} end block at the limit of the block line — "
                       f"allow a longer line ({limit}) to shorten the overhang")
        open_blocks = int(np.count_nonzero(np.asarray(evaluation.reactions) <= 0.0))
        if open_blocks:
            message = (f"{n} blocks with these strengths cannot be spaced within the target: "
                       f"{reason}; {open_blocks} block(s) carry no load because of gaps under "
                       f"the keel (e.g. welding deformation), so the gaps rather than the "
                       f"spacing govern and more blocks may not help — adjust the block "
                       f"heights to the keel line, or use stronger blocks where the keel "
                       f"bears (about {suggested} blocks would be needed if the load spread)")
        else:
            message = (f"{n} blocks with these strengths cannot be spaced within the target: "
                       f"{reason}; or try about {suggested} blocks or stronger blocks")
    return SpacingResult(status, message, current, final, first_eval, evaluation, history,
                         run.count, target, suggested, total_strength, total_load)
```

### B.3.15 src/shipsupport/analysis.py

195 lines

```python
"""Run a complete docking analysis for a Case."""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from .deformation import Identification, effective_model, gaps_from_lift, keel_lift
from .model import Block, Case, InputError, Ship, WeldingDeformation
from .solvers import OPEN, BlockSprings, SupportResult, solve_fem, solve_rigid
from .spacing import SpacingResult, design_spacing
from .strength import StrengthResult, compute_strength
from .units import MM, to_tonnes
from .weights import LoadModel, build_load_model, build_nodes

METHODS = ("rigid", "fem")


@dataclass(eq=False)
class AnalysisResult:
    case: Case
    load: LoadModel
    springs: BlockSprings
    supports: dict[str, SupportResult]
    strength: dict[str, StrengthResult]
    warnings: list[str] = field(default_factory=list)
    weld_gap: np.ndarray | None = None              # m per block, from welding deformation
    deformation: WeldingDeformation | None = None   # model actually used (after identification)
    identification: Identification | None = None
    spacing: SpacingResult | None = None

    @property
    def design_method(self) -> str:
        return "fem" if "fem" in self.supports else next(iter(self.supports))

    @property
    def block_supports(self) -> dict[str, SupportResult]:
        """Support solutions on the given block plan (the designed spacing is reported
        separately)."""
        return {m: r for m, r in self.supports.items() if not m.endswith("+spacing")}


@dataclass(frozen=True, eq=False)
class SpacingEvaluation:
    """One analysed candidate block layout (what the spacing design's evaluate() returns)."""

    load: LoadModel
    support: SupportResult
    springs: BlockSprings          # with the welding gaps at the candidate block positions

    @property
    def reactions(self) -> np.ndarray:
        return self.support.reactions


def _load_model(ship: Ship, blocks: list[Block], max_element_length: float) -> LoadModel:
    """Mesh with a node at every weight-item end, section end and block, and the load on it."""
    keys = ([w.xa for w in ship.weights] + [w.xb for w in ship.weights]
            + [s.xa for s in ship.sections] + [s.xb for s in ship.sections]
            + [b.x for b in blocks])
    return build_load_model(ship.weights,
                            build_nodes(keys, ship.x_aft, ship.x_fwd, max_element_length))


def _springs(blocks: list[Block], deformation: WeldingDeformation | None) -> BlockSprings:
    """Block springs; the welded hull's keel lift at the block positions adds to the gaps."""
    springs = BlockSprings.from_blocks(blocks)
    if deformation is not None:
        springs = dataclasses.replace(
            springs, gap=springs.gap + gaps_from_lift(keel_lift(deformation, springs.x), springs.x))
    return springs


def spacing_evaluator(case: Case, deformation: WeldingDeformation | None, method: str):
    """evaluate(blocks) for the spacing design: mesh, welding gaps at the candidate block
    positions and support solution of a candidate layout, exactly as run_case does for the
    given layout."""
    ship, settings = case.ship, case.settings

    def evaluate(blocks: list[Block]) -> SpacingEvaluation:
        load = _load_model(ship, blocks, settings.max_element_length)
        springs = _springs(blocks, deformation)
        if method == "rigid":
            support = solve_rigid(load, springs, settings.max_iterations)
        else:
            support = solve_fem(load, ship.sections, springs, settings.shear_deformation,
                                settings.max_iterations)
        return SpacingEvaluation(load, support, springs)

    return evaluate


def run_case(case: Case, methods: Iterable[str] | None = None) -> AnalysisResult:
    methods = tuple(methods or case.settings.methods)
    unknown = set(methods) - set(METHODS)
    if not methods or unknown:
        raise InputError(f"Unknown method(s) {sorted(unknown)}; choose from {METHODS}")

    ship, settings = case.ship, case.settings
    # sorted by x; blocks built in code may have no label, so give them B01, B02, … too
    blocks = [b if b.label else dataclasses.replace(b, label=f"B{i + 1:02d}")
              for i, b in enumerate(sorted(case.blocks, key=lambda b: b.x))]
    load = _load_model(ship, blocks, settings.max_element_length)
    warnings = _check_weights(case, load)

    x = np.array([b.x for b in blocks])
    if not x.min() <= load.lcg <= x.max():
        raise InputError(f"LCG {load.lcg:.3f} m lies outside the block line "
                         f"[{x.min():.3f}, {x.max():.3f}] m; the ship would tip")

    # Welding deformation: the welded hull's keel is not straight, which adds to the gaps.
    weld_gap = deformation = identification = None
    if case.deformation is not None:
        given = case.deformation
        if given.measured is not None and not given.measured.identify_joints and \
                (given.joints or given.zones):
            warnings.append("Welding deformation: the surveyed keel line is used as the total "
                            "shape; joints and zones given as well are ignored (set "
                            "identify_joint_angles to fit the joint angles to the survey)")
        deformation, identification = effective_model(given)
        weld_gap = gaps_from_lift(keel_lift(deformation, x), x)
    springs = _springs(blocks, deformation)

    def solve(method: str) -> SupportResult:
        if method == "rigid":
            return solve_rigid(load, springs, settings.max_iterations)
        return solve_fem(load, ship.sections, springs, settings.shear_deformation,
                         settings.max_iterations)

    supports: dict[str, SupportResult] = {m: solve(m) for m in methods}
    strength = {m: compute_strength(load, r, blocks, ship, settings) for m, r in supports.items()}

    for method, r in supports.items():
        if r.n_open:
            lifted = ", ".join(blocks[i].label for i in np.flatnonzero(r.states == OPEN))
            cause = " (welding end lift / open gap)" if weld_gap is not None else \
                " (lift-off or open gap)"
            warnings.append(f"{method}: {r.n_open} block(s) carry no load{cause}: {lifted}")
        if r.n_crushed:
            warnings.append(f"{method}: {r.n_crushed} block(s) loaded beyond crush force")
    if weld_gap is not None:
        warnings.append(f"Welding deformation: keel lift relative to the block line up to "
                        f"{weld_gap.max() / MM:.1f} mm (aft end {weld_gap[0] / MM:.1f} mm, "
                        f"forward end {weld_gap[-1] / MM:.1f} mm)")
    design = "fem" if "fem" in supports else methods[0]
    table = strength[design].blocks
    for i in [i for i, off in enumerate(table.off_structure) if off]:
        warnings.append(f"Block {table.label[i]} at x = {table.x[i]:.3f} m is "
                        f"{table.structure_offset[i]:.3f} m from the nearest structure")
    for method, s in strength.items():
        if abs(s.closure_shear) > 1e-6 or abs(s.closure_moment) > 1e-6:
            warnings.append(f"{method}: shear/moment closure error {s.closure_shear:.2e} / "
                            f"{s.closure_moment:.2e}")

    # Spacing design: positions of the given blocks (count and strengths kept); every
    # candidate is analysed with the welding gaps at its own block positions.
    spacing = None
    if case.spacing is not None:
        sd = case.spacing
        x_lo = max(ship.x_aft, sd.x_min if sd.x_min is not None else blocks[0].x)
        x_hi = min(ship.x_fwd, sd.x_max if sd.x_max is not None else blocks[-1].x)
        spacing = design_spacing(blocks, spacing_evaluator(case, deformation, sd.method),
                                 dataclasses.replace(sd, x_min=x_lo, x_max=x_hi), (x_lo, x_hi),
                                 ship.structure_positions)
        key = f"{sd.method}+spacing"
        supports[key] = dataclasses.replace(spacing.final_evaluation.support, method=key)
        strength[key] = compute_strength(spacing.final_evaluation.load, supports[key],
                                         spacing.final_blocks, ship, settings)
        if not spacing.fits:
            warnings.append(f"Spacing design: {spacing.message}")

    return AnalysisResult(case, load, springs, supports, strength, warnings,
                          weld_gap, deformation, identification, spacing)


def _check_weights(case: Case, load: LoadModel) -> list[str]:
    ship, settings = case.ship, case.settings
    warnings = []
    if abs(load.total - ship.total_weight) > 1e-9 * abs(ship.total_weight) or \
            abs(load.lcg - ship.lcg) > 1e-9 * ship.lpp:
        raise RuntimeError("Internal error: distributed load does not reproduce the weight items")
    if ship.declared_weight is not None:
        diff = (ship.total_weight - ship.declared_weight) / ship.declared_weight
        if abs(diff) > settings.weight_tolerance:
            warnings.append(f"Sum of weight items {to_tonnes(ship.total_weight):.1f} t differs from "
                            f"declared {to_tonnes(ship.declared_weight):.1f} t ({100 * diff:+.2f} %)")
    if ship.declared_lcg is not None:
        diff = ship.lcg - ship.declared_lcg
        if abs(diff) > settings.weight_tolerance * ship.lpp:
            warnings.append(f"LCG of weight items {ship.lcg:.3f} m differs from declared "
                            f"{ship.declared_lcg:.3f} m ({diff:+.3f} m)")
    return warnings
```

### B.3.16 src/shipsupport/report.py

711 lines

```python
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
```

### B.3.17 src/shipsupport/cli.py

71 lines

```python
"""Command line: shipsupport run CASE.yaml [-o OUT] [--methods rigid fem] [--no-plots]"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .analysis import METHODS, run_case
from .loader import load_case
from .model import InputError
from .solvers import ContactError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="shipsupport",
                                     description="Dry-docking block reaction analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="analyse a case file")
    run.add_argument("case", type=Path, help="YAML case file")
    run.add_argument("-o", "--out", type=Path, default=None,
                     help="output directory (default: <case name>_results next to the case)")
    run.add_argument("--methods", nargs="+", choices=METHODS, default=None,
                     help="override the methods in the case file")
    run.add_argument("--no-plots", action="store_true", help="skip PNG plots")
    run.add_argument("--strict", action="store_true",
                     help="exit with code 2 if any check fails")
    args = parser.parse_args(argv)

    from . import report

    # the summary uses characters such as Σ, δ and →; never fail on a narrow console encoding
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")

    try:
        case = load_case(args.case)
        result = run_case(case, args.methods)
    except (InputError, ContactError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    summary = report.format_summary(result)
    print(summary)

    out = args.out or args.case.with_name(f"{args.case.stem}_results")
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.txt").write_text(summary, encoding="utf-8")
    written = [out / "summary.txt", *report.write_csv(result, out)]
    if not args.no_plots:
        try:
            written.append(report.plot_results(result, out / "results.png"))
            written.append(report.plot_block_forces(result, out / "block_forces.png"))
            for extra in (report.plot_welding_data(result, out / "welding_data.png"),
                          report.plot_spacing(result, out / "spacing_design.png")):
                if extra:
                    written.append(extra)
        except ImportError:
            print("note: matplotlib not installed, plots skipped", file=sys.stderr)
    print("\nOutput files:")
    for p in written:
        print(f"  {p}")

    failed = any(s.passed is False for s in result.strength.values()) or \
        (result.spacing is not None and not result.spacing.fits)
    return 2 if args.strict and failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
```

## B.4 Graphical user interface

### B.4.1 src/shipsupport_gui/\_\_init\_\_.py

1 lines

```python
"""PyQt5 desktop front end for the shipsupport analysis package."""
```

### B.4.2 src/shipsupport_gui/\_\_main\_\_.py

3 lines

```python
from .app import main

raise SystemExit(main())
```

### B.4.3 src/shipsupport_gui/app.py

31 lines

```python
"""Entry point: shipsupport-gui [case.yaml]"""

from __future__ import annotations

import sys
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

from .main_window import MainWindow
from .style import apply_style


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv if argv is None else argv
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication.instance() or QApplication(argv)
    app.setApplicationName("Ship Support Force")
    apply_style(app)

    window = MainWindow()
    if len(argv) > 1:
        window.open_path(Path(argv[1]))
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
```

### B.4.4 src/shipsupport_gui/style.py

201 lines

```python
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
```

### B.4.5 src/shipsupport_gui/widgets.py

296 lines

```python
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
```

### B.4.6 src/shipsupport_gui/plan_panel.py

343 lines

```python
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
```

### B.4.7 src/shipsupport_gui/input_panels.py

601 lines

```python
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
```

### B.4.8 src/shipsupport_gui/results_panels.py

252 lines

```python
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
```

### B.4.9 src/shipsupport_gui/worker.py

34 lines

```python
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
```

### B.4.10 src/shipsupport_gui/main_window.py

446 lines

```python
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
```

## B.5 Automated tests

### B.5.1 tests/helpers.py

35 lines

```python
import numpy as np

from shipsupport.model import Block, BlockType, Section, WeightItem
from shipsupport.solvers import BlockSprings
from shipsupport.weights import build_load_model, build_nodes


def block_type(k=1e9, **kw):
    return BlockType("test", contact_width=1.0, contact_length=1.0, stiffness_override=k, **kw)


def make_blocks(xs, k=1e9, gaps=None, **kw):
    t = block_type(k, **kw)
    gaps = [0.0] * len(xs) if gaps is None else gaps
    return [Block(float(x), t, float(g), f"B{i + 1}") for i, (x, g) in enumerate(zip(xs, gaps))]


def uniform_setup(length, q, block_x, k=1e9, EI=1e12, kGA=None, gaps=None, max_len=0.25,
                  extra_items=(), **block_kw):
    """Uniform load q (N/m) over [0, length] on blocks at block_x."""
    items = [WeightItem("uniform", q * length, length / 2, 0.0, length), *extra_items]
    blocks = make_blocks(block_x, k, gaps, **block_kw)
    keys = [0.0, length, *block_x] + [w.xa for w in items] + [w.xb for w in items]
    nodes = build_nodes(keys, 0.0, length, max_len)
    load = build_load_model(items, nodes)
    sections = [Section(0.0, length, EI, kGA)]
    return load, sections, blocks, BlockSprings.from_blocks(blocks)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-30)


def interp_v(result, x):
    return float(np.interp(x, result.nodes, result.v))
```

### B.5.2 tests/test_weights.py

114 lines

```python
import numpy as np
import pytest

from shipsupport.model import InputError, WeightCurve, WeightItem
from shipsupport.weights import build_load_model, build_nodes, node_index, trapezoid_ends


def curve_load(curve, x_aft, x_fwd, max_len=0.5):
    items = curve.items()
    nodes = build_nodes([w.xa for w in items] + [w.xb for w in items], x_aft, x_fwd, max_len)
    return build_load_model(items, nodes)


def test_trapezoid_preserves_weight_and_lcg():
    item = WeightItem("a", 1000.0, 4.5, 2.0, 8.0)
    wa, wb = trapezoid_ends(item)
    L = item.length
    assert L * (wa + wb) / 2 == pytest.approx(1000.0)
    assert item.xa + L * (wa + 2 * wb) / (3 * (wa + wb)) == pytest.approx(4.5)


def test_trapezoid_limits_are_triangles():
    wa, wb = trapezoid_ends(WeightItem("a", 600.0, 2.0, 0.0, 6.0))   # x̄ = L/3
    assert wb == pytest.approx(0.0, abs=1e-12) and wa == pytest.approx(200.0)
    wa, wb = trapezoid_ends(WeightItem("a", 600.0, 4.0, 0.0, 6.0))   # x̄ = 2L/3
    assert wa == pytest.approx(0.0, abs=1e-12) and wb == pytest.approx(200.0)


def test_lcg_outside_middle_third_rejected():
    with pytest.raises(InputError):
        WeightItem("bad", 100.0, 1.0, 0.0, 6.0)


def test_random_items_integrate_exactly():
    rng = np.random.default_rng(1)
    items = []
    for i in range(30):
        xa = rng.uniform(0, 90)
        L = rng.uniform(1, 30)
        lcg = xa + L * rng.uniform(1 / 3, 2 / 3)
        items.append(WeightItem(f"w{i}", rng.uniform(1e4, 1e6), lcg, xa, xa + L))
    items.append(WeightItem("point", 5e5, 33.3, 33.3, 33.3))
    keys = [w.xa for w in items] + [w.xb for w in items]
    nodes = build_nodes(keys, 0.0, 125.0, 0.7)
    load = build_load_model(items, nodes)

    W = sum(w.weight for w in items)
    M = sum(w.weight * w.lcg for w in items)
    assert load.total == pytest.approx(W, rel=1e-12)
    assert load.moment == pytest.approx(M, rel=1e-12)
    assert np.all(np.diff(nodes) <= 0.7 + 1e-12)


def test_nodes_contain_key_points_and_lookup():
    nodes = build_nodes([3.3, 7.25, 7.25], 0.0, 10.0, 1.0)
    assert node_index(nodes, [0.0, 3.3, 7.25, 10.0]).size == 4
    with pytest.raises(InputError):
        node_index(nodes, 3.31)
    with pytest.raises(InputError):
        build_nodes([11.0], 0.0, 10.0, 1.0)


# --------------------------------------------------------------------------- weight curve

def test_curve_segments_are_exact_trapezoids():
    """Every segment [x_k, x_k+1] becomes one item with weight (q_k + q_k+1) L / 2 and
    LCG x_k + L (q_k + 2 q_k+1) / (3 (q_k + q_k+1)); the load model integrates it exactly."""
    x = [0.0, 10.0, 25.0, 40.0]
    q = [5e3, 20e3, 0.0, 12e3]
    curve = WeightCurve("c", x, q)
    items = curve.items()
    assert [w.name for w in items] == ["c 1/3", "c 2/3", "c 3/3"]
    assert items[0].weight == pytest.approx(12.5e4) and items[0].lcg == pytest.approx(
        10.0 * (5e3 + 40e3) / (3 * 25e3))
    assert items[1].lcg == pytest.approx(10.0 + 15.0 / 3)          # triangle, zero at the right
    assert items[2].lcg == pytest.approx(25.0 + 2 * 15.0 / 3)      # triangle, zero at the left
    # closed form of ∫q dx and ∫x q dx for the piecewise-linear curve
    W = sum(0.5 * (q[k] + q[k + 1]) * (x[k + 1] - x[k]) for k in range(3))
    M = sum((x[k + 1] - x[k]) * (q[k] * (2 * x[k] + x[k + 1]) + q[k + 1] * (x[k] + 2 * x[k + 1])) / 6
            for k in range(3))
    load = curve_load(curve, -5.0, 45.0)
    assert load.total == pytest.approx(W, rel=1e-12)
    assert load.moment == pytest.approx(M, rel=1e-12)
    # the intensity of the load model is the curve itself (and zero outside it)
    probe = np.array([-2.0, 3.0, 17.5, 30.0, 39.0, 43.0])
    expected = np.interp(probe, x, q, left=0.0, right=0.0)
    e = np.searchsorted(load.nodes, probe) - 1                      # element of each probe
    s = (probe - load.nodes[e]) / load.lengths[e]
    assert np.allclose(load.q1[e] + s * (load.q2[e] - load.q1[e]), expected, rtol=1e-12,
                       atol=1e-9)


def test_curve_with_a_jump():
    """The same x twice is a jump: 10 t/m on [0, 20], 30 t/m on [20, 30]."""
    curve = WeightCurve("jump", [0.0, 20.0, 20.0, 30.0], [10.0, 10.0, 30.0, 30.0])
    items = curve.items()
    assert len(items) == 2                                       # the zero-length segment is skipped
    load = curve_load(curve, 0.0, 30.0)
    assert load.total == pytest.approx(10 * 20 + 30 * 10, rel=1e-12)
    assert load.lcg == pytest.approx((200 * 10 + 300 * 25) / 500, rel=1e-12)
    i = int(np.searchsorted(load.nodes, 20.0))
    assert load.q2[i - 1] == pytest.approx(10.0) and load.q1[i] == pytest.approx(30.0)


@pytest.mark.parametrize("x, q, message", [
    ([0.0, 10.0, 5.0], [1.0, 1.0, 1.0], "point 3 \\(x = 5 m\\).*must not be smaller"),
    ([0.0, 10.0, 20.0], [1.0, -2.0, 1.0], "point 2 \\(x = 10 m\\).*negative"),
    ([0.0, float("nan")], [1.0, 1.0], "point 2"),
    ([0.0], [1.0], "at least two points"),
    ([0.0, 10.0], [0.0, 0.0], "carries no weight"),
])
def test_invalid_curve_points_are_named(x, q, message):
    with pytest.raises(InputError, match=message):
        WeightCurve("bad", x, q)
```

### B.5.3 tests/test_loader.py

196 lines

```python
"""Reading a case: structure positions, block plan, weight items (LCG or moment) and the
weight curve."""

from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, run_case
from shipsupport.loader import structure_positions
from shipsupport.model import InputError
from shipsupport.units import MM, TONNE_FORCE

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "example_case.yaml"


def example_dict():
    return yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- structure

def test_structure_start_end_spacing():
    assert np.allclose(structure_positions({"start": 0.0, "end": 120.0, "spacing": 0.8}),
                       0.8 * np.arange(151))


def test_structure_start_spacing_count_and_start_end_count():
    x = structure_positions({"start": 6.4, "spacing": 2.4, "count": 44})
    assert x.size == 44 and x[-1] == pytest.approx(109.6)
    x = structure_positions({"start": 10.0, "end": 100.0, "count": 31})
    assert x.size == 31 and np.allclose(np.diff(x), 3.0)


def test_structure_explicit_list_wins():
    assert structure_positions({"x": [3, 1, 2], "start": 0, "count": 9, "spacing": 1}).tolist() \
        == [3, 1, 2]


@pytest.mark.parametrize("spec, message", [
    ({"start": 0, "end": 10}, "exactly two"),
    ({"start": 0, "end": 10, "spacing": 1, "count": 11}, "exactly two"),
    ({"start": 0, "spacing": 1, "count": 2.5}, "whole number"),
    ({"start": 0, "spacing": -1, "count": 3}, "spacing"),
    ({"start": 5, "end": 1, "count": 3}, "end must be"),
])
def test_invalid_structure_specs(spec, message):
    with pytest.raises(InputError, match=message):
        structure_positions(spec)


# --------------------------------------------------------------------------- block plan

def test_block_plan_defines_count_spacing_strength_and_gap_of_every_block():
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "count": 4,
                              "spacing": [2.0, 3.0, 4.0], "strength_t": [150, None, 170, 180],
                              "gap_mm": [0.0, 0.5, 0.0, 1.0]}
    case = case_from_dict(data)
    assert [b.x for b in case.blocks] == [10.0, 12.0, 15.0, 19.0]
    assert [b.label for b in case.blocks] == ["B01", "B02", "B03", "B04"]
    # an empty strength is the type capacity (200 t)
    assert [round(b.effective_strength / TONNE_FORCE) for b in case.blocks] == [150, 200, 170, 180]
    assert [b.gap for b in case.blocks] == pytest.approx([0.0, 0.5 * MM, 0.0, 1.0 * MM])
    data["blocks"]["plan"] = {"type": "keel", "start_x": 6.4, "count": 44, "spacing": 2.4,
                              "strength_t": 190, "gap_mm": 0.2}
    case = case_from_dict(data)
    assert len(case.blocks) == 44 and case.blocks[-1].x == pytest.approx(109.6)
    assert all(b.strength == pytest.approx(190 * TONNE_FORCE) for b in case.blocks)
    assert all(b.gap == pytest.approx(0.2 * MM) for b in case.blocks)


def test_block_count_from_the_lists():
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "spacing": 3.0,
                              "strength_t": [100, 120, 140]}
    assert len(case_from_dict(data).blocks) == 3
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, "spacing": [3.0, 2.0, 1.0, 4.0]}
    assert len(case_from_dict(data).blocks) == 5


@pytest.mark.parametrize("plan, message", [
    ({"count": 4, "spacing": [2.0, 3.0]}, "3 expected"),
    ({"count": 3, "spacing": [2.0, -1.0]}, "> 0"),
    ({"count": 2.5, "spacing": 2.0}, "whole number"),
    ({"count": 1, "spacing": 2.0}, "whole number"),
    ({"count": 3, "spacing": 2.0, "strength_t": [1, 2]}, "3 expected"),
    ({"count": 3, "spacing": 2.0, "gap_mm": [0, 1, 2, 3]}, "3 expected"),
    ({"count": 3, "spacing": "wide"}, "expected a number"),
    ({"spacing": 2.0}, "count"),
    ({"count": 3, "spacing": 2.0, "type": "pier"}, "unknown block type"),
])
def test_invalid_block_plans_are_rejected(plan, message):
    data = example_dict()
    data["blocks"]["plan"] = {"type": "keel", "start_x": 10.0, **plan}
    with pytest.raises(InputError, match=message):
        case_from_dict(data)


def test_block_groups_are_no_longer_accepted():
    data = example_dict()
    data["blocks"].pop("plan")
    data["blocks"]["groups"] = [{"type": "keel", "start": 6.4, "spacing": 2.4, "count": 44}]
    with pytest.raises(InputError, match="blocks.plan"):
        case_from_dict(data)


def test_changing_count_and_spacing_changes_reactions():
    data = example_dict()
    data["structure"] = None
    data["blocks"]["plan"].update(count=33, spacing=3.2)
    result = run_case(case_from_dict(data))
    fem = result.supports["fem"]
    assert fem.reactions.size == 33
    assert np.allclose(np.diff(fem.block_x), 3.2)
    assert fem.reactions.sum() == pytest.approx(result.load.total, rel=1e-12)


# --------------------------------------------------------------------------- weights

def test_moment_gives_the_lcg():
    data = example_dict()
    data["weights"] = [{"name": "a", "weight_t": 400, "moment_tm": 22000, "xa": 40, "xb": 70},
                       {"name": "b", "weight_t": 100, "moment_tm": 3000}]          # point at 30 m
    items = case_from_dict(data).ship.weights
    assert items[0].lcg == pytest.approx(55.0) and (items[0].xa, items[0].xb) == (40.0, 70.0)
    assert items[1].is_point and items[1].lcg == pytest.approx(30.0)
    same = example_dict()
    same["weights"] = [{"name": "a", "weight_t": 400, "lcg": 55.0, "xa": 40, "xb": 70},
                       {"name": "b", "weight_t": 100, "x": 30.0}]
    assert [(w.weight, w.lcg, w.xa, w.xb) for w in case_from_dict(same).ship.weights] == \
        pytest.approx([(w.weight, w.lcg, w.xa, w.xb) for w in items])


def test_moment_and_lcg_must_agree():
    data = example_dict()
    data["weights"][1] = {"name": "main engine", "weight_t": 380, "lcg": 22.0,
                          "moment_tm": 380 * 22.0, "xa": 16.0, "xb": 28.0}
    assert case_from_dict(data).ship.weights[1].lcg == pytest.approx(22.0)
    data["weights"][1]["moment_tm"] = 380 * 23.0
    with pytest.raises(InputError, match="main engine.*does not agree"):
        case_from_dict(data)


@pytest.mark.parametrize("item, message", [
    ({"weight_t": 10, "xa": 0, "xb": 10}, "give lcg, moment_tm"),
    ({"weight_t": 10, "lcg": 5, "xa": 0}, "both xa and xb"),
    ({"weight_t": 10, "x": 5, "xa": 0, "xb": 10}, "no extent"),
    ({"weight_t": 0, "moment_tm": 5}, "must not be zero"),
])
def test_invalid_weight_items(item, message):
    data = example_dict()
    data["weights"] = [item]
    with pytest.raises(InputError, match=message):
        case_from_dict(data)


def test_weight_curve_is_read_as_items():
    data = example_dict()
    data["weight_curve"] = {"name": "lightship", "points": [
        {"x": 0.0, "q_t_per_m": 10.0}, {"x": 60.0, "q_t_per_m": 20.0},
        {"x": 60.0, "q_t_per_m": 30.0}, {"x": 120.0, "q_t_per_m": 0.0}]}
    ship = case_from_dict(data).ship
    curve_items = [w for w in ship.weights if w.name.startswith("lightship")]
    assert [w.name for w in curve_items] == ["lightship 1/3", "lightship 3/3"]
    assert ship.weight_curve.name == "lightship"
    added = (15.0 * 60 + 15.0 * 60) * TONNE_FORCE
    assert ship.total_weight == pytest.approx(4500 * TONNE_FORCE + added)
    as_list = example_dict()
    as_list["weight_curve"] = data["weight_curve"]["points"]          # a plain list of points
    assert case_from_dict(as_list).ship.total_weight == pytest.approx(ship.total_weight)
    only_curve = example_dict()
    only_curve.pop("weights")
    only_curve["weight_curve"] = data["weight_curve"]
    assert case_from_dict(only_curve).ship.total_weight == pytest.approx(added)


@pytest.mark.parametrize("points, message", [
    ([{"x": 0, "q_t_per_m": 1}, {"x": "aft", "q_t_per_m": 1}], r"weight_curve\.points\[1\]\.x"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 10}], r"weight_curve\.points\[1\].*q_t_per_m"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 10, "q_t_per_m": -1}], r"point 2 \(x = 10 m\)"),
    ([{"x": 0, "q_t_per_m": 1}, {"x": 130, "q_t_per_m": 1}], "outside the hull"),
])
def test_invalid_weight_curve_points_are_named(points, message):
    data = example_dict()
    data["weight_curve"] = {"points": points}
    with pytest.raises(InputError, match=message):
        run_case(case_from_dict(data))


def test_no_weights_at_all_is_rejected():
    data = example_dict()
    data.pop("weights")
    with pytest.raises(InputError, match="weight_curve"):
        case_from_dict(data)
```

### B.5.4 tests/test_rigid.py

68 lines

```python
import numpy as np
import pytest

from shipsupport.model import WeightItem
from shipsupport.solvers import OPEN, BlockSprings, UnstableSupportError, solve_rigid

from helpers import make_blocks, uniform_setup


def test_two_blocks_statically_determinate():
    # Uniform 1 kN/m over 20 m, blocks at 2 and 15: R2 = W·(10 − 2)/13
    load, _, _, springs = uniform_setup(20.0, 1e3, [2.0, 15.0], k=3e8)
    r = solve_rigid(load, springs)
    W = 2e4
    assert r.reactions[1] == pytest.approx(W * 8 / 13, rel=1e-12)
    assert r.reactions[0] == pytest.approx(W * 5 / 13, rel=1e-12)


def test_equal_blocks_symmetric_load():
    xs = np.linspace(1, 19, 10)
    load, _, _, springs = uniform_setup(20.0, 5e3, xs)
    r = solve_rigid(load, springs)
    assert np.allclose(r.reactions, 1e5 / 10, rtol=1e-12)


def test_eccentric_closed_form():
    xs = np.linspace(2, 18, 9)
    extra = (WeightItem("point", 4e4, 14.0, 14.0, 14.0),)
    load, _, _, springs = uniform_setup(20.0, 1e3, xs, extra_items=extra)
    r = solve_rigid(load, springs)
    W, xbar = load.total, xs.mean()
    e = load.lcg - xbar
    expected = W / xs.size + W * e * (xs - xbar) / np.sum((xs - xbar) ** 2)
    assert np.all(expected > 0)
    assert np.allclose(r.reactions, expected, rtol=1e-10)


def test_lift_off_matches_reduced_block_set():
    # Heavy point load near the forward end lifts the aft block
    xs = [1.0, 5.0, 9.0, 13.0, 17.0, 19.0]
    extra = (WeightItem("p", 4e5, 18.5, 18.5, 18.5),)
    load, _, _, springs = uniform_setup(20.0, 1e3, xs, extra_items=extra)
    r = solve_rigid(load, springs)
    assert r.states[0] == OPEN and r.reactions[0] == 0.0
    assert np.all(r.reactions >= 0)

    reduced = BlockSprings.from_blocks(make_blocks(xs[1:]))
    r2 = solve_rigid(load, reduced)
    assert np.allclose(r.reactions[1:], r2.reactions, rtol=1e-10)
    assert abs(r.force_residual) < 1e-12 and abs(r.moment_residual) < 1e-12


def test_crushing_two_symmetric_blocks():
    k, fy, ratio = 1e9, 3e5, 0.1
    load, _, _, springs = uniform_setup(10.0, 1e5, [2.0, 8.0], k=k, crush_stress=fy,
                                        post_crush_ratio=ratio)
    r = solve_rigid(load, springs)
    dy = fy / k
    expected_delta = dy + (5e5 - fy) / (ratio * k)
    assert np.allclose(r.delta, expected_delta, rtol=1e-12)
    assert np.allclose(r.reactions, 5e5, rtol=1e-12)
    assert r.n_crushed == 2


def test_lcg_outside_support_is_unstable():
    load, _, _, springs = uniform_setup(20.0, 1e3, [1.0, 4.0])
    with pytest.raises(UnstableSupportError):
        solve_rigid(load, springs)
```

### B.5.5 tests/test_fem.py

118 lines

```python
import numpy as np
import pytest

from shipsupport.model import Section, WeightItem
from shipsupport.solvers import OPEN, BlockSprings, solve_fem, solve_rigid
from shipsupport.weights import build_load_model, build_nodes

from helpers import interp_v, make_blocks, uniform_setup

RIGID_K = 1e15


def test_simply_supported_midspan_deflection_euler_bernoulli():
    L, q, EI = 10.0, 1e3, 1e8
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=RIGID_K, EI=EI, max_len=0.5)
    r = solve_fem(load, sections, springs)
    sag = interp_v(r, L / 2) - 0.5 * (interp_v(r, 0.0) + interp_v(r, L))
    assert sag == pytest.approx(5 * q * L**4 / (384 * EI), rel=1e-9)


def test_simply_supported_midspan_deflection_timoshenko():
    L, q, EI, kGA = 10.0, 1e3, 1e8, 1e6
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=RIGID_K, EI=EI, kGA=kGA,
                                               max_len=0.5)
    r = solve_fem(load, sections, springs, shear_deformation=True)
    sag = interp_v(r, L / 2) - 0.5 * (interp_v(r, 0.0) + interp_v(r, L))
    expected = 5 * q * L**4 / (384 * EI) + q * L**2 / (8 * kGA)
    assert sag == pytest.approx(expected, rel=1e-9)


def test_shear_flag_off_ignores_kGA():
    load, sections, _, springs = uniform_setup(10.0, 1e3, [0.0, 10.0], k=RIGID_K, EI=1e8,
                                               kGA=1e6)
    r = solve_fem(load, sections, springs, shear_deformation=False)
    sag = interp_v(r, 5.0) - 0.5 * (interp_v(r, 0.0) + interp_v(r, 10.0))
    assert sag == pytest.approx(5 * 1e3 * 1e4 / (384 * 1e8), rel=1e-9)


def test_two_span_continuous_beam_reactions():
    L, q = 20.0, 1e3
    load, sections, _, springs = uniform_setup(L, q, [0.0, 10.0, 20.0], k=RIGID_K, EI=1e8)
    r = solve_fem(load, sections, springs)
    W = q * L
    assert np.allclose(r.reactions / W, [3 / 16, 10 / 16, 3 / 16], rtol=1e-6)


def test_stiff_hull_converges_to_rigid_method():
    # The FEM–rigid difference is proportional to 1/EI (until round-off dominates)
    xs = np.arange(3.0, 98.0, 2.5)
    extra = (WeightItem("engine", 2e6, 25.0, 20.0, 30.0), WeightItem("p", 5e5, 80.0, 80.0, 80.0))
    errors = []
    for EI in (1e16, 1e17, 1e18):
        load, sections, _, springs = uniform_setup(100.0, 4e4, list(xs), k=2e9, EI=EI,
                                                   extra_items=extra, max_len=0.5)
        fem = solve_fem(load, sections, springs).reactions
        rigid = solve_rigid(load, springs).reactions
        errors.append(np.max(np.abs(fem - rigid)) / np.max(rigid))
    assert errors[0] / errors[1] == pytest.approx(10.0, rel=0.1)
    assert errors[1] / errors[2] == pytest.approx(10.0, rel=0.1)
    assert errors[2] < 5e-5


def test_flexible_hull_concentrates_load_under_weight():
    xs = np.arange(3.0, 98.0, 2.5)
    extra = (WeightItem("engine", 2e6, 25.0, 20.0, 30.0),)
    load, sections, _, springs = uniform_setup(100.0, 4e4, list(xs), k=2e9, EI=4e12,
                                               extra_items=extra)
    fem = solve_fem(load, sections, springs)
    rigid = solve_rigid(load, springs)
    i = int(np.argmin(np.abs(xs - 25.0)))
    assert fem.reactions[i] > 1.05 * rigid.reactions[i]
    assert abs(fem.force_residual) < 1e-9 and abs(fem.moment_residual) < 1e-9


def test_large_gap_equals_removed_block():
    xs = [2.0, 6.0, 10.0, 14.0, 18.0]
    load, sections, _, springs = uniform_setup(20.0, 1e4, xs, k=1e9, EI=1e9,
                                               gaps=[0, 0, 0.05, 0, 0])
    r = solve_fem(load, sections, springs)
    assert r.states[2] == OPEN and r.reactions[2] == 0.0

    reduced = BlockSprings.from_blocks(make_blocks([2.0, 6.0, 14.0, 18.0], k=1e9))
    r2 = solve_fem(load, sections, reduced)
    assert np.allclose(np.delete(r.reactions, 2), r2.reactions, rtol=1e-9)


def test_small_gap_partially_closes():
    xs = [2.0, 6.0, 10.0, 14.0, 18.0]
    load, sections, _, springs = uniform_setup(20.0, 1e5, xs, k=1e9, EI=1e9,
                                               gaps=[0, 0, 2e-4, 0, 0])
    r = solve_fem(load, sections, springs)
    assert 0.0 < r.reactions[2] < r.reactions[1]
    assert r.delta[2] == pytest.approx(r.reactions[2] / 1e9, rel=1e-9)


def test_mesh_independent_reactions():
    # Hermite elements with consistent loads are nodally exact for piecewise-linear
    # loads, so reactions must not depend on the mesh density.
    xs = list(np.arange(2.0, 59.0, 3.0))
    extra = (WeightItem("w", 1e6, 15.0, 10.0, 20.0),)

    def solve(max_len):
        load, sections, _, springs = uniform_setup(60.0, 2e4, xs, k=1e9, EI=2e11,
                                                   extra_items=extra, max_len=max_len)
        return solve_fem(load, sections, springs).reactions

    reference = solve(0.1)
    for max_len in (3.0, 1.5, 0.5):
        assert np.max(np.abs(solve(max_len) - reference)) < 1e-7 * np.max(reference)


def test_varying_sections_must_cover_hull():
    items = [WeightItem("u", 1e5, 5.0, 0.0, 10.0)]
    nodes = build_nodes([0.0, 10.0, 2.0, 8.0], 0.0, 10.0, 0.5)
    load = build_load_model(items, nodes)
    springs = BlockSprings.from_blocks(make_blocks([2.0, 8.0]))
    with pytest.raises(ValueError):
        solve_fem(load, [Section(0.0, 6.0, 1e9)], springs)
```

### B.5.6 tests/test_strength.py

51 lines

```python
import numpy as np
import pytest

from shipsupport.model import BendingEnvelope, Section, Settings, ShearEnvelope, Ship, WeightItem
from shipsupport.solvers import solve_fem
from shipsupport.strength import compute_strength, shear_moment_curves

from helpers import uniform_setup


def test_simply_supported_shear_and_moment():
    L, q = 10.0, 2e3
    load, sections, _, springs = uniform_setup(L, q, [0.0, L], k=1e15, EI=1e9)
    r = solve_fem(load, sections, springs)
    x, V, M = shear_moment_curves(load, r.block_nodes, r.reactions)
    assert M.max() == pytest.approx(q * L**2 / 8, rel=1e-9)        # sagging positive
    assert x[np.argmax(M)] == pytest.approx(L / 2)
    assert V.max() == pytest.approx(q * L / 2, rel=1e-9)
    assert abs(V[-1]) < 1e-6 and abs(M[-1]) < 1e-6


def test_overhangs_give_hogging_over_supports():
    L, q, a = 10.0, 1e3, 2.0
    load, sections, _, springs = uniform_setup(L, q, [a, L - a], k=1e15, EI=1e9)
    r = solve_fem(load, sections, springs)
    x, V, M = shear_moment_curves(load, r.block_nodes, r.reactions)
    at_support = M[np.argmin(np.abs(x - a))]
    assert at_support == pytest.approx(-q * a**2 / 2, rel=1e-9)
    mid = M[np.argmin(np.abs(x - L / 2))]
    assert mid == pytest.approx(q * L * (L / 2 - a) / 2 - q * (L / 2)**2 / 2, rel=1e-9)


def test_checks_and_utilisation():
    L, q = 10.0, 2e3
    load, _, blocks, springs = uniform_setup(L, q, [0.0, L], k=1e15, EI=1e9)
    sections = [Section(0.0, L, 1e9, None, z_deck=1e-4, z_keel=2e-4, shear_area=1e-3)]
    ship = Ship("t", L, [WeightItem("u", q * L, L / 2, 0.0, L)], sections,
                permissible_bm=BendingEnvelope([0.0, L], [5e4, 5e4], [5e4, 5e4]),
                permissible_sf=ShearEnvelope([0.0, L], [2e4, 2e4]))
    settings = Settings(allowable_bending_stress=200e6, allowable_shear_stress=100e6)
    r = solve_fem(load, sections, springs)
    s = compute_strength(load, r, blocks, ship, settings)
    checks = {c.name: c for c in s.checks}

    bm = checks["Hull bending moment"]
    assert bm.utilisation == pytest.approx(q * L**2 / 8 / 5e4, rel=1e-9)
    deck = checks["Deck bending stress"]
    assert deck.demand == pytest.approx(-(q * L**2 / 8) / 1e-4, rel=1e-9)   # compression
    assert deck.passed is False
    assert abs(checks["Hull shear stress"].demand) == pytest.approx(1e4 / 1e-3, rel=1e-9)
    assert s.passed is False
```

### B.5.7 tests/test_deformation.py

184 lines

```python
from pathlib import Path

import dataclasses

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, run_case
from shipsupport.deformation import (gaps_from_lift, identify_joint_angles, interpolate_profile,
                                     keel_lift, kink_shape, zone_shape)
from shipsupport.model import (InputError, MeasuredProfile, ShrinkageZone, WeldingDeformation,
                               WeldJoint)
from shipsupport.solvers import ELASTIC, OPEN, solve_fem, solve_rigid

from helpers import uniform_setup

ROOT = Path(__file__).resolve().parents[1]


def test_kink_changes_slope_by_angle():
    x = np.array([9.0, 10.0, 11.0])
    z = 2e-3 * kink_shape(x, 10.0)
    slope_left, slope_right = z[1] - z[0], z[2] - z[1]
    assert slope_right - slope_left == pytest.approx(2e-3)


def test_zone_has_unit_curvature_inside_and_is_straight_outside():
    x = np.linspace(-5, 25, 3001)
    z = zone_shape(x, 0.0, 20.0)
    curvature = np.gradient(np.gradient(z, x), x)
    inside = (x > 1) & (x < 19)
    outside = (x < -1) | (x > 21)
    assert np.allclose(curvature[inside], 1.0, atol=1e-6)
    assert np.allclose(curvature[outside], 0.0, atol=1e-6)
    assert zone_shape(np.array([20.0]), 0, 20)[0] == pytest.approx(200.0)


def test_gaps_ignore_straight_line_part():
    x = np.linspace(5, 95, 20)
    lift = 1e-3 * (x - 50) ** 2 / 100
    tilted = lift + 0.4 + 3e-3 * x
    assert np.allclose(gaps_from_lift(lift, x), gaps_from_lift(tilted, x), atol=1e-15)
    assert gaps_from_lift(lift, x).min() == 0.0


def test_reactions_unchanged_by_rigid_tilt_of_keel_line():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    gap = 1e-3 * ((x - 50) / 50) ** 2 * 10          # concave-up keel, 10 mm at the ends
    a = solve_fem(load, sections, dataclasses.replace(springs, gap=gap))
    b = solve_fem(load, sections, dataclasses.replace(springs, gap=gap + 0.01 + 2e-4 * x))
    assert np.allclose(a.reactions, b.reactions, rtol=1e-8, atol=1e-6 * load.total)


def test_large_end_lift_unloads_first_and_last_blocks():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    defo = WeldingDeformation(joints=tuple(WeldJoint(xj, 1.5e-3) for xj in (20, 40, 60, 80)))
    gap = gaps_from_lift(keel_lift(defo, x), x)
    fem = solve_fem(load, sections, dataclasses.replace(springs, gap=gap))
    assert fem.states[0] == OPEN and fem.states[-1] == OPEN
    assert fem.reactions[0] == 0.0 and fem.reactions[-1] == 0.0
    # gap terms k·g are ~40·W here, so round-off in the residual is larger than usual
    assert abs(fem.force_residual) < 1e-7 and abs(fem.moment_residual) < 1e-7
    # complementarity: open blocks have no compression, loaded blocks have positive compression
    assert np.all(fem.delta[fem.states == OPEN] <= 0)
    assert np.all(fem.delta[fem.states == ELASTIC] > 0)


def test_continuation_solves_rigid_case_with_large_gaps():
    xs = list(np.arange(3.0, 98.0, 2.5))
    load, sections, blocks, springs = uniform_setup(100.0, 4e4, xs, k=2e9, EI=4e12)
    x = springs.x
    gap = 0.04 * ((x - 50) / 47) ** 2                # 40 mm end lift: few blocks touch a rigid hull
    rigid = solve_rigid(load, dataclasses.replace(springs, gap=gap))
    assert rigid.n_open > 20
    assert abs(rigid.force_residual) < 1e-12 and abs(rigid.moment_residual) < 1e-12
    assert np.all(rigid.delta[rigid.states == OPEN] <= 0)
    assert np.all(rigid.delta[rigid.states == ELASTIC] > 0)


def test_identification_recovers_joint_angles_despite_datum():
    joints = (20.0, 45.0, 70.0, 95.0)
    angles = np.array([0.8e-3, 0.4e-3, 0.5e-3, 0.7e-3])
    xs = np.linspace(0, 120, 25)
    z = sum(a * kink_shape(xs, xj) for a, xj in zip(angles, joints)) + 0.123 - 4e-4 * xs
    ident = identify_joint_angles(MeasuredProfile(xs, z), joints)
    assert np.allclose([j.angle for j in ident.joints], angles, rtol=1e-9, atol=1e-12)
    assert ident.rms_residual < 1e-12


def test_identification_with_noise_and_regularisation():
    rng = np.random.default_rng(3)
    joints = (20.0, 45.0, 70.0, 95.0)
    angles = np.array([0.8e-3, 0.4e-3, 0.5e-3, 0.7e-3])
    xs = np.linspace(0, 120, 49)
    z = sum(a * kink_shape(xs, xj) for a, xj in zip(angles, joints))
    noisy = z + rng.normal(0, 0.5e-3, xs.size)                 # 0.5 mm survey noise
    ident = identify_joint_angles(MeasuredProfile(xs, noisy), joints, regularization=1e-6)
    assert np.allclose([j.angle for j in ident.joints], angles, atol=0.15e-3)
    assert ident.rms_residual == pytest.approx(0.5e-3, rel=0.35)


@pytest.mark.parametrize("method", ["pchip", "linear", "spline"])
def test_profile_interpolation_passes_through_survey_points(method):
    profile = MeasuredProfile([0, 10, 25, 40, 60], [5e-3, 2e-3, 0.0, 1e-3, 6e-3],
                              interpolation=method)
    assert np.allclose(interpolate_profile(profile, profile.x), profile.lift, atol=1e-15)


def test_loader_converts_shrinkage_to_curvature():
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    case = case_from_dict(data)
    zone = case.deformation.zones[0]
    assert zone.curvature == pytest.approx((0.24e-3 - 0.20e-3) / 10.0)
    assert case.deformation.joints[0].angle == pytest.approx(0.9e-3)
    assert case.spacing is None                                   # spacing_design.enabled: false


def test_invalid_deformation_input():
    with pytest.raises(InputError):
        WeldingDeformation()
    with pytest.raises(InputError):
        ShrinkageZone(10.0, 5.0, 1e-5)
    with pytest.raises(InputError):
        WeldingDeformation(measured=MeasuredProfile([0, 1], [0, 0], identify_joints=True))


def test_survey_csv_and_identification_with_known_zones(tmp_path):
    from shipsupport import load_case
    from shipsupport.deformation import effective_model
    from shipsupport.loader import read_survey_csv

    xs, zs = read_survey_csv(ROOT / "examples" / "welding_survey.csv")
    assert len(xs) == 25 and xs[0] == 0.0 and xs[-1] == 120.0
    case = load_case(ROOT / "examples" / "welding_survey_case.yaml")
    assert case.deformation.measured.identify_joints
    used, ident = effective_model(case.deformation)
    true = [0.9e-3, 0.5e-3, 0.4e-3, 0.5e-3, 0.6e-3]          # used to generate the sample survey
    assert np.allclose([j.angle for j in used.joints], true, atol=0.1e-3)
    assert ident.rms_residual < 0.6e-3                         # survey noise was 0.5 mm
    assert used.zones == case.deformation.zones                # zones kept, not double counted


def test_identification_subtracts_known_zone_part_exactly():
    from shipsupport.deformation import effective_model
    joints = (WeldJoint(30.0, 0.0), WeldJoint(70.0, 0.0))
    zones = (ShrinkageZone(0.0, 50.0, 8e-6), ShrinkageZone(50.0, 100.0, 3e-6))
    xs = np.linspace(0, 100, 21)
    truth = WeldingDeformation(joints=(WeldJoint(30.0, 0.6e-3), WeldJoint(70.0, 0.3e-3)),
                               zones=zones)
    survey = MeasuredProfile(xs, keel_lift(truth, xs) + 0.02, identify_joints=True)
    used, ident = effective_model(WeldingDeformation(joints=joints, zones=zones, measured=survey))
    assert np.allclose([j.angle for j in used.joints], [0.6e-3, 0.3e-3], rtol=1e-8)
    assert ident.rms_residual < 1e-12


def test_survey_without_identification_is_the_total_shape():
    from shipsupport.deformation import effective_model
    survey = MeasuredProfile([0, 50, 100], [5e-3, 0.0, 5e-3])
    used, ident = effective_model(WeldingDeformation(joints=(WeldJoint(50.0, 1e-3),),
                                                     measured=survey))
    assert ident is None and not used.joints and used.measured is survey


def test_bad_survey_file(tmp_path):
    from shipsupport.loader import read_survey_csv
    path = tmp_path / "s.csv"
    path.write_text("station,height\n0,1\n5,2\n", encoding="utf-8")
    with pytest.raises(InputError, match="header"):
        read_survey_csv(path)


def test_welding_example_reproduces_end_lift():
    result = run_case(case_from_dict(yaml.safe_load(
        (ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))))
    fem = result.supports["fem"]
    assert fem.reactions[0] == 0.0 and fem.reactions[-1] == 0.0      # first and last blocks lifted
    assert result.weld_gap[0] > 0.03 and result.weld_gap[-1] > 0.03
    assert np.allclose(result.springs.gap, result.weld_gap)          # no initial block gaps
    assert set(result.supports) == {"rigid", "fem"}
```

### B.5.8 tests/test_spacing.py

239 lines

```python
"""Block plan input and the spacing design for a given number of blocks (spacing.py)."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, load_case, run_case
from shipsupport.model import Block, BlockType, InputError, SpacingDesign
from shipsupport.spacing import design_spacing, equidistribution_step

ROOT = Path(__file__).resolve().parents[1]
TYPE = BlockType("t", 1.0, 1.0, stiffness_override=1e9)


@dataclass
class FakeEvaluation:
    reactions: np.ndarray


def tributary_load(q, a=0.0, b=30.0):
    """Substitute analysis of a hull from a to b: every block carries the load q(x) of its
    tributary cell, and the end blocks also carry the overhangs (the load is conserved)."""
    def evaluate(blocks):
        xs = np.array([b_.x for b_ in blocks])
        edges = np.concatenate([[a], 0.5 * (xs[1:] + xs[:-1]), [b]])
        grid = np.linspace(a, b, 40001)
        mid = 0.5 * (grid[1:] + grid[:-1])
        cell = np.clip(np.searchsorted(edges, mid, side="right") - 1, 0, xs.size - 1)
        return FakeEvaluation(np.bincount(cell, q(mid) * np.diff(grid), minlength=xs.size))
    return evaluate


def blocks_at(xs, strengths):
    return [Block(float(x), TYPE, strength=float(s)) for x, s in zip(xs, strengths, strict=True)]


def free_design(**kw):
    defaults = dict(min_spacing=0.5, max_spacing=10.0, snap_to_structure=False, max_iterations=30)
    defaults.update(kw)
    return SpacingDesign(**defaults)


# --------------------------------------------------------------------------- the method

def test_uniform_load_reaches_the_analytic_optimum():
    """Uniform load q on a hull of length L carried by N equal blocks of strength S: every
    block can carry exactly q L / N, at the cell centres x_i = (i − ½) L / N, so the
    smallest possible largest utilisation is U* = q L / (N S) and all utilisations are equal."""
    n, length, q, s = 11, 30.0, 950.0 / 30.0, 100.0
    xs = np.sort(np.random.default_rng(1).uniform(0.5, 29.5, n))
    result = design_spacing(blocks_at(xs, [s] * n), tributary_load(lambda x: np.full_like(x, q)),
                            free_design(), (0.0, length))
    u_star = q * length / (n * s)
    u = result.final_utilisation
    assert result.fits
    assert np.max(u) == pytest.approx(u_star, rel=0.01)
    assert np.max(u) - np.min(u) < 0.02 * u_star                 # every block equally used
    assert np.max(result.initial_utilisation) > 1.3 * u_star


def test_weak_blocks_get_a_smaller_spacing():
    n = 16
    strengths = [50.0] * 8 + [150.0] * 8                   # weak blocks aft, strong forward
    result = design_spacing(blocks_at(np.linspace(1, 29, n), strengths),
                            tributary_load(lambda x: np.full_like(x, 40.0)),
                            free_design(), (0.0, 30.0))
    gaps = result.spacing
    assert result.fits
    assert np.mean(gaps[1:7]) < 0.6 * np.mean(gaps[8:-1])
    before = result.initial_utilisation
    assert np.max(result.final_utilisation) < 0.75 * np.max(before)
    assert [b.effective_strength for b in result.final_blocks] == strengths   # order kept


def test_equidistribution_step_moves_blocks_towards_a_load_peak():
    xs = np.linspace(0.0, 20.0, 11)
    reactions = np.where(xs < 5.0, 300.0, 50.0)
    y = equidistribution_step(xs, reactions, np.full(11, 100.0))
    assert y[0] == 0.0 and y[-1] == 20.0
    assert np.all(np.diff(y) > 0)
    assert y[2] < xs[2]                                    # denser where the load is high


def test_spacing_limits_extent_and_grid_are_respected():
    q = lambda x: np.where(x < 6.0, 120.0, 15.0)            # noqa: E731
    design = SpacingDesign(min_spacing=0.8, max_spacing=3.0, x_min=-2.0, x_max=34.0,
                           grid=0.4, snap_to_structure=False)
    result = design_spacing(blocks_at(np.linspace(0, 30, 16), [100.0] * 16), tributary_load(q),
                            design, (0.0, 30.0))
    xs = np.array([b.x for b in result.final_blocks])
    gaps = np.diff(xs)
    assert np.all(gaps >= 0.8 - 1e-6) and np.all(gaps <= 3.0 + 1e-6)
    assert xs[0] >= -2.0 - 1e-6 and xs[-1] <= 34.0 + 1e-6
    assert np.allclose(xs / 0.4, np.round(xs / 0.4), atol=1e-6)    # on the 0.4 m grid
    assert [h.phase for h in result.history][0] == "start" and result.history[-1].phase == "grid"


def test_fixed_end_blocks_stay():
    result = design_spacing(blocks_at([1.0, 4.0, 5.0, 9.0, 13.0], [80.0] * 5),
                            tributary_load(lambda x: np.full_like(x, 20.0)),
                            free_design(fix_end_blocks=True), (0.0, 20.0))
    assert result.final_blocks[0].x == 1.0 and result.final_blocks[-1].x == 13.0


def test_reports_when_the_strength_is_not_enough():
    n = 6
    result = design_spacing(blocks_at(np.linspace(0, 20, n), [30.0] * n),
                            tributary_load(lambda x: np.full_like(x, 20.0)),   # 400 > 6 × 30
                            free_design(), (0.0, 20.0))
    assert result.status == "does_not_fit" and not result.fits
    assert result.suggested_count is not None and result.suggested_count > n
    assert "total strength" in result.message


def test_invalid_spacing_problems_are_rejected():
    load = tributary_load(lambda x: np.full_like(x, 1.0))
    with pytest.raises(InputError, match="strength"):
        design_spacing([Block(0.0, TYPE), Block(5.0, TYPE)], load, free_design(), (0.0, 5.0))
    with pytest.raises(InputError, match="minimum spacing"):
        design_spacing(blocks_at(np.linspace(0, 4, 6), [10.0] * 6), load,
                       free_design(min_spacing=1.0), (0.0, 4.0))
    with pytest.raises(InputError):
        SpacingDesign(min_spacing=2.0, max_spacing=1.0)


# --------------------------------------------------------------------------- integration

def test_spacing_example_fits_every_block_within_strength(tmp_path):
    result = run_case(load_case(ROOT / "examples" / "spacing_case.yaml"), ["fem"])
    sp = result.spacing
    assert sp.fits and len(sp.final_blocks) == 44
    assert np.max(sp.initial_utilisation) > 1.0 >= np.max(sp.final_utilisation)
    assert result.strength["fem+spacing"].passed is True
    assert list(result.block_supports) == ["fem"]
    gaps = sp.spacing
    assert gaps.min() >= 0.8 - 1e-6 and gaps.max() <= 4.8 + 1e-6
    structure = np.arange(0.0, 120.01, 0.8)
    for b in sp.final_blocks:                                          # under floors
        assert np.min(np.abs(structure - b.x)) < 1e-6
    assert [b.effective_strength for b in sp.final_blocks] == \
        [b.effective_strength for b in sp.initial_blocks]               # strengths kept in order
    final = result.supports["fem+spacing"]
    assert abs(final.force_residual) < 1e-9 and abs(final.moment_residual) < 1e-9
    from shipsupport.cli import main
    assert main(["run", str(ROOT / "examples" / "spacing_case.yaml"), "-o", str(tmp_path),
                 "--methods", "fem"]) == 0
    for name in ("spacing.csv", "spacing_history.csv", "spacing_design.png"):
        assert (tmp_path / name).exists()
    rows = (tmp_path / "spacing.csv").read_text(encoding="utf-8").splitlines()
    assert len(rows) == 45 and rows[0].startswith("block,x_m,spacing_to_next_m")
    assert "BLOCK SPACING DESIGN" in (tmp_path / "summary.txt").read_text(encoding="utf-8")


def test_initial_gaps_move_with_their_blocks():
    data = yaml.safe_load((ROOT / "examples" / "spacing_case.yaml").read_text(encoding="utf-8"))
    gap = [0.0] * 44
    gap[30] = 0.4
    data["blocks"]["plan"]["gap_mm"] = gap
    sp = run_case(case_from_dict(data), ["fem"]).spacing
    assert sp.final_blocks[30].gap == pytest.approx(0.4e-3)
    assert sp.final_evaluation.springs.gap[30] == pytest.approx(0.4e-3)
    assert np.count_nonzero(sp.final_evaluation.springs.gap) == 1


@pytest.fixture(scope="module")
def welded():
    """The welded hull of welding_case.yaml with the spacing design switched on."""
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    data["spacing_design"]["enabled"] = True
    return run_case(case_from_dict(data))


def test_spacing_design_on_the_welded_hull_recomputes_the_welding_gaps(welded):
    """Welding deformation is a constraint of the spacing design: every candidate layout is
    analysed with the welding gaps at its own block positions."""
    from shipsupport.deformation import gaps_from_lift, keel_lift
    sp = welded.spacing
    x0 = np.array([b.x for b in sp.initial_blocks])
    x1 = np.array([b.x for b in sp.final_blocks])
    assert not np.allclose(x0, x1)
    expected = gaps_from_lift(keel_lift(welded.deformation, x1), x1)
    assert np.allclose(sp.final_evaluation.springs.gap, expected, rtol=0, atol=1e-12)
    assert np.allclose(sp.initial_evaluation.springs.gap, welded.weld_gap, rtol=0, atol=1e-12)
    assert not np.allclose(expected, welded.weld_gap)                  # the gaps moved with x
    # the designed layout analysed from scratch gives the same reactions
    import dataclasses
    fresh = run_case(dataclasses.replace(welded.case, blocks=sp.final_blocks, spacing=None),
                     ["fem"])
    assert np.allclose(fresh.supports["fem"].reactions, welded.supports["fem+spacing"].reactions,
                       rtol=1e-12, atol=1e-6)
    assert np.allclose(fresh.weld_gap, expected, atol=1e-12)


def test_spacing_design_on_the_welded_hull_reports_its_outcome(welded):
    sp = welded.spacing
    u0, u1 = np.max(sp.initial_utilisation), np.max(sp.final_utilisation)
    assert u1 < u0                                                     # better than as built
    assert sp.total_load == pytest.approx(welded.load.total, rel=1e-9)
    if sp.fits:
        assert u1 <= 1.0
    else:
        assert sp.status == "does_not_fit" and sp.suggested_count > len(sp.final_blocks)
        assert any(w.startswith("Spacing design:") for w in welded.warnings)
    from shipsupport.report import format_summary
    text = format_summary(welded)
    assert "recomputed at the designed positions" in text


def scaled_welding(factor):
    data = yaml.safe_load((ROOT / "examples" / "welding_case.yaml").read_text(encoding="utf-8"))
    for joint in data["welding_deformation"]["joints"]:
        joint["angle_mrad"] *= factor
    for zone in data["welding_deformation"]["zones"]:
        zone["deck_shrinkage_mm_per_m"] *= factor
        zone["bottom_shrinkage_mm_per_m"] *= factor
    data["spacing_design"]["enabled"] = True
    return run_case(case_from_dict(data), ["fem"])


def test_moderate_welding_deformation_is_absorbed_by_the_spacing_design():
    """Half the example's welding deformation (end lift about 21 mm): the as-built spacing
    overloads blocks, the designed spacing brings every block within strength."""
    result = scaled_welding(0.5)
    sp = result.spacing
    assert 18e-3 < np.max(result.weld_gap) < 23e-3
    assert np.max(sp.initial_utilisation) > 1.4
    assert sp.fits and result.strength["fem+spacing"].passed is True


def test_severe_welding_deformation_is_diagnosed_as_gap_governed(welded):
    """The full welding deformation (end lift 41 mm) cannot be compensated by spacing alone:
    blocks are kept out of contact and the message says that the gaps govern."""
    sp = welded.spacing
    assert not sp.fits
    assert np.count_nonzero(sp.final_evaluation.reactions <= 0) > 0
    assert "gaps rather than the spacing govern" in sp.message
```

### B.5.9 tests/test_example_case.py

129 lines

```python
from pathlib import Path

import numpy as np
import pytest
import yaml

from shipsupport import case_from_dict, load_case, run_case
from shipsupport.cli import main
from shipsupport.report import format_summary
from shipsupport.units import TONNE_FORCE

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "example_case.yaml"


@pytest.fixture(scope="module")
def result():
    return run_case(load_case(EXAMPLE))


def test_example_equilibrium(result):
    assert result.load.total == pytest.approx(4500 * TONNE_FORCE, rel=1e-12)
    for r in result.supports.values():
        assert abs(r.force_residual) < 1e-9
        assert abs(r.moment_residual) < 1e-9
        assert np.all(r.reactions >= 0)
    for s in result.strength.values():
        assert abs(s.closure_shear) < 1e-6 and abs(s.closure_moment) < 1e-6


def test_summary_and_cli(result, tmp_path):
    text = format_summary(result)
    assert "SUPPORT SOLUTIONS" in text and "LANDING" not in text and "SHIM" not in text
    assert main(["run", str(EXAMPLE), "-o", str(tmp_path)]) == 0
    for name in ("summary.txt", "blocks.csv", "curves_fem.csv", "curves_rigid.csv",
                 "results.png", "block_forces.png"):
        assert (tmp_path / name).exists()
    for name in ("landing.csv", "landing.png", "layout.csv", "welding_and_shims.png",
                 "spacing.csv"):
        assert not (tmp_path / name).exists()
    header = (tmp_path / "blocks.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "shim" not in header and "strength_t" in header


def test_cli_survey_case_outputs(tmp_path):
    assert main(["run", str(EXAMPLE.with_name("welding_survey_case.yaml")), "-o", str(tmp_path),
                 "--methods", "fem"]) == 0
    assert (tmp_path / "welding_data.png").exists()
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "identified angles" in summary


def test_cli_welding_case_outputs(tmp_path):
    assert main(["run", str(EXAMPLE.with_name("welding_case.yaml")), "-o", str(tmp_path)]) == 0
    assert (tmp_path / "welding_data.png").exists()
    assert not (tmp_path / "spacing.csv").exists()               # spacing_design.enabled: false
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "WELDING DEFORMATION" in summary and "SHIM" not in summary


def test_weight_curve_case_runs(tmp_path):
    case = load_case(EXAMPLE.with_name("weight_curve_case.yaml"))
    ship = case.ship
    assert ship.total_weight == pytest.approx(4643.0 * TONNE_FORCE, rel=1e-12)
    assert ship.lcg == pytest.approx(56.639, abs=5e-4)
    curve = ship.weight_curve
    assert curve.x.size == 18
    assert sum(w.weight for w in curve.items()) == pytest.approx(4343.0 * TONNE_FORCE, rel=1e-12)
    assert main(["run", str(EXAMPLE.with_name("weight_curve_case.yaml")), "-o", str(tmp_path),
                 "--methods", "fem", "--no-plots"]) == 0
    summary = (tmp_path / "summary.txt").read_text(encoding="utf-8")
    assert "of which curve" in summary and "BLOCK SPACING DESIGN" in summary


def items_as_curve(items):
    """The piecewise-linear curve (with jumps) that is the sum of the trapezoidal items."""
    from shipsupport.weights import trapezoid_ends
    parts = [(w.xa, w.xb, *trapezoid_ends(w)) for w in items]
    breaks = sorted({p[0] for p in parts} | {p[1] for p in parts})

    def value(x, side):
        q = 0.0
        for xa, xb, wa, wb in parts:
            inside = xa < x <= xb if side < 0 else xa <= x < xb
            if inside:
                q += wa + (wb - wa) * (x - xa) / (xb - xa)
        return q

    points = []
    for x in breaks:
        left, right = value(x, -1), value(x, +1)
        points.append({"x": x, "q_t_per_m": left / TONNE_FORCE})
        if right != left:
            points.append({"x": x, "q_t_per_m": right / TONNE_FORCE})
    return points


def test_weight_curve_equivalent_to_the_example_items_gives_identical_reactions(result):
    """The example's distributed items as one weight curve (sum of trapezoids = a piecewise-
    linear curve with jumps) give the same load and the same block reactions."""
    data = yaml.safe_load(EXAMPLE.read_text(encoding="utf-8"))
    case = load_case(EXAMPLE)
    distributed = [w for w in case.ship.weights if not w.is_point]
    data["weight_curve"] = {"name": "items as curve", "points": items_as_curve(distributed)}
    data["weights"] = [w for w in data["weights"] if "x" in w]           # the point load stays
    curved = run_case(case_from_dict(data))
    assert curved.load.total == pytest.approx(result.load.total, rel=1e-12)
    assert curved.load.lcg == pytest.approx(result.load.lcg, rel=1e-12)
    for m in result.supports:
        assert np.allclose(curved.supports[m].reactions, result.supports[m].reactions,
                           rtol=1e-9, atol=1e-9 * result.load.total)


def test_blocks_built_in_code_get_labels():
    # regression: blocks without a label were shown with empty names
    import dataclasses
    from shipsupport.model import Block
    case = load_case(EXAMPLE)
    case.blocks = [Block(b.x, b.type, b.gap) for b in reversed(case.blocks)]
    case.blocks[0] = dataclasses.replace(case.blocks[0], gap=0.05)      # forward end block lifts off
    result = run_case(case)
    labels = result.strength["fem"].blocks.label
    assert labels[0] == "B01" and labels[-1] == f"B{len(labels):02d}"
    assert any("carry no load" in w and f"B{len(labels):02d}" in w for w in result.warnings)


def test_negative_weight_tolerance_rejected():
    from shipsupport.model import InputError, Settings
    with pytest.raises(InputError):
        Settings(weight_tolerance=-0.1)
```

### B.5.10 tests/test_gui.py

418 lines

```python
import os
import time
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("PyQt5")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import QEventLoop, QTimer  # noqa: E402
from PyQt5.QtWidgets import QApplication, QGroupBox  # noqa: E402

from shipsupport import case_from_dict, load_case, run_case  # noqa: E402
from shipsupport.units import TONNE_FORCE  # noqa: E402
from shipsupport_gui.input_panels import CaseEditor, format_layers, parse_layers  # noqa: E402
from shipsupport_gui.main_window import EXAMPLES, MainWindow  # noqa: E402
from shipsupport_gui.plan_panel import COL_GAP, COL_SPACING, COL_STRENGTH, COL_U, COL_X  # noqa: E402
from shipsupport_gui.worker import AnalysisThread, run_analysis  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "example_case.yaml"
SPACING = ROOT / "examples" / "spacing_case.yaml"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def reference():
    return run_case(load_case(EXAMPLE))


def example_dict(path=EXAMPLE):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def assert_same_result(a, b):
    assert a.supports.keys() == b.supports.keys()
    for m in a.supports:
        assert np.allclose(a.supports[m].reactions, b.supports[m].reactions, rtol=1e-12)


def wait_for(window, app, seconds=300):
    end = time.time() + seconds
    while window.thread is not None and window.thread.isRunning() and time.time() < end:
        app.processEvents()
        time.sleep(0.02)
    for _ in range(30):
        app.processEvents()


def tab_titles(tabs):
    return [tabs.tabText(i).replace("&&", "&") for i in range(tabs.count())]


def test_layers_text_round_trip():
    layers = parse_layers("0.9@30000; 0.2@1000;0.1@400", "t")
    assert layers[2] == {"thickness": 0.1, "E_MPa": 400.0}
    assert parse_layers(format_layers(layers), "t") == layers
    with pytest.raises(ValueError):
        parse_layers("0.9-30000", "t")


def test_tabs_as_built_and_no_removed_widgets(app):
    window = MainWindow()
    assert tab_titles(window.editor) == ["Ship", "Weights", "Hull sections", "Blocks",
                                         "Structure & limits", "Welding deformation", "Settings"]
    assert tab_titles(window.results) == ["Summary", "Block forces", "Plots", "Block reactions",
                                          "Checks", "Welding data", "Spacing design", "Warnings"]
    for name in ("landing_box", "landing_form", "hydro", "optimization_box", "optimization_form",
                 "layout_box", "layout_form", "strengths", "block_groups", "gaps",
                 "layout_preview", "strength_tab"):
        assert not hasattr(window.editor, name), name
    for name in ("landing_plot", "layout_plot", "weld_plot"):
        assert not hasattr(window.results, name), name
    plan = window.editor.plan_editor
    assert not hasattr(plan, "form") and not hasattr(plan, "copy_requested")
    assert not hasattr(window, "sample_layout_action")
    titles = [b.title() for b in window.findChildren(QGroupBox)]
    assert not any("ptimisation" in t or "landing" in t.lower() or "layout" in t.lower()
                   for t in titles)
    menu_texts = [a.text() for a in window.findChildren(type(window.run_action))]
    assert not any("strength &&" in t or "layout" in t.lower() for t in menu_texts)
    assert [path.name for _, path in EXAMPLES] == ["example_case.yaml", "spacing_case.yaml",
                                                   "weight_curve_case.yaml", "welding_case.yaml",
                                                   "welding_survey_case.yaml"]
    window.modified = False
    window.close()


@pytest.mark.parametrize("name", ["example_case", "welding_case", "welding_survey_case",
                                  "spacing_case", "weight_curve_case"])
def test_every_example_round_trips_through_the_editor(app, name):
    from shipsupport.loader import inline_data_files
    path = ROOT / "examples" / f"{name}.yaml"
    data = inline_data_files(example_dict(path), path.parent)
    editor = CaseEditor()
    editor.load_dict(data)
    via_gui = case_from_dict(editor.to_dict())
    direct = load_case(path)
    assert [b.x for b in via_gui.blocks] == [b.x for b in direct.blocks]
    assert [b.effective_strength for b in via_gui.blocks] == \
        [b.effective_strength for b in direct.blocks]
    assert [(w.name, w.weight, w.lcg, w.xa, w.xb) for w in via_gui.ship.weights] == \
        [(w.name, w.weight, w.lcg, w.xa, w.xb) for w in direct.ship.weights]
    assert (via_gui.spacing is None) == (direct.spacing is None)


def test_editor_round_trip_gives_identical_results(app, reference):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    assert_same_result(run_analysis(editor.to_dict()), reference)


def test_weights_tab_items_with_moment_and_weight_curve(app):
    editor = CaseEditor()
    editor.load_dict(example_dict(ROOT / "examples" / "weight_curve_case.yaml"))
    assert editor.curve.table.rowCount() == 18 and editor.curve_name.text() == "lightship curve"
    assert "Σ weight = 4643.0 t, LCG = 56.639 m" in editor.weight_summary.text()
    assert "curve 4343.0 t" in editor.weight_summary.text()
    data = editor.to_dict()
    assert data["weight_curve"]["name"] == "lightship curve"
    assert data["weights"][0]["moment_tm"] == 14000 and "lcg" not in data["weights"][0]
    assert data["weights"][2] == {"name": "stores", "weight_t": 40.0, "moment_tm": 1200.0}
    # an item given by its LCG on the example: the moment column gives the same result
    editor = CaseEditor()
    editor.load_dict(example_dict())
    table = editor.weights.table
    headers = [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]
    assert headers == ["Name", "Weight [t]", "LCG [m]", "Moment [t·m]", "xa [m]", "xb [m]"]
    table.item(1, 2).setText("")
    table.item(1, 3).setText(str(380 * 22.0))
    assert "Σ weight = 4500.0 t, LCG = 50.469 m" in editor.weight_summary.text()
    case = case_from_dict(editor.to_dict())
    assert case.ship.weights[1].lcg == pytest.approx(22.0)
    table.item(1, 2).setText("23.0")                               # disagrees with the moment
    with pytest.raises(Exception, match="does not agree"):
        case_from_dict(editor.to_dict())


def test_blocks_tab_types_plan_and_count_spacing_strength_gap(app):
    editor = CaseEditor()
    editor.load_dict(example_dict(SPACING))
    assert editor.indexOf(editor.blocks_tab) == 3
    plan = editor.plan_editor
    assert plan.count.value() == 44 and plan.block_type.currentText() == "keel"
    assert [plan.block_type.itemText(i) for i in range(plan.block_type.count())] == ["keel"]
    plan.count.setValue(46)                                  # rows added, copying the last row
    d = plan.plan()
    assert d["count"] == 46 and d["spacing"] == 2.4 and d["strength_t"][-1] == 190.0
    plan.count.setValue(10)
    plan.all_spacing.setText("3.0")
    plan._set_all(COL_SPACING, plan.all_spacing)
    plan.all_strength.setText("250")
    plan._set_all(COL_STRENGTH, plan.all_strength)
    plan.table.item(4, COL_SPACING).setText("1.6")           # edit one spacing
    plan.table.item(2, COL_GAP).setText("0.5")               # and one initial gap
    plan.table.item(7, COL_STRENGTH).setText("")             # empty = type capacity
    case = case_from_dict(editor.to_dict())
    xs = [b.x for b in case.blocks]
    assert len(xs) == 10 and np.allclose(np.diff(xs), [3, 3, 3, 3, 1.6, 3, 3, 3, 3])
    assert [round(b.effective_strength / TONNE_FORCE) for b in case.blocks] == \
        [250] * 7 + [200] + [250] * 2
    assert case.blocks[2].gap == pytest.approx(0.5e-3)
    assert float(plan.table.item(9, COL_X).text()) == pytest.approx(xs[-1])   # x column
    assert not plan.table.item(0, COL_X).flags() & 2                          # x read-only
    plan.table.item(3, COL_SPACING).setText("abc")
    with pytest.raises(ValueError, match="B04 spacing"):
        editor.to_dict()


def test_hidden_spacing_design_round_trips(app, tmp_path):
    window = MainWindow()
    assert window.open_path(SPACING)
    loaded = example_dict(SPACING)["spacing_design"]
    assert window.editor.to_dict()["spacing_design"] == loaded
    target = tmp_path / "saved.yaml"
    assert window._write(target)
    assert example_dict(target)["spacing_design"] == loaded
    assert window.editor.auto_spacing_design() == {**loaded, "enabled": True}
    window.open_path(EXAMPLE)                               # no spacing_design in the file
    assert "spacing_design" not in window.editor.to_dict()
    assert window.editor.auto_spacing_design() == {"enabled": True, "min_spacing": 0.8}
    window.editor.structure_box.setChecked(False)
    assert window.editor.auto_spacing_design() == {"enabled": True}
    window.modified = False
    window.close()


def test_f5_does_not_run_the_spacing_design(app):
    window = MainWindow()
    assert window.open_path(SPACING)
    window.run_analysis()
    wait_for(window, app)
    assert window.result is not None and window.result.spacing is None
    plan = window.editor.plan_editor
    assert plan.table.item(0, COL_U).text() != ""           # reactions shown in the plan
    window.modified = False
    window.close()


def test_auto_calculate_spacing_writes_the_design_into_the_plan(app):
    window = MainWindow()
    window.open_path(SPACING)
    window.editor.plan_editor.auto_button.click()
    wait_for(window, app)
    sp = window.result.spacing
    assert sp is not None and sp.fits
    case = case_from_dict(window.editor.to_dict())           # the plan now holds the design
    assert np.allclose([b.x for b in case.blocks], [b.x for b in sp.final_blocks], atol=1e-6)
    assert [b.effective_strength for b in case.blocks] == \
        pytest.approx([b.effective_strength for b in sp.final_blocks])
    plan = window.editor.plan_editor
    assert "Fits" in plan.status.text()
    u = [float(plan.table.item(r, COL_U).text()) for r in range(44)]
    assert max(u) <= 1.0 and max(u) == pytest.approx(np.max(sp.final_utilisation), abs=1e-3)
    assert not window.results.spacing_plot.action.isVisibleTo(window.results.spacing_plot)
    assert window.editor.to_dict()["spacing_design"] == example_dict(SPACING)["spacing_design"]
    window.modified = False
    window.close()


def test_auto_calculate_spacing_without_settings_uses_the_floor_spacing(app):
    window = MainWindow()
    window.open_path(EXAMPLE)
    plan = window.editor.plan_editor
    plan.count.setValue(30)
    plan.all_spacing.setText("3.2")
    plan._set_all(COL_SPACING, plan.all_spacing)
    window.auto_spacing()
    wait_for(window, app)
    sp = window.result.spacing
    assert sp is not None and window.result.case.spacing.min_spacing == pytest.approx(0.8)
    assert np.all(sp.spacing >= 0.8 - 1e-6)
    window.modified = False
    window.close()


def test_survey_case_opens_with_points_in_table(app):
    window = MainWindow()
    survey_case = EXAMPLE.with_name("welding_survey_case.yaml")
    assert window.open_path(survey_case)
    assert window.editor.weld_survey.table.rowCount() == 25
    via_gui = run_analysis(window.editor.to_dict())
    direct = run_case(load_case(survey_case))
    assert np.allclose(via_gui.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-12)
    window.show_result(via_gui)
    assert not window.results.weld_data_plot.message.isVisible()
    window.modified = False
    window.close()


def test_sample_welding_button_fills_and_shows_welding_tab(app):
    window = MainWindow()
    assert window.open_path(EXAMPLE)                  # example without welding data
    window.show_result(run_analysis(window.editor.to_dict()))
    assert window.results.weld_data_plot.action.isVisibleTo(window.results.weld_data_plot)
    window.editor.load_sample_welding(0.0, 120.0)     # what the button does before running
    result = run_analysis(window.editor.to_dict())
    window.show_result(result)
    assert result.deformation is not None
    assert not window.results.weld_data_plot.message.isVisible()
    assert not window.results.weld_data_plot.action.isVisibleTo(window.results.weld_data_plot)
    direct = run_case(load_case(EXAMPLE.with_name("welding_case.yaml")))
    assert np.allclose(result.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-12, atol=1e-6)
    window.modified = False
    window.close()


def test_sample_welding_is_mapped_onto_the_hull(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.load_sample_welding(-5.0, 115.0)
    joints = [row["x"] for row in editor.weld_joints.rows()]
    zones = editor.weld_zones.rows()
    assert joints == pytest.approx([15.0, 35.0, 55.0, 75.0, 95.0])
    assert zones[0]["xa"] == pytest.approx(-5.0) and zones[-1]["xb"] == pytest.approx(115.0)


def test_unresolved_survey_file_is_reported(app):
    data = example_dict()
    data["welding_deformation"] = {"measured": {"file": "survey.csv"}}
    with pytest.raises(ValueError):
        CaseEditor().load_dict(data)


def test_bad_cell_reports_location(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.weights.table.item(1, 1).setText("abc")
    with pytest.raises(ValueError, match="row 2"):
        editor.to_dict()


def test_point_load_needs_both_or_neither_extent(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    editor.weights.table.item(0, 5).setText("")
    with pytest.raises(ValueError, match="both xa and xb"):
        editor.to_dict()


def test_main_window_show_result_and_save(app, reference, tmp_path):
    window = MainWindow()
    assert window.open_path(EXAMPLE)
    window.show_result(run_analysis(window.editor.to_dict()))
    assert window.results.blocks.rowCount() == 44
    assert window.results.checks.rowCount() == 14
    assert window.export_action.isEnabled()
    assert "All checks passed" in window.state_label.text() or \
        "CHECKS FAILED" in window.state_label.text()
    assert "shims" not in window.state_label.text()

    target = tmp_path / "saved.yaml"
    assert window._write(target)
    assert not window.modified
    assert_same_result(run_case(load_case(target)), reference)
    window.close()


def test_open_every_example_and_run(app):
    window = MainWindow()
    for _, path in EXAMPLES:
        window.modified = False
        window.open_example(path)
        window.run_analysis()
        wait_for(window, app)
        assert window.result is not None, path
        assert window.result.spacing is None
        assert window.results.blocks.rowCount() == len(window.result.case.blocks)
    window.modified = False
    window.close()


def test_analysis_thread(app, reference):
    got = {}
    loop = QEventLoop()
    thread = AnalysisThread(example_dict())
    thread.succeeded.connect(lambda r: got.setdefault("result", r))
    thread.failed.connect(lambda msg: got.setdefault("error", msg))
    thread.finished.connect(loop.quit)
    QTimer.singleShot(30000, loop.quit)
    thread.start()
    loop.exec_()
    thread.wait()
    assert "error" not in got
    assert_same_result(got["result"], reference)


def test_analysis_thread_reports_input_error(app):
    data = example_dict()
    data["weights"][0]["lcg"] = 1.0   # outside the middle third
    got = {}
    loop = QEventLoop()
    thread = AnalysisThread(data)
    thread.failed.connect(lambda msg: got.setdefault("error", msg))
    thread.finished.connect(loop.quit)
    QTimer.singleShot(30000, loop.quit)
    thread.start()
    loop.exec_()
    thread.wait()
    assert "middle third" in got["error"]


# --------------------------------------------------------------------------- weight curve sample

def test_sample_weight_curve_matches_the_example_file(app):
    from shipsupport_gui.input_panels import SAMPLE_WEIGHT_CURVE
    data = yaml.safe_load((ROOT / "examples" / "weight_curve_case.yaml").read_text(encoding="utf-8"))
    points = [(p["x"], p["q_t_per_m"]) for p in data["weight_curve"]["points"]]
    assert SAMPLE_WEIGHT_CURVE["points"] == points
    assert [(i["name"], i["weight_t"], i["moment_tm"]) for i in SAMPLE_WEIGHT_CURVE["items"]] == \
        [(i["name"], i["weight_t"], i["moment_tm"]) for i in data["weights"]]


def test_sample_weight_curve_loads_and_runs(app):
    editor = CaseEditor()
    editor.load_dict(example_dict())
    total, lcg = editor.load_sample_weight_curve()
    assert total == pytest.approx(4643.0) and lcg == pytest.approx(56.639, abs=1e-3)
    result = run_analysis(editor.to_dict())
    assert result.load.total == pytest.approx(4643.0 * TONNE_FORCE, rel=1e-9)
    assert result.case.ship.weight_curve is not None and not result.warnings
    direct = run_case(load_case(ROOT / "examples" / "weight_curve_case.yaml"))
    assert np.allclose(result.supports["fem"].reactions, direct.supports["fem"].reactions,
                       rtol=1e-9)


def test_sample_weight_curve_is_mapped_onto_another_hull(app):
    editor = CaseEditor()
    data = example_dict()
    data["ship"].update(lpp=150.0, x_fwd=150.0)
    editor.load_dict(data)
    total, lcg = editor.load_sample_weight_curve()
    assert total == pytest.approx(4643.0, rel=1e-6)                 # weights kept
    assert lcg == pytest.approx(56.639 * 150.0 / 120.0, abs=2e-3)    # positions scaled
    assert editor.curve.rows()[-1]["x"] == pytest.approx(150.0)


def test_block_forces_and_spacing_charts_show_the_weight_distribution():
    from matplotlib.figure import Figure

    from shipsupport.report import draw_block_forces, draw_spacing
    result = run_case(load_case(ROOT / "examples" / "weight_curve_case.yaml"), ["fem"])
    for draw in (draw_block_forces, draw_spacing):
        fig = Figure()
        draw(fig, result)
        top = fig.axes[0]
        assert top.get_ylabel() == "load [t/m]"
        labels = [line.get_label() for line in top.get_lines()]
        assert "total distributed load w(x)" in labels
        assert any(label.startswith("weight curve") for label in labels)
```
