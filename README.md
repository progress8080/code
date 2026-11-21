# shipsupport

Dry-docking of a ship on keel blocks: the **reaction force of every block** and the
**automatic design of the block spacing** so that no block carries more than its own
strength — also on a hull that is not straight because of welding deformation.

1. **Input** – ship data, a detailed weight distribution (weight items with LCG *or*
   longitudinal moment, and/or a weight-distribution curve), hull girder sections, block
   types, the **block plan** (number of blocks, spacing, strength and initial gap of every
   block), structure (floors), permissible hull limits, welding deformation (joint kinks,
   shrinkage zones, keel survey with optional joint-angle identification) and settings.
2. **Block reactions** – Method A (rigid hull) and Method B (elastic Timoshenko beam FEM)
   on compression-only, optionally crushing block springs with initial gaps; exact
   active-set iteration with gap continuation for large gaps.
3. **Spacing design** – the user gives the strength of every block; the program computes
   the spacing [m] between the blocks so that no reaction exceeds its block strength
   (count and order of the blocks kept). Every candidate layout is analysed with the full
   support model, including the welding gaps at the candidate block positions.
4. **Checks** – hull girder shear force, bending moment, deck/keel and shear stress against
   the permissible envelopes; block load against strength, bearing pressure, block
   position against the structure.

Theory, methods, implementation and verification are documented in the bachelor's thesis
in [thesis/](thesis/) (built PDF: `thesis/output/thesis.pdf`).

## Install

```bash
pip install -e ".[dev]"      # numpy, scipy, pyyaml + matplotlib, PyQt5, pytest
```

Extras: `plots` (matplotlib only), `gui` (matplotlib, PyQt5), `dev` (GUI + pytest),
`thesis` (GUI + pypdfium2, for building the thesis).

## Command line

```bash
shipsupport run examples/example_case.yaml -o results
# or, without installing (bash):
PYTHONPATH=src python -m shipsupport run examples/example_case.yaml -o results
# or, without installing (PowerShell):
$env:PYTHONPATH = "src"; python -m shipsupport run examples/example_case.yaml -o results
```

Options: `--methods rigid fem`, `--no-plots`, `--strict` (exit code 2 if any check fails
or the spacing design does not fit). Without `-o` the results go to
`<case name>_results/` next to the case file. The spacing design runs in the CLI when the
case contains `spacing_design` with `enabled: true` (the default when the section is given).

## Input (YAML)

Units: lengths m (x from the aft perpendicular, positive forward), weights t, moments t·m,
distributed weights t/m, stresses MPa, moduli as suffixed (`_MPa`, `_GPa`), gaps mm.
All example values are illustrative placeholders, not design values.

| Section | Keys |
|---|---|
| `ship` | `name`, `lpp`, `x_aft`, `x_fwd`, `weight_t`, `lcg` (declared values, cross-checked) |
| `weights` | list of `{name, weight_t, lcg \| moment_tm, xa, xb}`; `x` (or no extent) = point load |
| `weight_curve` | `{name, points: [{x, q_t_per_m}, …]}` or a plain list of points |
| `sections` | list of `{xa, xb, I_m4, E_GPa, shear_area_m2, G_GPa, z_deck_m3, z_keel_m3}` (or `EI_Nm2`, `kGA_N`) |
| `blocks.types` | per name: `contact_width`, `contact_length`, `layers: [{thickness, E_MPa}]` or `stiffness_kN_per_mm`, `capacity_t`, `allowable_pressure_MPa`, `crush_stress_MPa`, `post_crush_ratio` |
| `blocks.plan` | `{type, start_x, count, spacing, strength_t, gap_mm}` (see below) |
| `structure` | `{start, end, spacing}` (start plus two of end/spacing/count) or `{x: [...]}` |
| `permissible` | `bm: [{x, sag_tm, hog_tm}]`, `sf: [{x, sf_t}]` |
| `welding_deformation` | `joints: [{x, angle_mrad}]`, `zones: [{xa, xb, curvature_mm_per_m2 \| deck_shrinkage_mm_per_m, bottom_shrinkage_mm_per_m, depth}]`, `measured: {x, lift_mm \| file, interpolation, identify_joint_angles, regularization}` |
| `spacing_design` | `enabled`, `method`, `target_utilisation`, `min_spacing`, `max_spacing`, `x_min`, `x_max`, `fix_end_blocks`, `snap_to_structure`, `grid`, `max_iterations` |
| `settings` | `methods`, `max_element_length`, `shear_deformation`, `allowable_bending_MPa`, `allowable_shear_MPa`, `structure_tolerance`, `weight_tolerance`, `max_iterations` |

### Weight items and weight curve

```yaml
weights:
  - {name: main engine,   weight_t: 380, lcg: 22.0,         xa: 16.0, xb: 28.0}
  - {name: residual fuel, weight_t: 200, moment_tm: 14000,  xa: 60.0, xb: 80.0}   # LCG 70 m
  - {name: stores,        weight_t: 40,  moment_tm: 1200}                         # point at 30 m
  - {name: rudder,        weight_t: 25,  x: 2.0}                                  # point load
weight_curve:
  name: lightship curve
  points:                  # piecewise linear, zero outside; the same x twice is a jump
    - {x: 0.0,  q_t_per_m: 8.0}
    - {x: 16.0, q_t_per_m: 36.0}
    - {x: 16.0, q_t_per_m: 64.0}
    - {x: 120.0, q_t_per_m: 0.0}
```

An item is spread as a trapezoid over `[xa, xb]` that keeps its weight and LCG (the LCG
must lie in the middle third of the extent). The LCG is given as `lcg` or as the moment
about AP `moment_tm` (LCG = moment / weight); if both are given they must agree. Every
segment of the weight curve is converted into one exact trapezoidal item
(weight (q_k + q_k+1) L / 2, LCG x_k + L (q_k + 2 q_k+1) / (3 (q_k + q_k+1))), so a curve
and items can be combined freely; x must not decrease and q ≥ 0 (errors name the point).

### Block plan

```yaml
blocks:
  plan:
    type: keel
    start_x: 6.4           # m, first block
    count: 44
    spacing: 2.4           # m to the next block: one value or a list of count-1 values
    strength_t: [190, 190, 190, 160, …]   # optional: one value or count values (default: type capacity)
    gap_mm: 0.0            # optional initial gap of every block (one value or count values)
```

The blocks are B01 … BNN from aft to forward. A block's gap moves with the block when the
spacing is designed.

### Spacing design

```yaml
spacing_design:
  enabled: true
  min_spacing: 0.8         # m (e.g. the floor spacing)
  max_spacing: 4.8         # m (default: twice the given mean spacing)
  x_min: 6.4               # m, allowed extent of the block line (default: first/last block)
  x_max: 109.6
  snap_to_structure: true  # final positions under floors
```

The number of blocks and their order (so which block has which strength) are kept; only
the spacing is designed, minimising the largest reaction/strength:

1. strength-weighted equidistribution (fixed point) as a start,
2. sequential linear programming (HiGHS) on the minimax problem with a trust region,
3. rounding to the floors (or a grid) and a 1-opt local search.

Welding deformation is a constraint: every candidate layout is analysed with the welding
gaps recomputed at its block positions. If the plan cannot fit, the summary says why
(total strength below the load, end block governing → longer block line) and suggests a
number of blocks.

## Examples

| File | Content |
|---|---|
| `examples/example_case.yaml` | 120 m, 4500 t example vessel on a plan of 44 keel blocks (6.4 m, 2.4 m spacing, type capacity 200 t) |
| `examples/spacing_case.yaml` | main design example: 44 blocks with their own strengths (190 t, 160 t for B04–B12 under the engine room) and spacing design: max reaction/strength 1.34 → 0.76, fits |
| `examples/weight_curve_case.yaml` | lightship as an 18-point weight curve with jumps (4343 t) plus deadweight items given by moment (300 t); 4643 t in all; spacing design fits (1.06 → 0.57) |
| `examples/welding_case.yaml` | the example vessel built from six welded construction blocks (joint kinks, shrinkage); `spacing_design` included but switched off |
| `examples/welding_survey_case.yaml` + `examples/welding_survey.csv` | joint angles identified from a (synthetic) keel survey |

## Outputs (results folder)

| File | Content |
|---|---|
| `summary.txt` | weights (curve and items), support solutions, reactions per block, hull and block checks, welding deformation, spacing design, warnings |
| `blocks.csv` | per block: gap, welding gap, stiffness, strength, reaction (each method), pressure, utilisation, compression, state, structure offset |
| `curves_<method>.csv` | V(x), M(x), deck/keel stress, shear stress, utilisation (`curves_fem_spacing.csv` for the designed spacing) |
| `spacing.csv`, `spacing_history.csv` | designed and given position, spacing, strength, reaction and reaction/strength of every block; convergence (only with a spacing design) |
| `results.png`, `block_forces.png` | load, reactions, settlement, V and M; weight distribution (top) and force on every block |
| `welding_data.png` | welding deformation input, keel lift and block gaps (only with welding input) |
| `spacing_design.png` | weight distribution with given and designed blocks, given vs designed spacing, reactions vs strength, convergence (only with a spacing design) |

## Desktop GUI (PyQt5)

```bash
pip install -e ".[gui]"
shipsupport-gui                       # or: shipsupport-gui examples/spacing_case.yaml
PYTHONPATH=src python -m shipsupport_gui   # without installing (bash)
```

- **Input tabs:** Ship · Weights (items with LCG *or* moment, the weight-curve points table
  with a *Load sample curve* button, live Σ weight / LCG of both) · Hull sections · Blocks (block types with layers written
  `t@E_MPa; …`, and the block plan: number of blocks, first block x, type, *Set all*
  spacing/strength, a per-block table with editable spacing, strength and gap and read-only
  x, reaction and reaction/strength, the **Auto-calculate spacing** button and its status) ·
  Structure & limits · Welding deformation (joints, zones, keel survey with *Import CSV…*,
  identification options) · Settings. Tables accept Ctrl+C / Ctrl+V / Delete like Excel.
- **Result tabs:** Summary · Block forces (weight distribution on top) · Plots (zoom/pan/save toolbar) · Block reactions
  (utilisation colour-coded) · Checks · Welding data (offers *Add sample welding data and
  run* when the case has none) · Spacing design (weight distribution with the given and designed
  blocks on top) · Warnings.
- **F5** runs the analysis in a background thread (analysis only); **Auto-calculate
  spacing** (button, toolbar or Ctrl+F5) runs the spacing design and writes the designed
  spacing into the block plan. The GUI has no spacing-design settings form: a
  `spacing_design` section of the opened case is used by the button and saved back
  unchanged; without one the defaults are used with the floor spacing as minimum spacing.
- *File → Open example* opens the five examples; *File → Add sample weight distribution
  curve* replaces the weights by a sample lightship curve (18 points, 5 jumps) and three
  deadweight items given by moments (the data of `weight_curve_case.yaml`, mapped onto the
  current hull, declared weight/LCG updated) and runs; *File → Add sample welding data* fills the
  Welding deformation tab for the current hull; *File → Export results* writes the summary,
  CSV tables and plots. Cases are saved as the same YAML the CLI reads (comments are not
  preserved; a survey CSV referenced by the case is stored inline).

The GUI (`src/shipsupport_gui/`) only edits the case and displays results; all
calculations go through the `shipsupport` package.

## Use from Python

```python
from shipsupport import load_case, run_case
from shipsupport.report import format_summary

result = run_case(load_case("examples/spacing_case.yaml"))
print(format_summary(result))
result.supports["fem"].reactions       # N, one per block of the plan (sorted by x)
result.spacing.final_blocks            # designed block positions (with their strengths)
result.spacing.spacing                 # m, designed distance from each block to the next
result.strength["fem"].M               # N·m, sagging positive
```

## Conventions and scope

SI internally. x from AP, positive forward. Displacement positive downward. Sagging bending
moment positive. Stresses positive in tension.

- Longitudinal (2-D) model of the keel-block line; transverse load sharing, side blocks,
  heel and lateral loads are not modelled.
- Dry-dock condition (no buoyancy); floating-dock flexibility can be included in the block
  stiffness.
- The example input contains placeholder numbers — not design values.

## Tests

```bash
pytest            # the GUI tests run offscreen and are skipped if PyQt5 is missing
```

Covers exact trapezoid and weight-curve integration (segments, jumps, invalid points),
moment ⇄ LCG input, the block plan, closed-form beam deflection (Euler–Bernoulli and
Timoshenko), continuous-beam reactions, rigid-limit convergence, lift-off/gap equivalence,
crushing, mesh independence, shear/moment closure, welding-deformation shapes and survey
identification, gap continuation, the spacing design (analytic optimum, weak blocks,
limits, grid, fixed ends, infeasible plans, the spacing example and the welded hull with
recomputed welding gaps), a weight curve equivalent to the example items (identical
reactions), GUI ↔ YAML round trips of every example and end-to-end CLI runs.
