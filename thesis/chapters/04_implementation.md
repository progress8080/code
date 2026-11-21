# 4 Software Implementation

This chapter describes how the models and algorithms of Chapter 3 were turned into a maintainable and verifiable program. The program consists of two Python packages. The calculation package `shipsupport` contains the data model, the input reader, the weight-distribution model, the welding-deformation model, the two support solvers with the shared contact iteration, the hull girder and block checks, the automatic spacing design and the reporting functions. The package `shipsupport_gui` is a desktop front end built on PyQt5.

The chapter treats the program as a software-engineering artefact in its own right (objective O5). It first describes the development process and the requirements, then justifies the technology choices, the architecture and the principal design decisions. The data model, the input format and the implementation of the solvers, the welding-deformation module and the spacing-design module follow. The chapter closes with the reporting functions, the two user interfaces, a requirements traceability matrix and code metrics.

## 4.1 Development process

The program was developed iteratively and incrementally [@sommerville2016]. Each increment added one capability on top of a working, tested core:

1. the longitudinal weight distribution from weight items and the rigid-hull method (Method A);
2. the elastic-hull finite element method (Method B) with compression-only, crushable blocks;
3. the hull girder and block checks;
4. reporting, the command-line interface and the graphical user interface;
5. the welding-deformation model of the assembled hull, including the identification of joint angles from a keel survey;
6. a robust contact solution for large gaps (gap continuation);
7. the block plan, in which the number, spacing, strength and initial gap of every block are given directly;
8. the automatic spacing design from the block strengths, developed test-first against a substitute analysis (Section 4.10);
9. weight items given by their longitudinal moment and the weight-distribution curve.

Within each increment the work followed a test-first discipline where this was practical. Before a numerical module was considered complete, at least one test compared it with an analytical solution or with an independent computation (Chapter 5). The complete test suite was run after every change, so that a later increment could not silently break an earlier one. Several defects were found in this way; for example, the first version of the contact iteration lost all support when large welding gaps were applied, which led to the gap continuation of Section 4.8.5. The static analyser ruff was used to keep the code free of unused names, shadowed variables and similar defects that tests do not reveal.

## 4.2 Requirements

The requirements were derived from the problem statement of Section 1.2 and from the intended use of the program: a transparent design tool for students and engineers, whose final result is the spacing of the keel blocks, and which a reader can check line by line against the theory. They are listed in Table 4.1. Functional requirements (F) describe what the program must calculate or provide; non-functional requirements (N) describe properties of the software itself.

| ID | Requirement |
|---|---|
| F1 | Build the longitudinal weight distribution from weight items given with their LCG or their longitudinal moment about AP (distributed as trapezoids, or as point loads), preserving the total weight and the LCG exactly |
| F2 | Accept a weight-distribution curve (piecewise-linear $q(x)$ with jumps), alone or together with weight items, and integrate it exactly |
| F3 | Define the blocks by a block plan: number of blocks, position of the first block, spacing to the next block, strength of every block (default: capacity of the block type) and initial gap of every block |
| F4 | Compute the block reactions with a rigid-hull model (Method A) |
| F5 | Compute the block reactions with an elastic hull girder model including shear deformation (Method B) |
| F6 | Treat blocks as compression-only springs with initial gaps and optional bilinear crushing, and solve the contact problem reliably also for gaps of several centimetres |
| F7 | Compute shear force, bending moment and hull girder stresses, and check them against permissible values |
| F8 | Check block bearing pressure, block load against the strength of each block, crushing, and the block position relative to the transverse structure |
| F9 | Model the welding deformation of the assembled hull from erection-joint kinks and construction-block shrinkage, and convert it into block gaps |
| F10 | Accept a surveyed keel line (inline or as a CSV file), interpolate it, and optionally identify the joint kink angles from it |
| F11 | Compute the spacing of the blocks of the plan (number and order of the blocks kept) at which no block reaction exceeds the strength of its block, within the spacing limits, the allowed extent and the floor positions, with the welding gaps recomputed for every candidate layout; diagnose the cause when no such spacing exists |
| F12 | Produce a text summary, CSV tables and charts, including the force on every block and the designed spacing |
| F13 | Offer a command-line interface with exit codes suitable for scripts |
| F14 | Offer a graphical user interface with a block plan editor and a one-click automatic spacing design that writes the designed spacing back into the plan |
| N1 | Use SI units internally and convert only at the input/output boundary |
| N2 | Reject invalid or inconsistent input with a message that names the offending item |
| N3 | Keep the calculation core independent of the user interfaces |
| N4 | Make results verifiable: equilibrium residuals and closure errors reported at run time |
| N5 | Provide an automated test suite with analytical benchmarks |
| N6 | Analyse a typical case in well under one second, and design its spacing in a few seconds, on a personal computer |
| N7 | Run on Windows, Linux and macOS with freely available libraries |
| N8 | Correctness of the design: every designed spacing is the result of a full analysis with the support model, not of a surrogate |
| N9 | Traceability: every requirement realised by identifiable modules and covered by tests |
| N10 | Maintainability: small modules with one responsibility each, readable against the theory |
| N11 | Testability of the design algorithm independently of the solvers |

Table: Table 4.1 – Functional (F) and non-functional (N) requirements.

The modules that realise each requirement and the tests that cover it are given in the traceability matrix of Section 4.14. Requirement N6 was checked by timing: a complete analysis of the example vessel with both methods took 0.017 s, and the analysis including the spacing design took 0.76 s for the spacing example and 1.97 s for the weight-curve example (Section 5.9).

## 4.3 Technology choices

Python was chosen as the implementation language. It is widely taught in engineering curricula, its syntax is close to the mathematical notation of Chapter 3, and it has mature, open-source libraries for every task required. The libraries and the reasons for choosing them are summarised in Table 4.2.

| Component | Library | Purpose in the program | Reason for the choice |
|---|---|---|---|
| Arrays and dense linear algebra | NumPy [@harris2020numpy] | load arrays, element matrices, rigid-hull 2×2 solve, least-squares identification, sensitivities | de facto standard; vectorised operations avoid slow Python loops |
| Sparse linear algebra | SciPy [@virtanen2020scipy] | assembly in coordinate format; sparse direct solution of the FEM system | provides `scipy.sparse` and `spsolve`, based on SuperLU [@li2005superlu] |
| Linear programming | SciPy `optimize.linprog` with HiGHS [@huangfu2018] | linearised minimax subproblem of the spacing design | high-performance open-source LP solver bundled with SciPy |
| Interpolation | SciPy `interpolate` | surveyed keel line (PCHIP [@fritsch1980], natural cubic spline) | shape-preserving interpolation without overshoot |
| Input files | PyYAML [@pyyaml] | reading and writing case files | human-readable, supports comments, maps directly to dictionaries |
| Plotting | Matplotlib [@hunter2007matplotlib] | result plots, block-force chart, welding data chart, spacing design chart | publication-quality output; embeddable in Qt |
| Graphical interface | PyQt5 [@pyqt5] | desktop application | mature cross-platform toolkit with table widgets and threading |
| Testing | pytest [@pytest] | automated verification suite | concise tests, fixtures, parametrisation |

Table: Table 4.2 – Libraries used and the reasons for their selection.

Several alternatives were considered. A spreadsheet implementation would have been easier to distribute, but it cannot reasonably handle a sparse finite element system with a contact iteration, let alone a design loop that performs hundreds of such analyses. A compiled language such as C++ or Fortran would run faster, but one analysis of the example vessel takes about 4 ms, so speed was not a deciding factor, and the code would be much less readable for the intended audience. A commercial optimisation package was not needed: the linearised subproblem of the spacing design is a linear programme with about $3N$ variables that HiGHS solves in milliseconds. A web front end was rejected in favour of a desktop application because docking plans and yard survey data are often confidential, and a local program does not transmit them anywhere.

## 4.4 Architecture

The program follows a layered architecture (Figure 4.1). Each layer uses only the layers below it, and the user interfaces sit at the top.

![Figure 4.1: Layered architecture of the program. The numerical core contains the weight model, the welding-deformation model, the solvers, the spacing design and the strength checks. The calculation package has no dependency on the graphical user interface.](figures/fig_architecture.png)

The layers are:

1. **Input and data model** (`model.py`, `loader.py`, `units.py`). Validated data classes describe the ship, its weight items and weight curve, hull sections, block types and blocks, permissible envelopes, welding deformation, the settings of the spacing design and the analysis settings. The loader converts a YAML file, or an equivalent dictionary, into these classes and converts all quantities to SI units.
2. **Numerical core.**
   - `weights.py` builds the finite element mesh and the piecewise-linear load intensity.
   - `deformation.py` evaluates the keel lift of the welded hull, identifies joint angles from a survey and converts the lift into block gaps.
   - `solvers/` contains the compression-only block springs with the contact iteration (`contact.py`), the rigid-hull solver (`rigid.py`), the FEM solver (`fem.py`) and the common result type (`result.py`).
   - `strength.py` computes shear force, bending moment, stresses and all checks.
   - `spacing.py` designs the spacing of the blocks of the plan. It does not import the solvers; the analysis of a candidate layout is passed to it as a function (Section 4.5).
3. **Orchestration and output** (`analysis.py`, `report.py`). The function `run_case` performs the complete analysis of a case in a fixed sequence, runs the spacing design when it is requested and collects warnings. It also creates the analysis function for the spacing design. The report module produces the text summary, the CSV files and the charts.
4. **User interfaces** (`cli.py`, `shipsupport_gui`). The command-line program and the desktop application.

The architecture follows the dependency rule of clean architecture [@martin2017]: source-code dependencies point only inwards, from the interfaces towards the numerical core and the data model. The numerical modules do not know how they are called, and in particular they import nothing from the GUI (requirement N3). Within the core, the spacing design sits beside the solvers rather than above them: it depends only on the data model, and the composition of the design algorithm with a concrete support model takes place in the orchestration layer.

The most important consequence of this structure is that the GUI never performs calculations of its own. It assembles a dictionary with exactly the schema of the YAML input file and passes it to the same function chain as the command-line interface: `case_from_dict` followed by `run_case`. The GUI and the command-line program therefore produce identical results for identical input, which was verified by automated tests (Section 5.8). The core can also be used from scripts, for example for the parameter studies of Chapter 6, without any graphical dependency.

## 4.5 Design decisions

The following design decisions shaped the implementation. Each was taken to support one or more non-functional requirements.

- **Validated data classes, frozen where they describe physical objects (N2, N10).** All input is represented by data classes that check their own consistency at construction. The classes of physical objects — weight items, weight curve, sections, layers, block types, blocks, envelopes, joints, zones and the survey — are declared `frozen`, so an object cannot change after it has been validated; a modified copy must be made explicitly with `dataclasses.replace`. The spacing design relies on this: a candidate layout is a list of copies of the given blocks with new positions, so a block keeps its type, strength and gap, and neither the user's input nor an earlier candidate can be altered by the design.
- **Validation at the boundary (N2).** Every check of the input happens once, when the data classes are constructed by the loader. Code inside the core can rely on consistent input, which keeps the numerical modules short; the core raises errors only for conditions that depend on the solution, such as a support that becomes unstable.
- **SI units only inside the core (N1).** Conversion happens in the loader and in the report. No numerical function has to know about tonnes, millimetres or MPa, which removes a frequent source of errors in engineering software.
- **Pure functions.** The numerical functions take arrays and data objects and return new objects. They hold no hidden state, which makes them easy to test in isolation and to call repeatedly, as the spacing design does hundreds of times.
- **A shared contact strategy (N10).** The contact iteration receives the support model as a function `linear_solve(kt, c, gap)`, which solves the linear problem for given tangent stiffnesses, constant terms and gaps. The rigid-hull and FEM solvers supply different implementations of this function, a form of the Strategy pattern. Gaps, lift-off, crushing and gap continuation are therefore implemented once and behave identically in both methods.
- **Dependency-injected analysis for the spacing design (N8, N11).** The design algorithm `design_spacing` receives a function `evaluate(blocks)` that returns an object with the block reactions of a candidate layout. In the program this function is created by `analysis.spacing_evaluator`, which rebuilds the mesh around the candidate block positions, recomputes the welding gaps at these positions and calls the chosen solver — exactly the analysis that `run_case` performs for the given plan. Every candidate is therefore judged with the full support model (N8). In the tests the function is replaced by a substitute analysis whose optimum is known in closed form, so the algorithm could be developed and tested before it was connected to the solvers, in milliseconds per test (N11).
- **Layered core, CLI and GUI (N3).** The user interfaces are thin: they translate between the user and the case dictionary and display the result objects. All logic that both interfaces need — reading survey files, computing structure positions, drawing charts — lives in the core and is called by both.
- **Verification built into the product (N4).** Every support solution reports its equilibrium residuals and every strength result the closure error of its shear and moment curves. The designed spacing is stored as an ordinary support solution and checked like any other.

## 4.6 Data model and input validation

All input data are represented by Python data classes defined in `model.py`. Table 4.3 lists the classes, their main fields with the units used internally, and the validation performed at construction.

| Class | Main fields (internal units) | Validation performed at construction |
|---|---|---|
| `WeightItem` | `name`, `weight` [N], `lcg`, `xa`, `xb` [m] | non-zero, finite weight; `xb ≥ xa`; LCG in the middle third of the extent; a point load requires `lcg = xa = xb` |
| `WeightCurve` | `name`, `x` [m], `q` [N/m] | at least two points of equal number; finite values; $q \ge 0$; `x` not decreasing (equal `x` = jump); at least one loaded segment |
| `Section` | `xa`, `xb` [m], `EI` [N·m²], `kGA` [N], `z_deck`, `z_keel` [m³], `shear_area` [m²] | `xb > xa`; positive stiffness and section properties |
| `Layer` | `thickness` [m], `E` [Pa] | positive values |
| `BlockType` | contact width and length [m], `layers`, stiffness override [N/m], capacity [N], allowable pressure, crush stress [Pa], post-crush ratio [–] | positive dimensions; layers or override required; post-crush ratio in [0, 1] |
| `Block` | `x` [m], `type`, `gap` [m], `label`, `strength` [N] (optional; default: capacity of the type) | positive strength if given |
| `BendingEnvelope`, `ShearEnvelope` | `x` [m], limits [N·m] or [N] | equal lengths; strictly increasing `x`; positive limits |
| `WeldJoint` | `x` [m], `angle` [rad] | finite values |
| `ShrinkageZone` | `xa`, `xb` [m], `curvature` [1/m] | `xb > xa`; finite curvature |
| `MeasuredProfile` | `x` [m], `lift` [m], interpolation, identification flag, regularisation | equal lengths (at least 2); strictly increasing `x`; known interpolation; regularisation ≥ 0 |
| `WeldingDeformation` | joints, zones, measured profile | at least one source; identification requires joint positions |
| `Ship` | name, $L_{pp}$, hull extent, weights, sections, declared weight and LCG, structure positions, envelopes, weight curve | positive $L_{pp}$; valid extent; at least one weight item and one section; positive total weight |
| `Settings` | methods, maximum element length, shear deformation flag, allowable stresses, tolerances, iteration limit | positive values; non-negative tolerances |
| `SpacingDesign` | support method, target utilisation [–], minimum and maximum spacing [m], allowed extent [m], fixed-end flag, snap-to-structure flag, rounding grid [m], SLP iteration limit | known method; positive target and minimum spacing; `max_spacing ≥ min_spacing`; `x_max > x_min`; positive grid |
| `Case` | ship, blocks, settings, welding deformation, spacing design | at least two blocks |

Table: Table 4.3 – Data classes of the calculation core, their main fields with internal units, and the validation performed.

The weight curve is not a separate load type inside the core. Its method `items` converts every loaded segment into one exact trapezoidal `WeightItem` (Section 4.7.1), and the loader appends these items to the weight items of the ship. The load model, the solvers and the checks are therefore unchanged by the curve; the `Ship` keeps a reference to the curve only for the report.

Two exception classes carry error information to the user. `InputError` (a subclass of Python's `ValueError`) is raised whenever input data are invalid or inconsistent. The message always identifies the item concerned, for example "Weight item 'main engine': LCG … must lie in the middle third of its extent", "weights[1] ('main engine'): moment_tm / weight_t = 23 m does not agree with lcg = 22 m", "Weight curve 'lightship curve', point 3 (x = 5 m): x must not be smaller than at point 2" or "blocks.plan.spacing: 2 values given, 3 expected". `ContactError` is raised when the block state iteration fails to converge, even with gap continuation. Its subclass `UnstableSupportError` is raised when fewer than two separated blocks remain in contact, which would leave the beam with a rigid-body mode (Section 3.7). Both interfaces catch these exceptions: the command-line program prints the message and exits with code 1, and the GUI shows it in a dialog box and in the summary pane.

## 4.7 Input format

A case is stored as a YAML file. Its top-level sections are listed in Table 4.4. Lengths are given in metres. Quantities in other units carry a suffix in the key name, such as `weight_t` (tonnes), `moment_tm` (tonne-metres), `q_t_per_m`, `E_MPa`, `capacity_t`, `strength_t`, `gap_mm` or `angle_mrad`. The suffix tells the user which unit to enter and tells the loader which conversion to apply. Five example cases are supplied (`examples/`): the example vessel, the spacing example with individual block strengths, the weight-curve example, the welded hull and the welded hull with a keel survey; they are analysed in Chapter 6.

| Section | Content | Required |
|---|---|---|
| `ship` | name, $L_{pp}$, hull extent, declared weight and LCG | yes |
| `weights` | weight items (distributed or point), with LCG or moment | yes, unless `weight_curve` is given |
| `weight_curve` | piecewise-linear weight distribution | no |
| `sections` | hull girder properties by longitudinal range | yes |
| `blocks` | block types and the block plan | yes |
| `structure` | positions of floors, bulkheads and web frames | no |
| `permissible` | permissible bending moment and shear force envelopes | no |
| `welding_deformation` | erection-joint kinks, construction-block shrinkage, surveyed keel line | no |
| `spacing_design` | settings of the automatic spacing design | no |
| `settings` | methods, mesh size, allowable stresses, tolerances | no |

Table: Table 4.4 – Top-level sections of the YAML case file.

### 4.7.1 Weights: items, moments and the weight curve

A weight item gives its weight and either its LCG (`lcg`) or its longitudinal moment about AP (`moment_tm`), from which the loader computes the LCG as $M/W$ (Section 3.3.1). Moments are what a weight list in a loading manual often provides, and adding moments is how the total LCG is obtained, so the program accepts them directly. If both are given they must agree to a relative tolerance of $10^{-6}$; otherwise an `InputError` names both values. An item with an extent `xa`, `xb` is distributed as a trapezoid that preserves its weight and LCG (Section 3.3.2); an item with `x`, or without an extent, is a point load. The keys are summarised in Table 4.5.

The weight curve is a list of points $(x_k, q_k)$ with $q$ in tonnes per metre, linear between the points and zero outside the first and last point; the same $x$ twice is a jump of $q$ (Section 3.3.4). Each segment $[x_k, x_{k+1}]$ of length $L_k$ is an exact trapezoid, and `WeightCurve.items` converts it into a weight item with

$$ W_k = \tfrac{1}{2}\,(q_k + q_{k+1})\,L_k, \qquad \bar{x}_k = x_k + \frac{L_k\,(q_k + 2\,q_{k+1})}{3\,(q_k + q_{k+1})} \tag{4.1} $$

Zero-length segments (jumps) and unloaded segments are skipped. Because $q_k, q_{k+1} \ge 0$, the centre $\bar{x}_k$ always lies in the middle third of the segment, so the middle-third condition of the trapezoidal spreading is satisfied automatically and the trapezoid of the item is the segment of the curve itself. Listing 4.1 shows the conversion, and Listing 4.2 the weight input of the weight-curve example, a lightship curve of 18 points with five jumps plus three deadweight items given by their moments.

**Listing 4.1 – Conversion of the weight curve into exact trapezoidal items (`model.py`, class `WeightCurve`).**

```python
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
```

**Listing 4.2 – Weight curve and items with moments in the weight-curve example (`examples/weight_curve_case.yaml`).**

```yaml
weight_curve:
  name: lightship curve
  points:
    - {x: 0.0,   q_t_per_m: 8.0}
    - {x: 4.0,   q_t_per_m: 16.0}
    - {x: 4.0,   q_t_per_m: 30.0}    # jump: aft peak bulkhead, superstructure starts
    - {x: 14.0,  q_t_per_m: 36.0}
    - {x: 16.0,  q_t_per_m: 36.0}
    - {x: 16.0,  q_t_per_m: 64.0}    # jump: engine room (main engine, auxiliaries)
    - {x: 28.0,  q_t_per_m: 64.0}
    - {x: 28.0,  q_t_per_m: 46.0}    # jump: forward end of the main engine
    - {x: 38.0,  q_t_per_m: 42.0}
    - {x: 38.0,  q_t_per_m: 34.0}    # jump: engine-room forward bulkhead
    - {x: 60.0,  q_t_per_m: 36.0}
    - {x: 80.0,  q_t_per_m: 35.0}
    - {x: 100.0, q_t_per_m: 30.0}
    - {x: 104.0, q_t_per_m: 28.0}
    - {x: 104.0, q_t_per_m: 36.0}    # jump: anchors, chain locker
    - {x: 112.0, q_t_per_m: 33.0}
    - {x: 118.0, q_t_per_m: 16.0}
    - {x: 120.0, q_t_per_m: 0.0}

# Deadweight items with the longitudinal moment about AP instead of the LCG
# (LCG = moment_tm / weight_t; giving lcg as well is allowed if it agrees).
weights:
  - {name: residual fuel, weight_t: 200, moment_tm: 14000, xa: 60.0, xb: 80.0}    # LCG 70 m
  - {name: fresh water,   weight_t: 60,  moment_tm: 6000,  xa: 94.0, xb: 106.0}   # LCG 100 m
  - {name: stores,        weight_t: 40,  moment_tm: 1200}                         # point load at 30 m
```

The curve of Listing 4.2 is converted into 12 trapezoidal items with a total of 4343.0 t and an LCG of 55.670 m; with the three items the docking weight is 4643.0 t at an LCG of 56.639 m, which the loader cross-checks against the declared values of the `ship` section.

| Key | Meaning | Unit |
|---|---|---|
| `weights[].name`, `weight_t` | name and weight of an item (negative for a deduction) | —, t |
| `weights[].lcg` or `moment_tm` | LCG, or moment about AP (LCG = moment / weight); both allowed if they agree | m, t·m |
| `weights[].xa`, `xb` | extent of a distributed item (both or neither) | m |
| `weights[].x` | position of a point load | m |
| `weight_curve.name`, `points[].x`, `points[].q_t_per_m` | piecewise-linear curve; a plain list of points is also accepted | —, m, t/m |

Table: Table 4.5 – Input keys of the weights and the weight curve.

### 4.7.2 Block plan

The blocks are defined only by the block plan `blocks.plan` (Section 3.4.4), because the number of blocks, their spacing and their individual strengths are exactly what a dock's block list provides and what the spacing design keeps or changes. The plan gives the block type, the position of the first block, the number of blocks and the spacing to the next block; optionally the strength and the initial gap of every block. A single spacing, strength or gap applies to all blocks, and a list gives one value per gap ($N-1$ values) or per block ($N$ values). An empty strength means the capacity of the block type. If `count` is omitted, it is taken from the length of a list. The function `loader._blocks` checks that the lists have exactly the required length, that every spacing is positive and that the count is a whole number of at least two, and labels the blocks B01, B02, … from aft to forward; an input with the former block groups is rejected with a message that points to `blocks.plan`. Table 4.6 lists the keys and Listing 4.3 the plan of the spacing example, with 190 t standard blocks and nine older blocks of 160 t under the engine room.

| Key | Meaning | Unit (default) |
|---|---|---|
| `blocks.types.<name>` | block type: contact width and length, layers (`thickness`, `E_MPa`) or `stiffness_kN_per_mm`, `capacity_t`, `allowable_pressure_MPa`, `crush_stress_MPa`, `post_crush_ratio` | m, MPa, t (ratio 0.05) |
| `blocks.plan.type` | block type of the blocks | — |
| `blocks.plan.start_x` | position of the first (aftmost) block | m |
| `blocks.plan.count` | number of blocks $N$ | — (length of a list) |
| `blocks.plan.spacing` | distance to the next block: one value or $N-1$ values | m |
| `blocks.plan.strength_t` | strength of every block: one value or $N$ values | t (type capacity) |
| `blocks.plan.gap_mm` | initial gap of every block (block top below the keel line): one value or $N$ values | mm (0) |

Table: Table 4.6 – Input keys of the block types and the block plan.

**Listing 4.3 – Block plan of the spacing example (`examples/spacing_case.yaml`).**

```yaml
  # Block plan: one spacing for all gaps or a list of count-1 values; one strength for all
  # blocks or a list of count values (aft to forward); optional gap_mm likewise.
  plan:
    type: keel
    start_x: 6.4                # m, first block
    count: 44
    spacing: 2.4                # m, distance to the next block
    strength_t: [               # t, strength of every block, B01 ... B44
      190, 190, 190, 160, 160, 160, 160, 160, 160, 160, 160,
      160, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190,
      190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190,
      190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190]
```

The loader evaluates the layer build-up of the block type as springs in series (Section 3.4.1); for the keel block of the examples this gives a stiffness of 2000 kN/mm. An initial gap belongs to its block, not to a position: when the spacing design moves a block, its gap moves with it.

### 4.7.3 Welding deformation

The keys of the welding-deformation section are summarised in Table 4.7 and Listing 4.4 shows the section of the welded-hull example (`examples/welding_case.yaml`), a hull built from six 20 m construction blocks.

| Key | Meaning | Unit |
|---|---|---|
| `welding_deformation.joints[].x`, `angle_mrad` | position and kink angle of an erection joint; positive lifts the hull ends (the angle may be omitted when it is identified) | m, mrad |
| `…zones[].xa`, `xb` | extent of a construction block | m |
| `…zones[].curvature_mm_per_m2` | shrinkage curvature given directly | mm/m² |
| `…zones[].deck_shrinkage_mm_per_m`, `bottom_shrinkage_mm_per_m`, `depth` | alternative: curvature $=(\varepsilon_{\text{deck}} - \varepsilon_{\text{bottom}})/D$ | mm/m, m |
| `…measured.x`, `lift_mm` | surveyed keel line, upward, arbitrary datum | m, mm |
| `…measured.file` | alternative: CSV file with columns `x_m`, `lift_mm` (path relative to the case file) | — |
| `…measured.interpolation` | `pchip`, `linear` or `spline` | — |
| `…measured.identify_joint_angles`, `regularization` | fit the joint angles to the survey; Tikhonov weight | —, — |

Table: Table 4.7 – Input keys of the welding-deformation section.

**Listing 4.4 – Welding deformation of the assembled hull in the welded-hull example (`examples/welding_case.yaml`).**

```yaml
welding_deformation:
  joints:                       # erection butt joints between construction blocks
    - {x: 20.0,  angle_mrad: 0.9}
    - {x: 40.0,  angle_mrad: 0.5}
    - {x: 60.0,  angle_mrad: 0.4}
    - {x: 80.0,  angle_mrad: 0.5}
    - {x: 100.0, angle_mrad: 0.6}
  zones:                        # curvature inside each construction block
    # curvature = (deck shrinkage − bottom shrinkage) / depth; or give curvature_mm_per_m2
    - {xa: 0.0,   xb: 20.0,  deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
    - {xa: 20.0,  xb: 40.0,  deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
    - {xa: 40.0,  xb: 60.0,  deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
    - {xa: 60.0,  xb: 80.0,  deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
    - {xa: 80.0,  xb: 100.0, deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
    - {xa: 100.0, xb: 120.0, deck_shrinkage_mm_per_m: 0.24, bottom_shrinkage_mm_per_m: 0.20, depth: 10.0}
```

Keel surveys are normally delivered as a list of stations, for example from a total station or a laser tracker. The loader therefore also reads the survey from a CSV file (`read_survey_csv`): a header row must contain an `x_m` (or `x`) and a `lift_mm` (or `lift`) column, lines beginning with `#` are comments, and malformed lines are reported with their line number. The function `inline_data_files` replaces the `file` key by the lists it contains before the case is built, so that the rest of the program, and the GUI, work with the same in-memory structure as for an inline survey; the core itself never reads files. A sample survey (`examples/welding_survey.csv`, 25 synthetic stations every 5 m) and a case that identifies the joint angles from it (`examples/welding_survey_case.yaml`) are supplied. A survey used without identification is taken as the complete keel shape: joints or zones given as well are then ignored with a warning. When the joint angles are identified, the known zone-shrinkage part is subtracted from the survey before the fit (Section 3.10).

### 4.7.4 Spacing design

The section `spacing_design` requests the automatic spacing design and holds its settings (Table 4.8). If the section is present and not disabled with `enabled: false`, `run_case` designs the spacing of the blocks of the plan. Listing 4.5 shows the settings of the spacing example. In the welded-hull example the same section is present but disabled, so that the as-built condition is analysed by default.

| Key | Meaning | Unit (default) |
|---|---|---|
| `enabled` | run the spacing design | — (true if the section is present) |
| `method` | support model of the design: `fem` or `rigid` | — (`fem`) |
| `target_utilisation` | largest admissible ratio $u$ of reaction to strength | — (1.0) |
| `min_spacing`, `max_spacing` | limits $s_{\min}$, $s_{\max}$ of the distance between neighbouring blocks | m (1.0; twice the mean spacing of the given plan) |
| `x_min`, `x_max` | allowed extent of the block line | m (first and last block of the plan) |
| `fix_end_blocks` | keep the first and last block in place | — (false) |
| `snap_to_structure` | final positions at the structure positions (floors) | — (true) |
| `grid` | rounding step when not snapping to structure | m (none) |
| `max_iterations` | iteration limit of the sequential linear programming | — (20) |

Table: Table 4.8 – Input keys of the spacing design.

**Listing 4.5 – Spacing-design settings of the spacing example (`examples/spacing_case.yaml`).**

```yaml
spacing_design:
  enabled: true
  method: fem
  target_utilisation: 1.0       # reaction / strength must not exceed this
  min_spacing: 0.8              # m
  max_spacing: 4.8              # m (default: twice the given mean spacing)
  x_min: 6.4                    # m, aft limit of the block line (default: first block)
  x_max: 109.6                  # m, forward limit (default: last block)
  fix_end_blocks: false
  snap_to_structure: true       # final positions under floors (structure positions)
  max_iterations: 20            # sequential-LP iterations
```

## 4.8 Solver implementation

### 4.8.1 Mesh and load model

The mesh is built by `weights.build_nodes`. Every *key point* becomes a node: weight-item boundaries (including the segment ends of the weight curve), hull-section boundaries, block positions and the hull ends. The intervals between key points are then subdivided so that no element is longer than the user-specified maximum (0.5 m by default). Because every item boundary is a node, the load on every element is exactly linear, and the consistent load vector integrates it without approximation (Section 3.3.6). Because every block is a node, its spring acts directly on a vertical degree of freedom. The example vessel produced 256 nodes, 255 elements and 512 degrees of freedom. Since block positions are key points, a new mesh is built for every candidate layout of the spacing design; this costs well under a millisecond and keeps the load integration exact for every candidate.

### 4.8.2 Vectorised assembly

The element stiffness matrices of all elements are computed at once as a three-dimensional NumPy array of shape (number of elements, 4, 4), using the Timoshenko element of Section 3.6. Euler–Bernoulli elements follow as the special case $\Phi = 0$. The global matrix is then assembled without a Python loop over elements (Listing 4.6). For element $e$ the global degrees of freedom are $2e, 2e+1, 2e+2, 2e+3$. Repeating and tiling this index array produces the row and column index of every one of the 16 entries per element. The entries are passed to SciPy in coordinate (COO) format, and converting to compressed sparse row (CSR) format sums the duplicate entries at shared nodes, which is exactly the assembly operation of the direct stiffness method [@przemieniecki1968; @cook2002].

**Listing 4.6 – Vectorised assembly of the global stiffness matrix and load vector (`solvers/fem.py`).**

```python
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
```

The function `np.add.at` is used for the load vector instead of ordinary fancy-index addition because it accumulates repeated indices correctly.

### 4.8.3 Springs, gaps and the linear solve

For a given set of block states, each block force is written in the linearised form

$$ R_i = k_{t,i}\,\delta_i + c_i, \qquad \delta_i = v_i - g_i \tag{4.2} $$

where $(k_{t,i}, c_i)$ equals $(0, 0)$ for an open block, $(k_i, 0)$ for an elastic block and $(k_{2,i}, (k_i - k_{2,i})\,\delta_{y,i})$ for a crushed block (Section 3.4.3). Moving the terms that depend on the unknown displacement to the left-hand side gives the system that is solved in each iteration:

$$ \left(\mathbf{K}_0 + \operatorname{diag}(k_{t,i})\right)\mathbf{u} = \mathbf{F}_0 + \sum_i \left(k_{t,i}\,g_i - c_i\right)\mathbf{e}_{v,i} \tag{4.3} $$

where $\mathbf{K}_0$ and $\mathbf{F}_0$ are the hull stiffness matrix and load vector from Listing 4.6, and $\mathbf{e}_{v,i}$ is the unit vector of the vertical degree of freedom at block $i$. The beam itself is free–free: without springs $\mathbf{K}_0$ is singular, with the two rigid-body modes as its null space. Listing 4.7 shows the implementation of Equation (4.3). The matrix is symmetric and banded, and it is solved with SciPy's sparse direct solver, which is based on SuperLU [@li2005superlu; @davis2006].

**Listing 4.7 – Insertion of block springs and gaps and the sparse solve (`solvers/fem.py`).**

```python
    def linear_solve(kt, c, gap):
        kdiag = np.zeros(ndof)
        np.add.at(kdiag, vdofs, kt)
        F = F0.copy()
        np.add.at(F, vdofs, kt * gap - c)
        u = spla.spsolve((K0 + sp.diags(kdiag)).tocsc(), F)
        if not np.all(np.isfinite(u)):
            raise UnstableSupportError("FEM system is singular")
        return u[vdofs] - gap, u
```

The function returns the block compressions $\delta_i$ and the full displacement vector. The hull matrix $\mathbf{K}_0$ is assembled only once per analysis, and each contact iteration merely adds a diagonal matrix. The gap vector is an explicit argument because the gap continuation of Section 4.8.5 must solve the same support model for a sequence of intermediate gap vectors.

### 4.8.4 State iteration

The linearisation of Equation (4.2) is computed by the method `BlockSprings.linearize` (Listing 4.8). The states are integer codes (open, elastic, crushed). For blocks without a crush limit, the crush deflection is infinite, and it is replaced by zero before the multiplication to avoid the undefined product $0 \cdot \infty$.

**Listing 4.8 – Linearisation of the piecewise-linear block law (`solvers/contact.py`).**

```python
    def linearize(self, states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(kt, c) such that R = kt·δ + c for the given states."""
        dy = np.where(np.isfinite(self.crush_deflection), self.crush_deflection, 0.0)
        kt = np.select([states == ELASTIC, states == CRUSHED], [self.k, self.k_crushed], 0.0)
        c = np.where(states == CRUSHED, (self.k - self.k_crushed) * dy, 0.0)
        return kt, c
```

The iteration (Listing 4.9) starts from given block states, normally all elastic. It solves the linear system, determines the state that each block's compression actually implies, and repeats until the assumed and implied states agree. A converged solution then satisfies equilibrium, compatibility and the block law exactly. The iteration is an active-set method for the complementarity problem formed by compression-only supports [@cottle1992; @wriggers2006] (Section 3.7). Three safeguards were implemented:

1. Before every solve, `check_support` verifies that at least two separated blocks have positive stiffness, and raises `UnstableSupportError` otherwise.
2. Every visited set of states is recorded. If a set would be visited again, which indicates cycling, the method switches to changing only the block with the largest violation of its assumed state per iteration.
3. The same one-at-a-time update is used if a full update would leave fewer than two separated blocks in contact. This happens with large, uneven gaps: from an all-elastic start, the first solution may show negative compression at almost every block.

**Listing 4.9 – Main loop of the block state iteration (`solvers/contact.py`, function `_iterate`).**

```python
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
```

The rigid-hull solver (`solvers/rigid.py`) passes a different `linear_solve` function to the same loop. It solves the $2 \times 2$ system of Section 3.5 in coordinates centred on the block centroid, which improves its conditioning. Sharing the state iteration means that lift-off, gaps, crushing and the continuation of the next section are handled identically in both methods, and that the two methods can be compared directly, as in Chapter 6.

### 4.8.5 Gap continuation

The safeguards of Section 4.8.4 did not suffice for all welding cases. The one-at-a-time update can itself cycle, which is a known weakness of simple active-set methods. During development this occurred for the rigid-hull method of the example vessel with joint kinks of 0.2 mrad at five erection joints: the direct iteration did not converge within 200 iterations (Section 5.4.3). The function `solve_contact` was therefore extended with a gap continuation (homotopy, Section 3.11), shown in Listing 4.10.

If the direct iteration fails, the problem is first solved with all gaps equal to the smallest gap $g_0$. This is a rigid shift of the block line, so every block starts in contact and the iteration converges as for a straight keel. The gaps are then increased along

$$ \mathbf{g}(\lambda) = g_0 + \lambda\,(\mathbf{g} - g_0), \qquad 0 \le \lambda \le 1 \tag{4.4} $$

in steps of $\lambda$. Each step starts from the contact states of the previous one, so only a few blocks change state per step. A failed step is halved; a successful one doubles the next step, up to 0.5. If the step falls below $1/4096$, the original error is reported.

**Listing 4.10 – Gap continuation in `solve_contact` (`solvers/contact.py`).**

```python
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
```

The continuation is a fallback. The direct iteration is always tried first, and for the welded-hull example it converged by itself, in 5 iterations for the rigid hull and 3 for the FEM. The continuation changes only the path to the solution: the final state set is accepted only when it is consistent with the full gap vector, so the result satisfies the same complementarity conditions as a directly converged solution. For the spacing design this robustness is essential, because the design analyses hundreds of layouts with welding gaps at positions that the user never chose, and a single failed analysis would abort the design.

### 4.8.6 Results

Both solvers return a `SupportResult` object. It contains the block positions, the reactions, the compressions, the final states, the number of iterations and the nodal displacements and rotations. It also provides two derived properties, the relative force residual $(\sum R_i - W)/W$ and the relative moment residual, which are reported in every summary (requirement N4).

## 4.9 Welding-deformation module

The module `deformation.py` implements the model of Section 3.9 for the initial keel shape of the welded hull. The two shape functions are given in Listing 4.11. A kink $\varphi_j$ at an erection joint contributes $\varphi_j\,|x - x_j|/2$ to the keel lift. A construction block $[x_a, x_b]$ with shrinkage curvature $\kappa$ contributes $\kappa\,\psi(x)$, where $\psi$ has unit second derivative inside the block, is zero aft of it and is straight forward of it.

**Listing 4.11 – Shape functions of a joint kink and of a shrinkage zone (`deformation.py`).**

```python
def kink_shape(x: np.ndarray, xj: float) -> np.ndarray:
    """Lift per radian of a joint kink at xj (symmetric form)."""
    return 0.5 * np.abs(np.asarray(x, dtype=float) - xj)


def zone_shape(x: np.ndarray, xa: float, xb: float) -> np.ndarray:
    """Lift per unit curvature of a zone [xa, xb]: z'' = 1 inside, 0 outside."""
    x = np.asarray(x, dtype=float)
    inside = (x - xa) ** 2 / 2.0
    after = (xb - xa) * (x - 0.5 * (xa + xb))
    return np.where(x < xa, 0.0, np.where(x <= xb, inside, after))
```

A surveyed keel line is interpolated with a monotonicity-preserving piecewise cubic Hermite interpolant (PCHIP) [@fritsch1980] by default; linear interpolation and a natural cubic spline are available as alternatives. Optionally the survey is not interpolated but used to identify the joint angles: `identify_joint_angles` solves the linear least-squares problem of Section 3.10 with NumPy, including a straight line for the unknown survey datum and an optional Tikhonov term [@hansen1998]. The function `effective_model` decides which model is used for the gaps — the given joints and zones, the survey as the total shape, or the zones plus the identified joint angles — and returns the identification with its rms and largest residual for the report.

The keel lift $z(x)$ is defined only up to a straight line, because a straight line is a rigid tilt that both support solvers absorb. The function `gaps_from_lift` converts the lift at the blocks into gaps by subtracting the least-squares line through the block positions and shifting the result so that the lowest point has zero gap:

$$ g_i^{\text{weld}} = z(x_i) - (\alpha + \beta x_i) - \min_j \left[z(x_j) - (\alpha + \beta x_j)\right] \tag{4.5} $$

where $\alpha + \beta x$ is the least-squares line. The normalisation makes the reported gaps easy to interpret ("end lift relative to the block line") but has no effect on the reactions (Section 5.5). The gaps depend on the block positions through $z(x_i)$ and through the line $\alpha + \beta x$. This is why the deformation model, not a list of gaps, is passed to the analysis: the function `analysis._springs` adds the welding gaps of the current block positions to the initial gaps of the blocks with `dataclasses.replace`, both for the given plan and for every candidate layout of the spacing design (Section 4.10.5). The welding deformation is thereby a constraint of the design, as required by Section 3.12.2.

## 4.10 Spacing-design module

The module `spacing.py` implements the automatic spacing design of Section 3.12 — the final result of the program. Given $N$ blocks in a fixed order with strengths $S_i$, it finds positions $x_1 < \dots < x_N$ within the allowed extent and the spacing limits that minimise the largest utilisation $U(\mathbf{x}) = \max_i R_i(\mathbf{x})/S_i$, and reports whether $U \le u$, the target utilisation. This section describes how the method was engineered; its algorithm is given in Section 3.12.8.

### 4.10.1 Structure of the module

The module has 509 lines and is organised by the phases of the method (Table 4.9). It imports only the data model and SciPy's `linprog`; in particular it imports no solver, which is the precondition for testing it with a substitute analysis.

| Part | Function or class | Responsibility |
|---|---|---|
| helpers | `utilisation`, `relabel` | $R_i/S_i$ (NaN without strength); sort and label B01, B02, … |
| | `_admissible` | one aft-to-forward sweep that enforces $s_{\min}$, $s_{\max}$, fixed ends and, optionally, a set of admissible positions |
| | `_Counter` | wraps the injected `evaluate`, counts analyses, returns $(\text{evaluation}, U_i)$ |
| | `_score` | ordering of layouts: largest utilisation, then $\sum_i U_i^8$ as tie-break |
| results | `SpacingStep`, `SpacingResult` | frozen record of one recorded step; result with initial and final blocks and evaluations, history, analysis count, diagnosis |
| phase 1 | `equidistribution_step` | strength-weighted equidistribution of the estimated load density (Section 3.12.3) |
| phase 2 | `_sensitivities`, `_slp_step` | finite-difference Jacobian $\partial U_i/\partial x_k$; assembly and solution of the linearised minimax LP (Section 3.12.4) |
| phase 3 | `_grid_neighbour`, `_round_to_grid`, `_local_moves` | rounding to admissible positions; neighbourhood of the local search (Section 3.12.5) |
| driver | `design_spacing` | input checks, limits, the three phases, diagnosis (Section 3.12.6) |

Table: Table 4.9 – Structure of the module `spacing.py`.

The interface is the single function

```text
design_spacing(blocks, evaluate, design, x_bounds, structure=None,
               equidistribution_iterations=15, fd_step=0.05, local_search_evaluations=None)
```

in which `evaluate(blocks)` is any function that returns an object with an attribute `reactions`. The type of `evaluate` is declared as `Evaluate = Callable[[Sequence[Block]], object]`; Python's structural typing makes the substitute analysis of the tests and the `SpacingEvaluation` of the program interchangeable without a common base class.

The wrapper `_Counter` (Listing 4.12) isolates the design from the representation of blocks. The phases work on a plain position vector $\mathbf{x}$; `_Counter` turns it into blocks by copying the given blocks with new positions, which preserves the order, the type, the strength and the initial gap of every block. It also counts the analyses, which is the machine-independent measure of the cost of the design (Section 5.9), and clips negative reactions, which cannot occur in a converged contact solution but may occur in a substitute analysis.

**Listing 4.12 – Wrapper of the injected analysis (`spacing.py`).**

```python
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
```

Positions are rounded to 1 µm, so that a block position is always a reproducible mesh node.

### 4.10.2 Limits and the admissible sweep

Before the phases start, `design_spacing` checks that every block has a strength and fixes the limits. The extent $[x_{lo}, x_{hi}]$ is taken from the design settings, from the given block line, or from the given end blocks if they are fixed. The maximum spacing defaults to twice the mean spacing of the given plan; an unlimited spacing would let the minimax optimum leave long unsupported spans amidships, where the block reactions are small but the hull girder spans long. Two infeasible combinations are rejected with an `InputError`: more blocks than fit into the extent at $s_{\min}$, and, with fixed end blocks, an extent that $N-1$ spacings of $s_{\max}$ cannot span. With free ends, the second case shortens the line instead.

All three phases need to map an arbitrary position vector to an admissible one. The function `_admissible` does this in one sweep from aft to forward with the end blocks fixed: block $i$ is kept in the window

$$ \max\left(x_{i-1} + s_{\min},\; x_N - (N-i)\,s_{\max}\right) \;\le\; x_i \;\le\; \min\left(x_{i-1} + s_{\max},\; x_N - (N-i)\,s_{\min}\right) \tag{4.6} $$

that its predecessor and the room needed by the remaining blocks leave for it, and moved to the nearest point of the window, or to the nearest admissible position in it when a set of positions (floors or a grid) is given. The window is never empty if the line length lies between $(N-1)s_{\min}$ and $(N-1)s_{\max}$, so the sweep succeeds in $O(N)$ operations without a search; with a discrete set of positions it may fail, and returns `None`.

### 4.10.3 Phases 1 and 2

Phase 1 starts from the given plan mapped proportionally onto the extent and applies the under-relaxed fixed-point iteration of `equidistribution_step` (Listing 4.13) up to 15 times. The reactions of the current layout are spread uniformly over the tributary cells of the blocks, which gives a piecewise-constant load density with cumulative $Q(x)$; inverting $Q$ at the cumulative interval weights $(S_i + S_{i+1})/2$ places the blocks so that the load between two neighbours is proportional to their mean strength (Section 3.12.3). Blocks that are out of contact are given 5 % of the mean load, so that the density never vanishes and $Q$ remains invertible. Each new layout is moved halfway towards the equidistributed one, made admissible and analysed; the best layout found is kept.

**Listing 4.13 – Strength-weighted equidistribution step (`spacing.py`).**

```python
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
```

Phase 2 is sequential linear programming with a trust region (Section 3.12.4). The function `_sensitivities` computes the Jacobian $J_{ik} = \partial U_i/\partial x_k$ by one-sided finite differences with a step of $h = 0.05$ m. For each free block the forward step is tried first and the backward step if the forward one would violate a limit; a block that cannot move in either direction keeps a zero column. One Jacobian therefore costs at most $N$ analyses. The function `_slp_step` assembles the linearised minimax problem as a dense linear programme in the variables $\Delta^{+}, \Delta^{-} \ge 0$ (the move of each free block, $\Delta = \Delta^{+} - \Delta^{-}$), $a_i \ge 0$ (the magnitude of the second difference of the new positions) and $t$:

$$ \begin{aligned} \min\;& t + \mu \sum_{i=1}^{N-2} a_i + \varepsilon \sum_k \left(\Delta_k^{+} + \Delta_k^{-}\right) \\ \text{s.t.}\;& U_i + \mathbf{J}_i\,\boldsymbol{\Delta} \le t, \qquad s_{\min} \le x'_{i+1} - x'_i \le s_{\max}, \\ & \left|x'_{i+2} - 2x'_{i+1} + x'_i\right| \le a_i, \qquad 0 \le \Delta_k^{\pm} \le \rho, \qquad x'_1 \ge x_{lo},\; x'_N \le x_{hi} \end{aligned} \tag{4.7} $$

with $\mathbf{x}' = \mathbf{x} + \boldsymbol{\Delta}$. The trust radius $\rho$ enters as a bound on the variables, and the extent only as a bound on the moves of the two end blocks, which keeps the constraint matrix small. The smoothing weight $\mu = 0.002/\bar{s}$, with $\bar{s}$ the mean spacing, trades a utilisation of 0.002 for one metre of spacing variation; it selects, among layouts with practically the same largest utilisation, the one whose spacing changes gradually, which is easier to lay out in the dock. The tiny weight $\varepsilon = 10^{-6}$ suppresses needless moves. For $N = 44$ the LP has 131 variables and 214 inequality rows and is solved by HiGHS through `scipy.optimize.linprog`.

The driver accepts a step only if the analysed layout has a better score than the current one; the trust radius then doubles, up to the mean spacing, and the Jacobian is recomputed. Otherwise the radius is halved, and the iteration stops when it falls below 5 mm, when the LP predicts an improvement of less than $10^{-5}$, or after `max_iterations` iterations (Listing 4.14). Because every accepted layout was analysed with the full support model, the linearisation affects only the efficiency of the search, never the correctness of its result (requirement N8).

**Listing 4.14 – Sequential linear programming with trust region (`spacing.py`, function `design_spacing`).**

```python
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
```

### 4.10.4 Phase 3: admissible positions and local search

Blocks must stand under floors (Section 3.12.5). The admissible positions are the structure positions when `snap_to_structure` is set, otherwise the multiples of `grid`, if given; without either, phase 3 is skipped. `_round_to_grid` moves free end blocks to the nearest admissible position inside the extent and places the interior blocks by the admissible sweep of Equation (4.6). Rounding a continuous optimum independently block by block could violate the minimum spacing; the sweep cannot.

The rounded layout is then improved by a best-improvement local search on the grid (Listing 4.15). The neighbourhood is generated by `_local_moves`: the governing blocks — those within 2 % of the largest utilisation, at most three — and their neighbours are each moved by one or two grid positions forward or backward, either alone or together with the next block. Only admissible layouts are returned, and duplicates are removed. All trials are analysed, the best one is accepted if it improves the score, and the search stops at a local optimum on the grid or when its budget of $25N$ analyses is used. Moving pairs of blocks matters because, at the minimum spacing, a single block often cannot move without violating $s_{\min}$, whereas a block and its neighbour can.

**Listing 4.15 – Best-improvement local search on the admissible positions (`spacing.py`, function `design_spacing`).**

```python
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
```

### 4.10.5 Result object and diagnosis

`design_spacing` returns a `SpacingResult` with the status `fits` or `does_not_fit`, a message, the initial and final blocks with their evaluations, the history of `SpacingStep` records (phase, iteration, largest and smallest utilisation, analyses so far), the number of analyses, the target, and the total strength and load. Its properties give the initial and final utilisations and the given and designed spacing ($N-1$ values). When the layout does not fit, the result carries the diagnosis of Section 3.12.6: the suggested number of blocks $\max(N+1, \lceil N U/u \rceil)$; whether the total strength is less than the load; whether an end block at the limit of the line governs, with the advice to allow a longer line; and, if blocks carry no load in the best layout, the statement that the gaps rather than the spacing govern and that more blocks may not help.

### 4.10.6 Integration into the analysis

`run_case` runs the spacing design after the analysis of the given plan (Listing 4.16). The extent defaults to the given block line and is always limited to the hull. The injected analysis is created by `spacing_evaluator` (Listing 4.17), a closure over the case, the welding-deformation model actually used (after identification) and the support method. For every candidate it builds the mesh and the load model, computes the springs with the welding gaps at the candidate positions, and solves the contact problem — exactly as `run_case` does for the given plan. The result is a frozen `SpacingEvaluation` that keeps the load model, the support result and the springs, so that the final layout needs no further analysis. The designed blocks are then stored as the additional support solution `fem+spacing` (or `rigid+spacing`) and checked by the strength module like any other solution, so that all hull girder and block checks are applied to the designed spacing. If the design does not fit, its message is added to the warnings.

**Listing 4.16 – Spacing design in `run_case` (`analysis.py`).**

```python
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
```

**Listing 4.17 – The analysis injected into the spacing design (`analysis.py`).**

```python
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
```

Because the designed layout has the same number of blocks as the plan but different positions, the report must not mix it with the solutions on the given plan. The property `AnalysisResult.block_supports` returns only the solutions on the given plan, which the block table of the summary, `blocks.csv` and the GUI's block table show; the designed spacing is reported separately (Section 4.11).

The cost of the design is the number of analyses times the cost of one analysis. For the 44-block examples one analysis — mesh, load model and FEM contact solution — took about 3.9 ms, and the designs needed between 78 and 501 analyses (Section 5.9). Every analysis is independent, so the finite differences of phase 2 and the trials of phase 3 could be evaluated in parallel; this was not necessary at the present problem sizes.

## 4.11 Post-processing and reporting

### 4.11.1 Strength curves and checks

The function `strength.compute_strength` converts a support result into shear force and bending moment curves. Because the load on each element is linear, $V(x)$ is quadratic and $M(x)$ cubic within each element, and both are evaluated exactly at eight points per element (Section 3.8). At every node the curves are stored twice, once just aft and once just forward of the node, so that the jumps caused by block reactions and point loads appear as vertical steps in the plots. The deck and keel stresses, the shear stress and the utilisation against the permissible envelopes are evaluated at the same points.

The checks are represented uniformly by a `Check` object with a name, a demand, a capacity, a utilisation and a location. A single helper function, `governing`, finds the location of highest utilisation for any combination of demand and capacity arrays; where no limit is defined, it reports the largest absolute demand with the status "no limit". Seven checks are produced for each support solution: hull bending moment, hull shear force, deck stress, keel stress, hull shear stress, block bearing pressure, and block load against the strength of each block (`Block.effective_strength`, the block's own strength or else the capacity of its type). Crushing appears as the block state and in the warnings. The function `run_case` issues warnings for:

- a sum of the weight items or an LCG that differs from the declared values by more than the tolerance;
- blocks that carry no load, listed by name (with a note when welding deformation is present);
- blocks that have crushed;
- blocks positioned away from transverse structure;
- closure errors of the shear and moment curves;
- a survey used as the total shape while joints or zones are also given;
- the welding end lift at both ends of the block line;
- a spacing design that does not fit, with its diagnosis.

### 4.11.2 Output files

The reporting module writes the files listed in Table 4.10; the command-line interface and the GUI's *Export results* command call the same functions. The text summary contains, in this order: the weight check (total weight and LCG, and with a weight curve its share and the number of items), the mesh and block line; a table of all support solutions with the sum, maximum and minimum of the reactions, open and crushed blocks, iterations and the two residuals; the block table of the given plan; the hull girder and block checks of every support solution with an overall verdict; the section **Welding deformation of the assembled hull** (model sources, identified angles with their fit residual, end lift and the blocks without contact as built) when welding deformation is present; the section **Block spacing design** when the spacing design ran; and the warnings. The spacing section gives the status and message, the block line and the range of the spacing before and after, the largest and mean ratio of reaction to strength before and after, the total strength against the load, the phases with their largest utilisation and the number of analyses, one row per block (position, spacing to the next block, strength, reaction and ratio), and, for a welded hull, the largest recomputed welding gap and the blocks without contact in the designed layout.

| File | Content | Written |
|---|---|---|
| `summary.txt` | text summary described above | always |
| `blocks.csv` | per block of the given plan: position, total and welding gap, stiffness, strength, reaction for each solution, pressure, utilisations, compression, state, offset from structure | always |
| `curves_rigid.csv`, `curves_fem.csv`, `curves_fem_spacing.csv` | $x$, $V$, $M$, deck, keel and shear stress, utilisations — one file per support solution | per solution |
| `spacing.csv` | per block: designed position and spacing to the next block, given position and spacing, strength, reaction and utilisation designed and given | with spacing design |
| `spacing_history.csv` | per recorded step: phase, iteration, analyses so far, largest and smallest utilisation | with spacing design |
| `results.png` | load and reactions, settlement, shear force, bending moment | unless `--no-plots` |
| `block_forces.png` | weight distribution and force on every block (Figure 4.2) | unless `--no-plots` |
| `welding_data.png` | keel shape with its components and survey points, joint angles and zone curvature, gap at every block (Section 6.7) | with welding deformation |
| `spacing_design.png` | given and designed spacing along the hull, reactions against strengths, convergence over the analyses (Section 6.8) | with spacing design |

Table: Table 4.10 – Files written by the command-line interface.

All plots are drawn into Matplotlib `Figure` objects without the `pyplot` state machine. The same drawing functions (`draw_results`, `draw_block_forces`, `draw_welding_data`, `draw_spacing`) therefore save PNG files from the command line and draw into the GUI's embedded canvases.

The block-force chart (Figure 4.2) answers the most common practical question, "what is the force on each block?", in a single picture. The upper panel shows a schematic hull with the blocks, the total weight and the LCG. The lower panel shows one bar per block, labelled with its force in tonnes and its block number. The bars are coloured by utilisation (green up to 90 %, orange from 90 % to 100 %, red above 100 % of the governing limit, grey for blocks without contact). The other support method is shown as markers for comparison, together with the block strength and the average load $W/N$. The spacing-design chart has three panels: the given and the designed spacing as bars along the block line, labelled where runs of equal spacing are long enough; the reactions of the designed layout as bars coloured by $R_i/S_i$ against the stepped strength line, with the reactions of the given plan as markers; and the largest utilisation over the number of analyses, coloured by phase, with the target.

![Figure 4.2: Block-force chart produced by the program for the example vessel (elastic hull FEM). Top: weight distribution w(x) with W, LCG and the block positions; middle: hull schematic with the blocks; bottom: force on every block, bars coloured by utilisation, diamonds showing the rigid-hull result for comparison.](figures/fig_block_forces_example.png)

## 4.12 Command-line interface

The command-line interface is invoked as `shipsupport run CASE.yaml`, or as `python -m shipsupport run CASE.yaml` without installation. Its options are given in Table 4.11. It prints the text summary to the console and writes the files of Table 4.10 to the output directory.

| Option | Effect |
|---|---|
| `-o OUT`, `--out OUT` | output directory (default: `<case name>_results` next to the case file) |
| `--methods {rigid,fem} …` | override the methods selected in the case file |
| `--no-plots` | skip the PNG charts |
| `--strict` | return exit code 2 if any check of any support solution fails or the spacing design does not fit |

Table: Table 4.11 – Options of the command-line interface.

The exit code makes the program usable in scripts and automated workflows. It is 0 after a successful run, 1 if the input could not be read or was invalid, or if the solver failed, and 2 if `--strict` was given and the result is not acceptable. Under `--strict` the result is not acceptable if any support solution — including the solutions on the given plan — has a failed check, or if a requested spacing design does not fit. A designed spacing that passes all checks therefore does not mask a given plan that fails: the user is told that the plan in the file must be replaced by the designed one. The behaviour was confirmed on the example cases:

- **Example vessel:** a normal run returned 0 and wrote six files. A run with `--strict` returned 2, because the FEM block bearing pressure and block strength checks fail. A run with `--strict --methods rigid` returned 0, because all rigid-hull checks pass. A run with a non-existent input file returned 1.
- **Spacing and weight-curve examples:** a normal run returned 0 and wrote ten files, among them `spacing.csv`, `spacing_history.csv`, `spacing_design.png` and `curves_fem_spacing.csv`. With `--strict` the spacing example returned 2, because the given plan fails the block checks, although the summary reports "ALL CHECKS PASSED" for the designed spacing `fem+spacing`.
- **Welded hull and survey case:** a normal run returned 0 and wrote seven files, including `welding_data.png`; with `--strict` the welded hull returned 2.

## 4.13 Graphical user interface

### 4.13.1 Overall design

The graphical user interface was built with PyQt5 [@pyqt5]. The main window is divided by a movable splitter: the left side holds the input editor, the right side the results. Its appearance is defined in one place, the module `style.py`. It provides a light palette with one accent colour, an application stylesheet (tabs, group boxes, tables with alternating rows and bold headers, rounded inputs, and primary buttons for the two analysis commands), the program icon drawn in code, and the colours of the verdict badge in the status bar. The charts share one Matplotlib style, set in the reporting module, so that the GUI, the exported PNG files and the figures of this thesis look alike. They use Matplotlib's constrained layout, which is recomputed whenever the canvas is resized; below a readable minimum size, a chart scrolls instead of being squeezed. A menu and a toolbar provide the file operations (new, open, save, save as, export results) and the two analysis commands: **Run analysis** (F5) analyses the case as it is, and **Auto-calculate spacing** (Ctrl+F5) runs the analysis with the spacing design and writes the designed spacing into the block plan. An *Open example* submenu loads any of the five supplied cases as an untitled copy, so that a shipped example cannot be overwritten by accident. The action *Add sample welding data* fills the welding input of the current case with the sample deformation of the welded hull, mapped onto the extent of the current hull, and runs the analysis. In the same way, *Add sample weight distribution curve* (also a button *Load sample curve* on the Weights tab) replaces the weights, after a confirmation, by a sample lightship curve of 18 points with five jumps and three deadweight items given by their moments (the data of `weight_curve_case.yaml`, 4643.0 t at LCG 56.639 m), mapped onto the current hull with the intensities scaled so that every weight is kept; it also sets the declared weight and LCG to the new totals and runs the analysis. The weight distribution is drawn as the top panel of the Block forces and Spacing design charts: the total distributed load w(x), the weight curve on its own (dashed), point loads, W and LCG, and the block positions (the given and the designed ones in the spacing chart). The status bar shows the overall verdict after each run, extended by the outcome of the spacing design when it ran, and reports when the results no longer match the input. For the example vessel it reads "CHECKS FAILED — max R 254.6 t at x = 6.40 m (fem); 0 warning(s)" (Figure 4.4).

### 4.13.2 Input editor

The input editor (class `CaseEditor`) is a set of tabs that corresponds to the sections of the YAML file (Table 4.12). It has two methods. `load_dict` fills all widgets from a case dictionary, and `to_dict` collects a case dictionary from the widgets. The GUI therefore works with exactly the same data structure as the file format: saving a case writes this dictionary as YAML, and running an analysis passes it to the calculation core.

| Tab | Content |
|---|---|
| Ship | name, $L_{pp}$, hull extent, declared weight and LCG |
| Weights | weight items with LCG or moment and optional extent; weight curve (name and points); live total weight and LCG of items and curve together (Figure 4.3) |
| Hull sections | hull girder sections ($I$, $E$, shear area, $G$, section moduli, optional $EI$ and $\kappa GA$) |
| Blocks | block types (layers written as `thickness@E_MPa; …`) and the block plan editor with the *Auto-calculate spacing* button (Figure 4.4) |
| Structure & limits | positions of transverse structure (start plus two of end, spacing, count, or an explicit list); permissible bending moment and shear force tables |
| Welding deformation | erection joints, construction-block shrinkage, surveyed keel line with an *Import CSV…* button, interpolation and identification options (Figure 4.7) |
| Settings | methods, mesh size, shear deformation, allowable stresses, tolerances, iteration limit |

Table: Table 4.12 – Tabs of the input editor.

Two reusable widget classes were developed for the editor. `RecordTable` edits a list of records, such as weight items, curve points, erection joints or survey points, in a table with buttons to add, duplicate and remove rows. Its table supports copying and pasting tab-separated blocks with Ctrl+C and Ctrl+V, so that data can be pasted directly from a spreadsheet, which is how weight lists, block lists and keel surveys are usually available. Empty cells are omitted from the records, so the loader's default values apply. `FormEditor` edits single values and shows each default value as greyed placeholder text. Unreadable entries are reported with their table, row and column.

The Weights tab (Figure 4.3) shows the weight items and the weight curve one above the other, and below them a live summary, for example "Σ weight = 4643.0 t, LCG = 56.639 m (3 items; curve 4343.0 t, LCG 55.670 m (18 points))" for the weight-curve example. The summary is computed exactly for the trapezoidal segments of the curve (function `curve_totals`), and an item given by its moment contributes that moment, so the user sees the effect of every edit on the docking weight and LCG before running an analysis.

![Figure 4.3: Weights tab of the input editor with the example vessel's weight items, the (empty) weight curve and the live total weight and LCG.](figures/fig_gui_weights.png)

### 4.13.3 Block plan editor and the automatic spacing design

The Blocks tab holds the block types and, below them, the block plan editor (class `BlockPlanEditor`, module `plan_panel.py`; Figure 4.4). It has a spin box for the number of blocks, a field for the position of the first block, a selector for the block type, and fields with buttons that set all spacings or all strengths at once. The table has one row per block with the columns Block, x, Spacing to next, Strength, Gap, Reaction and Reaction / strength. The spacing, the strength and the gap are edited directly or pasted from a spreadsheet; an empty strength means the capacity of the type. The label, the position — the cumulative sum of the spacings — and the two result columns are read-only and shaded. Changing the count adds rows at the forward end, which copy the last spacing and strength, or removes rows there. A line below the table summarises the plan, for example "44 blocks, block line 6.40 – 109.60 m (103.20 m), spacing 2.40 – 2.40 m, 44 block(s) with the type capacity as strength". The method `plan` returns the `blocks.plan` dictionary of Section 4.7.2, with a single value instead of a list where all values are equal.

The table is also the main view of the result. After every analysis the window passes the positions, reactions and ratios of reaction to strength to `show_reactions`, which fills the two result columns, coloured green, orange or red at 90 % and 100 % of the target, provided the result belongs to the plan as it is in the table. A result of an older plan is therefore never shown next to a changed plan.

![Figure 4.4: Blocks tab with the block type and the block plan editor (left), showing the reaction and the ratio of reaction to strength of every block of the example vessel, and the Block reactions result tab (right).](figures/fig_gui_blocks.png)

The **Auto-calculate spacing** button, the toolbar button and Ctrl+F5 all call `MainWindow.auto_spacing`, which carries out the following sequence:

1. The Blocks tab is shown and the status line under the button reports that the spacing is being designed.
2. `run_analysis(spacing_design=True)` collects the case dictionary from the editor, removes any `spacing_design` section and inserts the settings returned by `CaseEditor.auto_spacing_design` (Listing 4.18): the settings of the loaded case, always enabled, or — if the case has none — default settings with the floor spacing of the Structure tab as the minimum spacing.
3. Both analysis commands and the button are disabled, and the analysis runs in the background thread of Section 4.13.4.
4. When the result arrives, `show_result` displays it in all result tabs, and `_apply_spacing` writes the designed blocks back into the plan: the first position, the $N-1$ designed spacings, and the unchanged type, strengths and gaps. The status line reports "Fits" or "Does not fit", the largest ratio of reaction to strength before and after, the range of the designed spacing and the number of analyses, and, if the design does not fit, the diagnosis.
5. The designed plan is an ordinary input: it can be edited further, analysed with F5 and saved.

**Listing 4.18 – Settings of the spacing design for the Auto-calculate button (`shipsupport_gui/input_panels.py`).**

```python
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

The spacing-design settings of a loaded case are not edited in the GUI; they are kept by the editor and written back unchanged when the case is saved. F5 deliberately runs the analysis *without* the spacing design, even if the loaded case contains an enabled `spacing_design` section. This separates the two questions that the user asks — "what are the reactions of my plan?" and "what spacing should the plan have?" — and prevents a plain analysis from rewriting the plan.

### 4.13.4 Background analysis

The analysis runs in a background thread so that the window remains responsive, which matters for the spacing design, which takes up to a few seconds. The class `AnalysisThread` (Listing 4.19) receives the case dictionary and runs the same pipeline as the command-line interface. It reports the result or an error message to the main window through Qt signals. Signals emitted from a worker thread are delivered to the main thread through Qt's event loop, so all widget updates take place in the main thread, as Qt requires.

**Listing 4.19 – Background analysis thread of the GUI (`shipsupport_gui/worker.py`).**

```python
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

### 4.13.5 Result views

The results side has the eight tabs of Table 4.13. The plots are embedded Matplotlib canvases with a navigation toolbar for zooming, panning and saving. In the tables, utilisations above 1.0 are highlighted in red, utilisations between 0.9 and 1.0 in orange, and blocks without contact in grey. If a case contains no welding deformation or no spacing design, the corresponding tab explains this and offers a button that adds the sample welding data, or starts the automatic spacing design, respectively.

| Tab | Content |
|---|---|
| Summary | the text summary of Section 4.11.2 |
| Block forces | the chart of Figure 4.2: weight distribution on top, then the force on every block (Figure 4.5) |
| Plots | load and reactions, settlement, shear force, bending moment |
| Block reactions | one row per block of the given plan, with a welding-gap column when present (Figure 4.4) |
| Checks | all checks of all support solutions, including `fem+spacing` (Figure 4.6) |
| Welding data | the welding data chart (Figure 4.7) |
| Spacing design | the spacing design chart (Section 6.8) |
| Warnings | the warnings, with their number in the tab title |

Table: Table 4.13 – Tabs of the result view.

![Figure 4.5: Main window after an analysis of the example vessel, with the Blocks input tab on the left and the Block forces result tab on the right.](figures/fig_gui_block_forces.png)

![Figure 4.6: Settings input tab (left) and Checks result tab (right), listing every hull girder and block check for each support solution with utilisation and status.](figures/fig_gui_checks.png)

![Figure 4.7: Welding deformation input tab (left) and Welding data result tab (right) for the welded-hull example.](figures/fig_gui_welding.png)

## 4.14 Requirements traceability

Table 4.14 traces every requirement of Table 4.1 to the modules that realise it and to the test modules that verify it (requirement N9). The test modules are described in Chapter 5.

| Requirement | Realised by | Verified by |
|---|---|---|
| F1 weight items with LCG or moment | `model.WeightItem`, `loader._weight`, `weights.py` | `test_weights.py`, `test_loader.py` |
| F2 weight curve | `model.WeightCurve`, `loader._weight_curve` | `test_weights.py`, `test_loader.py`, `test_example_case.py` |
| F3 block plan | `loader._blocks`, `model.Block` | `test_loader.py`, `test_gui.py` |
| F4 rigid-hull method | `solvers/rigid.py` | `test_rigid.py`, `test_fem.py` (rigid limit) |
| F5 elastic-hull FEM | `solvers/fem.py` | `test_fem.py` |
| F6 compression-only, crushable blocks; large gaps | `solvers/contact.py` | `test_rigid.py`, `test_fem.py`, `test_deformation.py` |
| F7 hull girder checks | `strength.py` | `test_strength.py` |
| F8 block checks against block strength | `strength.py`, `model.Block.effective_strength` | `test_strength.py`, `test_example_case.py`, `test_spacing.py` |
| F9 welding-deformation model | `deformation.py`, `analysis._springs` | `test_deformation.py` |
| F10 survey import, interpolation, identification | `deformation.py`, `loader.read_survey_csv`, GUI import | `test_deformation.py`, `test_example_case.py`, `test_gui.py` |
| F11 automatic spacing design | `spacing.py`, `analysis.spacing_evaluator`, `analysis.run_case` | `test_spacing.py` |
| F12 reporting | `report.py` | `test_example_case.py`, `test_spacing.py` |
| F13 command-line interface | `cli.py` | `test_example_case.py`, `test_spacing.py` |
| F14 GUI with block plan editor and Auto-calculate spacing | `shipsupport_gui` (`plan_panel.py`, `main_window.py`) | `test_gui.py` |
| N1 SI units | `units.py`, `loader.py` | `test_loader.py`, `test_deformation.py` (unit conversion), all benchmarks |
| N2 input validation | `model.py`, `loader.py` | `test_weights.py`, `test_loader.py`, `test_deformation.py`, `test_spacing.py`, `test_gui.py` |
| N3 core independent of GUI | package structure | `test_gui.py` (identical results) |
| N4 run-time self-checks | `solvers/result.py`, `strength.py` | `test_example_case.py`, `test_spacing.py` |
| N5 automated tests | `tests/` | 131 tests |
| N6 performance | sparse solvers, vectorised assembly, analysis count of the design | measured (Section 5.9) |
| N7 portability | Python ecosystem | — |
| N8 designed spacing from full analyses | `spacing._Counter`, `analysis.spacing_evaluator` | `test_spacing.py` (designed layout re-analysed from scratch) |
| N9 traceability | this table | — |
| N10 maintainability | module structure | code metrics (Section 4.15) |
| N11 testability of the design | injected `evaluate` | `test_spacing.py` (substitute analyses) |

Table: Table 4.14 – Requirements traceability matrix: requirements, realising modules and verifying tests.

## 4.15 Code metrics

Table 4.15 lists the size of each module, measured as total lines and as lines that are neither blank nor full-line comments (docstrings are counted as code). The calculation core comprises 3283 lines, of which 2736 are non-blank, non-comment lines. The GUI comprises 2161 lines (1876). The test suite comprises 1552 lines (1236) in ten files — nine test modules and a helper module — and contains 131 tests (Chapter 5).

| Package / module | Lines | Non-blank, non-comment |
|---|---:|---:|
| `shipsupport/report.py` | 713 | 628 |
| `shipsupport/spacing.py` | 509 | 429 |
| `shipsupport/model.py` | 461 | 352 |
| `shipsupport/loader.py` | 420 | 367 |
| `shipsupport/strength.py` | 205 | 176 |
| `shipsupport/analysis.py` | 195 | 158 |
| `shipsupport/deformation.py` | 169 | 137 |
| `shipsupport/solvers/contact.py` | 169 | 134 |
| `shipsupport/weights.py` | 118 | 92 |
| `shipsupport/solvers/fem.py` | 116 | 97 |
| `shipsupport/cli.py` | 71 | 57 |
| `shipsupport/solvers/rigid.py` | 54 | 44 |
| `shipsupport/solvers/result.py` | 42 | 33 |
| `shipsupport/units.py` | 14 | 10 |
| other (`__init__`, `__main__`, `solvers/__init__`) | 27 | 22 |
| **Calculation core total** | **3283** | **2736** |
| `shipsupport_gui/input_panels.py` | 601 | 524 |
| `shipsupport_gui/main_window.py` | 444 | 389 |
| `shipsupport_gui/plan_panel.py` | 343 | 305 |
| `shipsupport_gui/widgets.py` | 296 | 256 |
| `shipsupport_gui/results_panels.py` | 252 | 220 |
| `shipsupport_gui/style.py` | 156 | 133 |
| `shipsupport_gui/worker.py` | 34 | 24 |
| `shipsupport_gui/app.py` | 31 | 22 |
| other (`__init__`, `__main__`) | 4 | 3 |
| **GUI total** | **2161** | **1876** |
| `tests/test_gui.py` | 418 | 348 |
| `tests/test_spacing.py` | 239 | 192 |
| `tests/test_loader.py` | 196 | 153 |
| `tests/test_deformation.py` | 184 | 144 |
| `tests/test_example_case.py` | 129 | 103 |
| `tests/test_fem.py` | 118 | 89 |
| `tests/test_weights.py` | 114 | 90 |
| `tests/test_rigid.py` | 68 | 51 |
| `tests/test_strength.py` | 51 | 42 |
| `tests/helpers.py` | 35 | 24 |
| **Tests total** (10 files, 131 tests) | **1552** | **1236** |

Table: Table 4.15 – Size of the modules of the program.

The numerical heart of the program is small. The two solvers with the contact iteration and the gap continuation amount to 381 lines, the welding-deformation model to 169 lines and the load model to 118 lines. The spacing design, the final result of the program, is with 509 lines the largest numerical module; its module docstring states the problem and the three phases, so that the code can be read against Section 3.12. The largest module is the reporting module, and the data model and the loader together account for 881 lines, reflecting the effort spent on clear presentation of results and on input validation with precise messages. The test code amounts to 47 % of the core, and the 239 lines of `test_spacing.py` verify a module of 509 lines.

## 4.16 Summary

The models of Chapter 3 were implemented in Python as a layered calculation package and a separate PyQt5 user interface, developed in tested increments. The calculation core uses validated data classes with SI units, frozen for all physical objects, and validates all input once, at the boundary. Weight items given by their LCG or their moment and a weight-distribution curve, whose segments are converted exactly into trapezoidal items, feed one load model on a mesh that contains every key point, so that the load integration is exact. The blocks are defined by a block plan with the number of blocks and the spacing, strength and initial gap of every block. Sparse finite element assembly is vectorised, and one active-set iteration for compression-only, optionally crushing blocks, with a gap continuation for large welding gaps, is shared by the rigid-hull and the elastic-hull methods through a common `linear_solve` interface.

The module `spacing.py` designs the spacing of the blocks of the plan from their strengths. It combines strength-weighted equidistribution, sequential linear programming with a trust region and a smoothing term, and an admissible rounding to the floors with a best-improvement local search, and it diagnoses why a plan cannot fit. The analysis of a candidate layout is injected as a function; in the program this function rebuilds the mesh and recomputes the welding gaps at the candidate positions, so that every designed spacing is the result of a full analysis and welding deformation acts as a constraint, while in the tests a substitute analysis with a known optimum takes its place.

The command-line interface supports scripted use through meaningful exit codes. The GUI offers spreadsheet-compatible tables, a live weight summary, a block plan editor that also shows the reaction and the ratio of reaction to strength of every block, and an *Auto-calculate spacing* button that runs the design in a background thread and writes the designed spacing back into the plan. A traceability matrix links every requirement to its modules and tests. The correctness of the implementation is examined in the next chapter.
