# 5 Verification

A numerical program is useful only if its results can be trusted, and a design program even more so: its final result, the spacing of the keel blocks, is laid out in the dock without the engineer having seen the hundreds of analyses behind it. This chapter documents the evidence that the program described in Chapter 4 solves the equations of Chapter 3 correctly. It first defines the verification strategy and gives an overview of the automated test suite. It then presents the analytical benchmark cases for the load, support and strength models, including the weight curve and weight items given by their moment, followed by convergence studies. Separate sections verify the welding-deformation model with the identification from a survey, and the automatic spacing design. The run-time self-checks built into every analysis, the consistency between the graphical user interface and the calculation core, and the performance of the program follow. The chapter closes with the limitations of the verification, in particular the absence of validation against measurements.

## 5.1 Verification strategy

Following @oberkampf2010, three activities that are often confused are distinguished:

- **Code verification** establishes that the software implements the chosen mathematical model without programming errors. Its main tool is comparison with exact analytical solutions of problems that exercise the relevant parts of the code.
- **Solution verification** estimates the numerical error of a particular computed solution, for example discretisation and round-off error, and confirms that the solution has converged.
- **Validation** establishes to what extent the mathematical model represents physical reality, by comparison with experimental or full-scale measurements.

This chapter deals with the first two activities. The models of Chapter 3 contain assumptions — the hull treated as a one-dimensional beam, linear-elastic or bilinear block behaviour, a kinematic description of welding deformation and no transverse load sharing — whose adequacy can only be assessed by validation. No measured block loads or keel surveys were available for this work, so validation remains outside its scope (Section 5.10).

The verification strategy had six elements:

1. **Analytical benchmarks.** Small problems with known closed-form solutions were solved with the program and the results compared with the closed form. Together, the benchmarks exercise every computational module: load integration from weight items and from the weight curve, the rigid-hull solver, the FEM solver with Euler–Bernoulli and Timoshenko elements, the contact iteration with gaps and crushing, the shear-force and bending-moment integration, and the shape functions of the welding-deformation model.
2. **Convergence and limit behaviour.** Properties that follow from the theory were tested: the FEM solution must approach the rigid-hull solution as the hull bending stiffness tends to infinity; the reactions must not depend on the mesh density; the contact solution for large gaps must satisfy the complementarity conditions whichever path the iteration takes.
3. **Invariance and equivalence tests.** Properties that must hold exactly were tested where no closed form exists: the invariance of the reactions to a straight-line tilt of the keel, the equivalence of a block with an unclosable gap and a removed block, and the equivalence of the example's weight items with the piecewise-linear curve that is their sum.
4. **Verification of the design algorithm against known optima.** The spacing design was tested with substitute analyses whose optimal layout is known in closed form, and it was checked on the full model that the designed layout, analysed again from scratch, gives the same reactions — the designed spacing is never the output of a surrogate.
5. **Run-time self-checks.** Every analysis computes equilibrium residuals and the closure error of the shear and moment curves, and reports them with the results.
6. **Automation.** All benchmarks and consistency checks were coded as automated tests with pytest [@pytest]. They are repeated after every change to the program, which protects the verified state against regressions.

The accuracy of a result $y$ compared with an exact reference value $y_{\text{ref}}$ is expressed throughout as the relative error

$$ \varepsilon = \frac{\left| y - y_{\text{ref}} \right|}{\left| y_{\text{ref}} \right|} \tag{5.1} $$

For calculations in IEEE double precision, the machine epsilon is approximately $2.2 \times 10^{-16}$. A relative error of order $10^{-16}$ to $10^{-12}$ therefore means agreement to the limits of floating-point arithmetic. All numbers in this chapter were computed with the final version of the program.

## 5.2 Test suite overview

The test suite consists of nine test modules and a helper module, and contains 131 tests in total. The count was determined with `pytest --collect-only`; parametrised tests count once per parameter set. Table 5.1 summarises what each module verifies. All 131 tests passed on the final version of the program. The GUI tests run Qt in its off-screen mode, so they need no display; if PyQt5 is not installed, they are skipped automatically.

| Test module | Tests | What is verified |
|---|---:|---|
| `test_weights.py` | 12 | trapezoid end intensities preserving weight and LCG; triangle limits; rejection of an LCG outside the middle third; exact integration of 31 random items; mesh construction and node lookup; weight-curve segments as exact trapezoids (weight, LCG, $\int q\,dx$, $\int xq\,dx$ and the intensity itself); a curve with a jump; five invalid curves, each reported with the offending point |
| `test_loader.py` | 33 | structure positions by the three rules and by an explicit list, five invalid specifications; block plan with lists and single values (positions, labels, strength default and override, gaps); count taken from a list; nine invalid plans, each with a specific message; rejection of the former block groups; reactions after changing count and spacing; LCG from the moment of an item and of a point load; moment and LCG must agree; four invalid weight items; weight curve read as items (also as a plain list and without items); four invalid curve points named in the message; rejection of a case without any weights |
| `test_rigid.py` | 6 | two-block determinate case; equal blocks under symmetric load; eccentric closed form; lift-off equivalence; crushing closed form; instability detection |
| `test_fem.py` | 10 | Euler–Bernoulli and Timoshenko midspan deflection; shear-deformation switch; two-span reactions; convergence to the rigid method; load concentration of a flexible hull; gap equivalence; partially closed gap; mesh independence; section coverage check |
| `test_strength.py` | 3 | simply supported shear and moment; hogging over supports with overhangs; stresses and utilisations of all hull girder checks |
| `test_deformation.py` | 18 | kink slope change; zone curvature; gaps and reactions invariant to a keel tilt; end lift unloading the first and last blocks; continuation with large gaps; exact and noisy identification of joint angles; interpolation through the survey points (three methods); unit conversion of shrinkage; invalid input; survey CSV and the sample survey case; identification with known zone shrinkage subtracted; survey used directly as the total shape; rejection of a malformed survey file; welding example with end lift |
| `test_spacing.py` | 13 | closed-form minimax optimum for a uniform load; smaller spacing for weaker blocks with the order of strengths kept; equidistribution step towards a load peak; spacing limits, extent and grid; fixed end blocks; report when the strength is not enough; invalid problems; spacing example within strength, under floors, in equilibrium, with command-line outputs; initial gaps moving with their blocks; welding gaps recomputed at the designed positions and the designed layout re-analysed from scratch; outcome on the welded hull; moderate welding deformation absorbed; severe welding deformation diagnosed as gap-governed |
| `test_example_case.py` | 8 | equilibrium, non-negative reactions and closure for the example vessel; summary and command-line run with exactly the expected output files; command-line runs of the welded hull and the survey case with the welding data chart; weight-curve example; the example's items as one weight curve giving identical reactions; labelling of blocks created in code; rejection of a negative weight tolerance |
| `test_gui.py` | 28 | layer text format; input and result tabs exactly as built; every example round-tripping through the editor (5 cases); identical results via the editor; Weights tab with moments and the weight curve, live Σ; Blocks tab with plan editing (count, set all, single cells, empty strength, gap, read-only x, bad cell); spacing-design settings kept and saved; F5 not running the spacing design; *Auto-calculate spacing* writing the design into the plan; default minimum spacing from the floors; survey case with its CSV points; sample welding data and its mapping onto the hull; unresolved survey file; error location; point-load rule; result display and save; every example opened and run; background thread and its error reporting; sample weight curve identical to `weight_curve_case.yaml`, loaded and analysed with the same reactions, and mapped onto a 150 m hull with the weights kept; weight-distribution panel on top of the Block forces and Spacing design charts |
| **Total** | **131** | |

Table: Table 5.1 – Modules of the automated test suite, number of tests and verified properties.

The complete test suite ran in 33 s on the workstation of Section 5.9.1, much of it spent in the tests that perform complete spacing designs, directly and through the GUI's background thread. It was run after every change during development, and the benchmarks of the following sections are part of it.

## 5.3 Analytical benchmark cases

The benchmark results are collected in Table 5.2 at the end of this section. The individual cases are described first.

### 5.3.1 Trapezoidal load integration

A single weight item of 1000 N spread over $[2, 8]$ m with its centre at 4.5 m was converted into a trapezoid. With $L = 6$ m and $\bar{x} = 2.5$ m, the end intensities of Section 3.3.2,

$$ w_a = \frac{2W}{L}\left(2 - \frac{3\bar{x}}{L}\right), \qquad w_b = \frac{2W}{L}\left(\frac{3\bar{x}}{L} - 1\right) \tag{5.2} $$

give $w_a = 250$ N/m and $w_b = 83.333$ N/m. The program returned exactly these values, and integrating the trapezoid reproduced the weight and the centre of 4.5 m. Items with their centre at $L/3$ and $2L/3$ gave triangles with one end intensity zero to $10^{-12}$. In a more demanding test, 31 items (30 random trapezoids and one point load) were superimposed on a common mesh of 212 nodes. The total load and its first moment computed from the mesh agreed with the sums over the items to relative errors of $1.3 \times 10^{-16}$ and $1.5 \times 10^{-16}$. This confirms that integration on a mesh containing every item boundary is exact.

### 5.3.2 Weight curve and weight items given by their moment

**Curve segments.** A curve through $(0, 5)$, $(10, 20)$, $(25, 0)$ and $(40, 12)$ kN/m was converted into three items. The closed forms of the weight and the first moment of a piecewise-linear curve are

$$ W = \sum_k \frac{(q_k + q_{k+1})\,L_k}{2}, \qquad M = \sum_k \frac{L_k\left[q_k\,(2x_k + x_{k+1}) + q_{k+1}\,(x_k + 2x_{k+1})\right]}{6} \tag{5.3} $$

which give $W = 365$ kN and $M = 6150$ kN·m. The three items had the weights 125, 150 and 90 kN with the centres 6.0, 15.0 and 35.0 m expected from Equation (4.1); the second and third segments are triangles with their centres at one third and two thirds of their length. The load model built from the items reproduced $W$ and $M$ with a relative error of zero, and its intensity at six probe points — including points outside the curve, where it must vanish — equalled the curve with an absolute error of zero. A curve with a jump from 10 to 30 kN/m at 20 m gave two items (the zero-length segment is skipped), the exact weight of 500 and the exact centre of 19.0 m, and the intensity on the two elements adjacent to the jump was 10 and 30.

**Equivalence with weight items.** The sum of the eight distributed items of the example vessel is itself a piecewise-linear curve with jumps. A test builds this curve (26 points) and analyses the example with the curve and the remaining point load instead of the items. The total weight and the LCG agreed with the item input exactly, and the block reactions agreed to $0$ for the rigid hull and $1.1 \times 10^{-14}$ of the largest reaction for the FEM. Both input forms therefore lead to the same load model, as intended by the design of Section 4.6.

**Moments.** An item of 400 t over $[40, 70]$ m with a moment of 22 000 t·m was given the LCG of 55.0 m, and a point load of 100 t with 3000 t·m was placed at 30.0 m; the same items given with `lcg` and `x` produced identical items. A moment that disagrees with a given LCG, or a weight of zero (which makes the LCG undefined), was rejected with a message naming the item. For the weight-curve example, the loader produced the documented totals: 4343.0 t for the curve (12 loaded segments, LCG 55.670 m) and 4643.0 t at an LCG of 56.639 m in total.

### 5.3.3 Two blocks: the statically determinate case

With exactly two blocks the reactions follow from equilibrium alone and do not depend on any stiffness. A uniform load of 1 kN/m over a 20 m hull was supported by blocks at 2 m and 15 m. Taking moments about each block gives $R_1 = 20 \cdot 5/13 = 7.6923$ kN and $R_2 = 20 \cdot 8/13 = 12.3077$ kN. The rigid-hull solver reproduced these values with relative errors of $1.2 \times 10^{-16}$ and $1.5 \times 10^{-16}$.

### 5.3.4 Equal blocks under an eccentric load

For $N$ blocks of equal stiffness under a rigid hull, the reactions have the closed form

$$ R_i = \frac{W}{N} + \frac{W e \,(x_i - \bar{x})}{\sum_j (x_j - \bar{x})^2}, \qquad e = x_G - \bar{x} \tag{5.4} $$

where $\bar{x}$ is the block centroid. Nine equal blocks between 2 m and 18 m were loaded with 1 kN/m over 20 m and a point load of 40 kN at 14 m, giving $W = 60$ kN and $x_G = 12.667$ m. The reactions ranged from 1.333 kN to 12.000 kN. The largest relative deviation from Equation (5.4) over all nine blocks was $1.4 \times 10^{-15}$.

### 5.3.5 Simply supported beam: Euler–Bernoulli and Timoshenko deflection

The FEM solver was tested on a beam of length $L = 10$ m with uniform load $q = 1$ kN/m, supported at both ends by very stiff springs ($k = 10^{15}$ N/m) and with $EI = 10^{8}$ N·m². The midspan deflection relative to the supports has the classical value [@gere2013]

$$ v_{\text{mid}} = \frac{5 q L^4}{384\, EI} \tag{5.5} $$

which equals 1.302083 mm. The program returned this value with a relative error of $1.7 \times 10^{-13}$. When the shear stiffness $\kappa G A = 10^{6}$ N was included, the Timoshenko solution [@timoshenko1921] adds a shear deflection:

$$ v_{\text{mid}} = \frac{5 q L^4}{384\, EI} + \frac{q L^2}{8\, \kappa G A} \tag{5.6} $$

which equals 13.802083 mm. The program agreed to $1.8 \times 10^{-15}$. A separate test confirmed that switching shear deformation off restores Equation (5.5), with the same error of $1.7 \times 10^{-13}$, even when $\kappa G A$ is given.

These results are exact to round-off, although the beam was divided into twenty elements of 0.5 m. This is expected: for cubic Hermite elements with consistent load vectors, the nodal displacements of a beam under polynomial loading are exact [@cook2002]. The benchmark therefore checks the element stiffness matrix, the consistent load vector and the assembly together.

### 5.3.6 Two-span continuous beam

A beam of 20 m under 1 kN/m on three equally spaced supports (at 0, 10 and 20 m) is statically indeterminate to the first degree. Its exact reactions are $3/16$, $10/16$ and $3/16$ of the total load [@gere2013], that is 3.75 kN, 12.5 kN and 3.75 kN. The program, with springs of $10^{15}$ N/m and $EI = 10^{8}$ N·m², returned 3750.0000026 N, 12499.999995 N and 3750.0000026 N. The relative deviations of $4$–$7 \times 10^{-10}$ are not a numerical error. They are the physically correct effect of the very small, but finite, compression of the supporting springs; with truly rigid supports the deviation would vanish.

### 5.3.7 Shear force and bending moment

The integration of $V(x)$ and $M(x)$ was checked on two cases.

- **Simply supported beam** with $L = 10$ m and $q = 2$ kN/m. The maximum moment is $qL^2/8 = 25{,}000$ N·m at midspan [@gere2013]. The program returned this value at $x = 5$ m with a relative error of $1.7 \times 10^{-12}$, and with a positive sign, confirming the sagging-positive convention. The largest shear force $qL/2$ agreed to $8.4 \times 10^{-13}$.
- **Beam with overhangs.** A beam of 10 m under $q = 1$ kN/m rests on supports at $a = 2$ m from each end. Over the supports the cantilevers produce the hogging moment

$$ M_{\text{support}} = -\frac{q a^2}{2} \tag{5.7} $$

which equals $-2000$ N·m. At midspan the moment is $R(L/2 - a) - q(L/2)^2/2 = 2500$ N·m. The program reproduced the support moment exactly (relative error zero) and the midspan moment with a relative error of $1.1 \times 10^{-12}$.

In both cases the shear force and bending moment returned to zero at the free end, to within $3.3 \times 10^{-8}$ N and $8.4 \times 10^{-8}$ N·m, about $10^{-12}$ of the largest values. A third test checked the stresses and utilisations of all hull girder checks of a section with given moduli and envelopes against hand calculation, including the sign of the deck stress (compression in sagging) and the failed verdict when a limit is exceeded.

### 5.3.8 Block crushing

Two identical blocks with stiffness $k = 10^{9}$ N/m, crush force $F_y = 3 \times 10^{5}$ N and post-crush stiffness ratio $0.1$ carry a symmetric load of $W = 10^{6}$ N. Each block carries $W/2 = 5 \times 10^{5}$ N, which exceeds $F_y$. The bilinear block law of Section 3.4.2 then gives the compression

$$ \delta = \frac{F_y}{k} + \frac{W/2 - F_y}{k_2}, \qquad k_2 = 0.1\,k \tag{5.8} $$

which equals 2.3 mm. The program identified both blocks as crushed and returned exactly 2.3 mm and $5 \times 10^{5}$ N per block. This test exercises the crushed branch of the block law, the constant term $c_i$ of the linearisation and the state iteration.

### 5.3.9 Gaps and lift-off

For gaps no simple closed form exists, but an exact equivalence can be tested. A block whose gap is too large to close carries no load, so the solution must equal that of the same case with the block removed. In the test, a gap of 50 mm was applied to the middle block of five. The program returned zero load on that block, and the other four reactions agreed with the four-block solution to round-off (the test tolerance is $10^{-9}$). The same equivalence was tested for a block that lifts off under an eccentric load in the rigid-hull method. A further test confirmed that a gap of 0.2 mm, which closes partially, gives a reaction $0 < R = 375.2$ kN $< R_{\text{neighbour}} = 398.5$ kN that is consistent with $R = k\,\delta$ to round-off.

### 5.3.10 Summary of the benchmarks

| Benchmark | Quantity | Exact value | Program | Relative error |
|---|---|---|---|---|
| Trapezoid integration (31 items) | total load / first moment | sums over items | — | $1.3 \times 10^{-16}$ / $1.5 \times 10^{-16}$ |
| Weight curve (Eq. 5.3) | $W$, $M$, intensity | 365 kN, 6150 kN·m, $q(x)$ | same | 0 |
| Weight curve with a jump | $W$, LCG | 500, 19.0 m | same | 0 |
| Items as one curve (example vessel) | reactions | item input | same | 0 (rigid), $1.1 \times 10^{-14}$ (FEM) |
| Moment input | LCG $= M/W$ | 55.0 m, 30.0 m | same | round-off |
| Two blocks | $R_1$, $R_2$ | 7692.31 N, 12307.69 N | 7692.31 N, 12307.69 N | $1.2 \times 10^{-16}$, $1.5 \times 10^{-16}$ |
| Eccentric, 9 equal blocks | all $R_i$ (Eq. 5.4) | 1333.3 … 12000 N | 1333.3 … 12000 N | $\le 1.4 \times 10^{-15}$ |
| Simply supported, Euler–Bernoulli | $v_{\text{mid}}$ (Eq. 5.5) | 1.302083 mm | 1.302083 mm | $1.7 \times 10^{-13}$ |
| Simply supported, Timoshenko | $v_{\text{mid}}$ (Eq. 5.6) | 13.802083 mm | 13.802083 mm | $1.8 \times 10^{-15}$ |
| Two-span beam | $R$ (3/16, 10/16, 3/16) | 3750, 12500, 3750 N | 3750.0000026, 12499.999995, 3750.0000026 N | $4$–$7 \times 10^{-10}$ (finite spring stiffness) |
| Simply supported moment | $M_{\max} = qL^2/8$ | 25 000 N·m | 25 000 N·m | $1.7 \times 10^{-12}$ |
| Overhang moment | $M$ at support (Eq. 5.7) / midspan | −2000 / 2500 N·m | −2000 / 2500 N·m | 0 / $1.1 \times 10^{-12}$ |
| Crushing | $\delta$ (Eq. 5.8) | 2.3 mm | 2.3 mm | 0 |
| Unclosable gap | reactions vs block removed | equal | equal | round-off |

Table: Table 5.2 – Analytical benchmark cases: exact values, program results and relative errors.

Every benchmark was reproduced to the limit of double-precision arithmetic, apart from the two-span case, whose small deviation is explained by the finite support stiffness used in the model, not by any error in the program.

## 5.4 Convergence studies

### 5.4.1 Convergence to the rigid-hull method

The rigid-hull method (Method A) is the limit of the elastic-hull method (Method B) as the hull bending stiffness tends to infinity. Comparing the two methods for increasing $EI$ therefore tests both solvers against each other, under loading that neither closed form above covers. The test case was a 100 m hull on 38 equal blocks of $k = 2 \times 10^{9}$ N/m at 2.5 m spacing. It carried a uniform load of 40 kN/m, a 2 MN item over 20–30 m and a 0.5 MN point load at 80 m, giving a total weight of 6.5 MN. The difference between the methods, listed in Table 5.3, was measured as

$$ \Delta = \frac{\max_i \left| R_i^{\text{FEM}} - R_i^{\text{rigid}} \right|}{\max_i R_i^{\text{rigid}}} \tag{5.9} $$

| $EI$ [N·m²] | $\Delta$ |
|---|---|
| $10^{13}$ | $3.548 \times 10^{-1}$ |
| $10^{14}$ | $8.765 \times 10^{-2}$ |
| $10^{15}$ | $1.397 \times 10^{-2}$ |
| $10^{16}$ | $1.506 \times 10^{-3}$ |
| $10^{17}$ | $1.517 \times 10^{-4}$ |
| $10^{18}$ | $1.569 \times 10^{-5}$ |
| $10^{19}$ | $3.569 \times 10^{-5}$ |
| $10^{20}$ | $3.095 \times 10^{-4}$ |

Table: Table 5.3 – Difference between FEM and rigid-hull reactions as a function of the hull bending stiffness.

From $EI = 10^{15}$ to $10^{18}$ N·m², the difference decreases by a factor of very nearly ten for each tenfold increase in $EI$. This is the first-order behaviour $\Delta \propto 1/EI$ expected from a perturbation of the rigid solution, and it is strong evidence that both solvers are correct. Above $10^{18}$ N·m², the difference rises again. This is not a discretisation effect but round-off error: the ratio between the element bending stiffness, of order $12EI/L_e^3$, and the spring stiffness becomes so large that the springs' contribution to the system matrix is lost in the last digits of double-precision arithmetic. The realistic range of hull stiffness, around $10^{12}$–$10^{13}$ N·m² for the example vessel, lies far below this round-off floor. An automated test checks the factor-of-ten reduction between $10^{16}$, $10^{17}$ and $10^{18}$ N·m², and requires $\Delta < 5 \times 10^{-5}$ at $10^{18}$ N·m².

### 5.4.2 Mesh independence

For the Euler–Bernoulli beam with consistent loads, the finite element solution is nodally exact for piecewise-polynomial loads, as noted in Section 5.3.5. Because every block and every load discontinuity is a node, the block reactions should not depend on the element size at all. This was tested on a 60 m hull on 19 blocks at 3 m spacing ($k = 10^{9}$ N/m, $EI = 2 \times 10^{11}$ N·m²), loaded with 20 kN/m and a 1 MN item over 10–20 m. Table 5.4 compares the reactions for maximum element lengths from 3.0 m to 0.25 m with a reference solution with 0.1 m elements and 601 nodes.

| Maximum element length [m] | Nodes | Max. relative difference from reference |
|---|---:|---|
| 3.0 | 23 | $7.873 \times 10^{-9}$ |
| 1.5 | 43 | $7.873 \times 10^{-9}$ |
| 0.5 | 121 | $7.873 \times 10^{-9}$ |
| 0.25 | 241 | $7.885 \times 10^{-9}$ |

Table: Table 5.4 – Block reactions for different mesh densities compared with a 0.1 m reference mesh.

All meshes, including the coarsest, in which each element spans the full block spacing, gave the same reactions to within $8 \times 10^{-9}$. The difference from the reference is practically identical for all meshes, so it originates in the reference solution: the finest mesh, with the most unknowns, accumulates the most round-off. The mesh therefore affects only the resolution of the displacement and moment curves between nodes, not the computed reactions. This property matters for the spacing design: each candidate layout has its own mesh (Section 4.8.1), and a mesh-dependent reaction would add numerical noise to the comparison of layouts.

### 5.4.3 Contact solution for large gaps

For large, uneven gaps there is no closed form, but every converged contact solution must satisfy the complementarity conditions of Section 3.7: blocks reported as open have non-positive compression and zero reaction, loaded blocks have positive compression, and the reactions balance the weight and its moment. Two cases were examined (Table 5.5).

- **Rigid hull, 40 mm parabolic end lift.** The 100 m test hull of Section 5.4.1, with the uniform load only, was given gaps of $0.04\,((x - 50)/47)^2$ m, rising to 40 mm at the end blocks. This case is part of the automated test suite.
- **Rigid hull, example vessel with joint kinks.** The example vessel was given kinks of 0.2 mrad at five erection joints (20, 40, 60, 80 and 100 m), which lifts the ends of the block line by 13.7 mm. For this case the direct state iteration did *not* converge: it cycled until the limit of 200 iterations. This was the failure that led to the gap continuation of Section 4.8.5.

| Case | Direct iteration | Iterations (total) | Open blocks | Max. $\delta$ of open blocks | Min. $\delta$ of loaded blocks | Force / moment residual |
|---|---|---:|---:|---|---|---|
| Rigid, 40 mm parabolic lift (38 blocks) | converged | 6 | 33 | $-0.272$ mm | $+0.094$ mm | $1.2 \times 10^{-16}$ / 0 |
| Rigid, example vessel, 0.2 mrad kinks (44 blocks) | failed (cycling) | 11 (continuation) | 28 | $-0.217$ mm | $+0.161$ mm | $3.4 \times 10^{-16}$ / $1.8 \times 10^{-16}$ |

Table: Table 5.5 – Contact solutions for large gaps: convergence, complementarity and equilibrium.

In both cases the complementarity conditions are satisfied with clear margins: the open blocks would need a further 0.2–0.3 mm of settlement to touch, and every loaded block is compressed. The residuals are at round-off level. The elastic-hull solution of the second case converged directly in one iteration with all 44 blocks loaded, because the flexible hull closes gaps of this size (Chapter 6). The continuation changes only the path to the solution and, because it accepts a solution only when it is consistent with the full gap vector, it has no effect on the result.

## 5.5 Verification of the welding-deformation model and the identification

### 5.5.1 Shape functions

A kink $\varphi$ at an erection joint $x_j$ must change the slope of the keel line by exactly $\varphi$:

$$ z'(x_j^{+}) - z'(x_j^{-}) = \varphi \tag{5.10} $$

With $\varphi = 2$ mrad at $x_j = 10$ m, the difference quotients on either side of the joint gave a slope change equal to $\varphi$ with a relative error of zero. For a construction block from 0 to 20 m with unit curvature, the numerical second derivative of the zone shape function was 1 inside the block and 0 outside it, with deviations of $2.5 \times 10^{-10}$ and $2.8 \times 10^{-10}$ (the accuracy of the repeated difference quotient), and the shape function at the forward end of the zone equalled its closed form $(x_b - x_a)^2/2 = 200$ m² exactly.

### 5.5.2 Invariance to a straight-line tilt of the keel

The keel lift is defined only up to a straight line, because a straight line is a rigid tilt that the support solvers absorb (Section 4.9). This was verified at two levels:

- **Gaps.** Adding $0.4 + 0.003\,x$ m to a parabolic lift changed the gaps of Equation (4.5) by at most $1.6 \times 10^{-16}$ m.
- **Reactions.** The elastic-hull solution for a 10 mm concave-up gap distribution was compared with the solution for the same gaps plus $0.01 + 2 \times 10^{-4}\,x$ m. The reactions agreed to a relative difference of $3.1 \times 10^{-7}$. This is larger than machine precision because the added rigid-body displacement (10–30 mm) is about a hundred times the elastic deformation, so its round-off is amplified by the condition number of the system matrix. The difference is nevertheless negligible for engineering purposes; the automated test uses a tolerance of $10^{-8}$ relative plus $10^{-6} W$ absolute.

### 5.5.3 End lift and lift-off of the end blocks

The phenomenon that motivates the welding model — the first and last blocks carrying no load because the hull ends have lifted — was reproduced in a controlled test. A 100 m hull of $EI = 4 \times 10^{12}$ N·m² on 38 blocks was given kinks of 1.5 mrad at four erection joints, which lifts the end blocks by 78.9 mm. The elastic-hull solution showed both end blocks open with zero reaction (30 open blocks in total). All open blocks had $\delta \le 0$, all loaded blocks $\delta > 0$, and the equilibrium residuals were $1.7 \times 10^{-9}$ and $1.2 \times 10^{-9}$. The residual is larger than usual because the largest gap term $k\,g$ is about forty times the ship weight here. For the welded-hull example, the same behaviour was confirmed on the full vessel: with welding gaps of up to 41.4 mm the first and last blocks carry no load as built, and the welding gaps are the only gaps of the case.

### 5.5.4 Identification of joint angles from a survey

**Noise-free survey.** A synthetic keel survey was generated from known joint angles (0.8, 0.4, 0.5 and 0.7 mrad at 20, 45, 70 and 95 m) at 25 stations over 120 m. An arbitrary datum offset of 0.123 m and a tilt of $-4 \times 10^{-4}$ were added. The identification recovered all four angles with a largest relative error of $3.8 \times 10^{-13}$ and an rms residual of $7.5 \times 10^{-16}$ m. This confirms that the least-squares formulation, including the straight-line terms for the survey datum, is implemented correctly.

**Noisy survey.** With noise the identification is an ill-posed inverse problem [@hansen1998], and its accuracy depends on the survey. With 0.5 mm noise at 49 stations and a regularisation weight of $10^{-6}$, the angles were recovered to within 0.048 mrad, and the rms residual of 0.525 mm matched the noise level.

**Known zone shrinkage.** When the shrinkage curvature of the construction blocks is known, it is subtracted from the survey before the joint angles are fitted (Section 3.10). A noise-free survey was generated from two joint kinks (0.6 and 0.3 mrad) and two shrinkage zones with different curvatures ($8 \times 10^{-6}$ and $3 \times 10^{-6}$ m⁻¹), plus a datum offset. The two angles were recovered with a relative error of $2.3 \times 10^{-15}$ and a residual of $1.2 \times 10^{-17}$ m. Without the subtraction the zone curvature would be absorbed into the angles; this test guards against that double counting.

**Sample survey file.** The sample survey supplied with the program (`examples/welding_survey.csv`, 25 stations, 0.5 mm noise) was read from its CSV file and the five joint angles were identified as 0.938, 0.477, 0.384, 0.497 and 0.670 mrad, within 0.070 mrad of the values 0.9, 0.5, 0.4, 0.5 and 0.6 mrad used to generate it, with an rms residual of 0.295 mm and a largest residual of 0.698 mm. The zones of the case were kept unchanged. A survey file without the required columns was rejected with an error message naming the header, and a survey used without identification was confirmed to replace the joints and zones as the total keel shape.

The surveyed keel line itself is reproduced exactly by all three interpolation methods (PCHIP, linear and natural spline): the interpolant passes through every survey point with zero error. These results show that the implementation is correct and that the identified angles are as accurate as the survey allows.

## 5.6 Verification of the spacing design

The spacing design is the final result of the program, and its verification has two parts. The algorithm of Section 3.12 is first verified in isolation with substitute analyses whose optimum is known (Sections 5.6.1–5.6.5); this is possible because the analysis is injected (Section 4.5). The integration with the full support model, the welding deformation and the GUI is then verified on the example cases (Sections 5.6.6–5.6.8). The results are summarised in Table 5.6.

### 5.6.1 Substitute analysis

The substitute analysis represents a hull from $a$ to $b$ (by default 0 to 30 in the units of the test) in which every block carries the load $q(x)$ of its tributary cell — from the midpoint to its aft neighbour to the midpoint to its forward neighbour — and the end blocks also carry the overhangs. The load is integrated numerically on 40 000 intervals and is conserved exactly. Under this model a block's reaction depends only on its own and its neighbours' positions, which makes the optimum computable by hand, while the design algorithm sees the same interface as with the FEM: a function from blocks to reactions.

### 5.6.2 Closed-form optimum for a uniform load

For a uniform load $q$ on a hull of length $L$ carried by $N$ blocks of equal strength $S$, every block can carry exactly $qL/N$ if it stands at the centre of its cell, $x_i = (i - \tfrac{1}{2})L/N$. No layout can do better, because the reactions sum to $qL$, so the smallest possible largest utilisation is (Section 3.12.7)

$$ U^\ast = \frac{qL}{N S}, \qquad U_i = U^\ast \text{ for all } i \tag{5.11} $$

The test uses $N = 11$, $L = 30$, $qL = 950$ and $S = 100$, so that $U^\ast = 0.8636$, and starts from eleven random positions with a largest utilisation of 1.305. The design reached a largest utilisation of 0.8638 after 82 analyses, a relative deviation of $1.8 \times 10^{-4}$ from $U^\ast$; the spread of the utilisations was $2.8 \times 10^{-4}\,U^\ast$, and every block was within 0.73 mm of its optimal position. Phase 1 (equidistribution) brought the largest utilisation from 1.305 to 0.950 and phase 2 (sequential linear programming) from 0.950 to 0.8638. The remaining deviation is of the order of the finite-difference step of 0.05 m and the stopping tolerances of the trust region, not a systematic error. This test verifies that the method finds the global optimum where it is known, and that the minimax formulation equalises the governing utilisations.

### 5.6.3 Weak blocks get a smaller spacing

Sixteen blocks, eight of strength 50 aft and eight of strength 150 forward, carry a uniform load of 40 over 30 (total 1200; total strength 1600). For equal utilisation the load per block must be proportional to its strength, so the spacing of the weak blocks should be one third of that of the strong blocks, and the best possible utilisation is $1200/1600 = 0.75$. The design fitted with a largest utilisation of 0.766, down from 1.547 for the evenly spaced plan. The mean spacing of the interior weak blocks was 0.952 and that of the interior strong blocks 2.862, a ratio of 0.333, and fifteen of the sixteen blocks had utilisations between 0.759 and 0.766; only the first strong block next to the weak group, at 0.629, did not govern. The strengths stayed in their order, confirming that the design changes only the positions.

### 5.6.4 Equidistribution step

A single call of `equidistribution_step` with eleven blocks at 2 m spacing on $[0, 20]$, reactions of 300 for $x < 5$ and 50 elsewhere, and equal strengths kept the end blocks at 0 and 20, kept the order, and moved the blocks towards the load peak: the first five interior blocks moved to 0.433, 0.867, 1.600, 2.467 and 3.333 m, and the forward spacing grew to 5.2 m. This verifies the inversion of the cumulative load of Section 3.12.3.

### 5.6.5 Limits, grid, fixed ends, infeasibility and invalid input

**Limits and grid.** Sixteen blocks of strength 100 under a load of 120 on the first 6 m and 15 elsewhere were designed with $s_{\min} = 0.8$, $s_{\max} = 3.0$, an extent of $-2$ to 34 and a grid of 0.4 without structure. The result fitted with a largest utilisation of 0.960 (from 2.400): every spacing lay between 0.800 and 2.400, the line from 0.4 to 24.8 within the extent, and every position on the 0.4 grid to $7 \times 10^{-15}$. The history started with the phase *start* and ended with *grid*; the rounding and the local search of phase 3 improved the largest utilisation from 0.974 after phase 2 to 0.960.

**Fixed ends.** With `fix_end_blocks`, a line of five blocks with its end blocks at 1 and 13 kept both end blocks exactly in place while the interior blocks were moved.

**Infeasibility.** Six blocks of strength 30 under a load of 600 (total strength 180) were reported as `does_not_fit`, with the best largest utilisation of 6.83, the message that the total strength is less than the load, and a suggested count of 42 blocks, $\lceil 6 \times 6.83 \rceil$.

**Invalid input.** Blocks without strength (message naming B01 and B02), six blocks at a minimum spacing of 1 m on a 4 m line ("need 5.00 m, but the block line is only 4.00 m long"), and a maximum spacing below the minimum spacing were rejected with an `InputError`.

### 5.6.6 The spacing example with the full support model

The spacing example (`examples/spacing_case.yaml`: 44 blocks with 190 t and 160 t strengths, FEM) was designed with the full support model. The given plan has a largest ratio of reaction to strength of 1.340; the designed spacing fits with 0.761 after 193 analyses. The test confirms that the designed layout has 44 blocks with the strengths in their original order, that every spacing lies between 0.8 and 4.8 m, that every block stands at a floor (multiples of 0.8 m) to $10^{-6}$ m, that all hull girder and block checks of `fem+spacing` pass, and that the force and moment residuals of the designed solution are below $10^{-9}$ (they were $-9.6 \times 10^{-15}$ and $-2.3 \times 10^{-15}$). The command-line program wrote `spacing.csv` with a header and 44 rows, `spacing_history.csv`, `spacing_design.png` and a summary with the section *Block spacing design*.

**Gaps move with their blocks.** An initial gap of 0.4 mm was given to block B31 of the plan. After the design, B31 had moved from 78.4 m to 72.0 m and still had its gap of 0.4 mm, the springs of the final evaluation contained exactly one non-zero gap — that of B31 — and the design fitted (0.745). This verifies that the design moves blocks, not gaps, as required by Section 3.4.4.

### 5.6.7 Welded hull: welding gaps recomputed at every candidate

The welded-hull example was analysed with the spacing design switched on. The test verifies the property that makes welding deformation a constraint of the design (Section 3.12.2):

- the gaps of the final evaluation equal the welding gaps of Equation (4.5) evaluated at the *designed* positions, with an absolute difference of zero, and those of the initial evaluation equal the as-built welding gaps of the given plan;
- the two gap vectors differ (by up to 9.2 mm), so the gaps did move with the blocks;
- the designed blocks, analysed again from scratch by `run_case` as an ordinary plan, give the same reactions as the design reported (relative difference zero) and the same welding gaps.

The outcome is also verified. With the full deformation (end lift 41.4 mm) the design improved the largest ratio of reaction to strength from 1.971 to 1.350 but did not fit: five blocks carry no load in the best layout, the governing block B06 stands at 20.8 m next to the erection joint at 20 m, and the message states that the gaps rather than the spacing govern and that more blocks may not help, with a suggested count of 60 blocks if the load spread. The warning appears in the analysis and the summary reports the recomputed gaps (up to 40.3 mm in the designed layout). With half the deformation (end lift 20.7 mm), the given plan is overloaded (1.547), and the design fitted with 0.917 after 104 analyses, with every hull girder and block check of the designed layout passed. These two tests together verify both branches of the diagnosis on the full model: a moderate welding deformation is absorbed by the spacing, a severe one is identified as gap-governed.

### 5.6.8 GUI round trip

The *Auto-calculate spacing* button was pressed on the spacing example in the main window, and the analysis ran in the background thread. The tests confirm that the designed positions and the unchanged strengths were written into the block plan (to $10^{-6}$ m), so that collecting the case from the editor reproduces the designed blocks; that the status line reports "Fits"; that the ratio of reaction to strength shown in the plan table equals that of the result to $10^{-3}$ and is at most 1; and that the spacing-design settings of the case are kept unchanged. Running the same dictionary through the core directly gave the same 193 analyses, identical positions and reactions identical to the command-line result (relative difference zero). A plain F5 analysis of the same case did not run the spacing design but filled the result columns of the plan. For the example vessel, which has no spacing-design settings, changed to 30 blocks at 3.2 m, the button used the floor spacing of 0.8 m as the minimum spacing; every designed spacing respected it. This case did not fit (2.005 to 1.256), and the diagnosis named the forward end block at the limit of the block line, 99.2 m, as governing and recommended a longer line — an example of the end-block branch of the diagnosis.

### 5.6.9 Summary of the spacing-design verification

| Verification case | Quantity | Expected | Program |
|---|---|---|---|
| Uniform load, 11 equal blocks (Eq. 5.11) | largest utilisation, spread | $U^\ast = 0.8636$, 0 | 0.8638 ($1.8 \times 10^{-4}$), spread $2.8 \times 10^{-4} U^\ast$ |
| 8 weak + 8 strong blocks | spacing ratio weak/strong; $U$ | 1/3; $\ge 0.75$ | 0.333; 0.766 (from 1.547) |
| Equidistribution step | ends, order, density at peak | kept, kept, higher | as expected |
| Limits and 0.4 grid | spacing, extent, grid | within 0.8–3.0, −2–34, on grid | 0.8–2.4, 0.4–24.8, $7 \times 10^{-15}$ |
| Fixed ends | end positions | 1, 13 | 1, 13 |
| Strength < load | status, suggested count | does not fit, $> N$ | does not fit, 42 |
| Spacing example (FEM, 44 blocks) | $\max R_i/S_i$, floors, checks, residuals | $\le 1$, all, pass, $< 10^{-9}$ | 0.761, all, pass, $\le 9.6 \times 10^{-15}$ |
| Initial gap of B31 | gap after the design | moves with B31 | 0.4 mm at 72.0 m |
| Welded hull, gaps at the design | gap difference; re-analysis | 0; identical | 0; 0 |
| Welded hull, half deformation | status | fits | 0.917, all checks pass |
| Welded hull, full deformation | status, diagnosis | does not fit, gap-governed | 1.350, 5 open blocks, gap-governed |
| GUI *Auto-calculate spacing* | plan, status, design | design in plan, "Fits", as CLI | as expected, identical |

Table: Table 5.6 – Verification of the spacing design.

## 5.7 Run-time self-checks

The benchmarks verify the program on simple cases. To give confidence in every individual analysis, the program also checks its own results at run time. Each support solution reports two equilibrium residuals: the force residual $(\sum R_i - W)/W$, and the moment residual $(\sum R_i x_i - W x_G)/(W \cdot \text{span})$. Each strength result reports the closure error of the shear force and bending moment at the forward end of the hull, relative to their maximum values, and a warning is issued if a closure error exceeds $10^{-6}$. Table 5.7 lists the values for the example cases of Chapter 6, including the designed spacings.

| Case and method | Force residual | Moment residual | Shear closure | Moment closure | Iterations |
|---|---|---|---|---|---:|
| Example vessel, rigid | $-1.7 \times 10^{-16}$ | 0 | $9.2 \times 10^{-16}$ | $2.3 \times 10^{-15}$ | 1 |
| Example vessel, FEM | 0 | $5.1 \times 10^{-15}$ | $2.7 \times 10^{-15}$ | $-1.6 \times 10^{-12}$ | 1 |
| Spacing example, `fem+spacing` | $-9.6 \times 10^{-15}$ | $-2.3 \times 10^{-15}$ | $-1.5 \times 10^{-13}$ | $-2.5 \times 10^{-12}$ | 1 |
| Weight-curve example, rigid | $1.6 \times 10^{-16}$ | 0 | $2.2 \times 10^{-16}$ | $1.5 \times 10^{-14}$ | 1 |
| Weight-curve example, FEM | $8.2 \times 10^{-16}$ | $-1.7 \times 10^{-16}$ | $1.6 \times 10^{-14}$ | $5.3 \times 10^{-13}$ | 1 |
| Weight-curve example, `fem+spacing` | $-1.5 \times 10^{-14}$ | $-7.2 \times 10^{-15}$ | $-2.9 \times 10^{-13}$ | $-3.7 \times 10^{-12}$ | 1 |
| Welded hull, rigid | 0 | $-9.0 \times 10^{-17}$ | $7.1 \times 10^{-16}$ | $9.5 \times 10^{-16}$ | 5 |
| Welded hull, FEM | $6.3 \times 10^{-14}$ | $9.2 \times 10^{-14}$ | $3.5 \times 10^{-13}$ | $-1.2 \times 10^{-12}$ | 3 |
| Welded hull, `fem+spacing` | $-2.2 \times 10^{-13}$ | $-2.6 \times 10^{-14}$ | $-1.0 \times 10^{-12}$ | $-8.0 \times 10^{-12}$ | 3 |
| Welded hull (half deformation), `fem+spacing` | $-3.2 \times 10^{-14}$ | $-2.0 \times 10^{-14}$ | $-2.6 \times 10^{-13}$ | $-1.0 \times 10^{-12}$ | 1 |
| Survey case, rigid | 0 | $9.0 \times 10^{-17}$ | $8.9 \times 10^{-16}$ | $1.8 \times 10^{-15}$ | 5 |
| Survey case, FEM | $-8.1 \times 10^{-15}$ | $9.3 \times 10^{-14}$ | $-4.4 \times 10^{-14}$ | $-4.3 \times 10^{-12}$ | 4 |

Table: Table 5.7 – Run-time self-checks for the example cases, including the designed spacings.

All residuals are at the level of round-off. For the rigid method, global equilibrium is satisfied exactly by construction, because the two equilibrium equations are solved directly. For the FEM, equilibrium is a consequence of solving the global system, so a small residual confirms that the system was solved accurately. The welded-hull residuals are somewhat larger than those of the example because the gap terms are large compared with the weight, but they remain below $10^{-11}$. The designed spacings are solved with the same accuracy as the given plans, which is expected because they are ordinary analyses. The self-checks do not prove that the input is correct, but they do show that the reactions printed for any case balance the weight and its moment, the first property an engineer would check by hand.

Additional run-time checks guard against invalid model states rather than numerical error. Before every solve, the program verifies that at least two separated blocks are in contact, and raises an error if not. Before the analysis it verifies that the LCG lies within the block line, and that the load model reproduces the total weight and LCG of the weight items to $10^{-9}$; a failure of the latter would be an internal error and stops the analysis.

## 5.8 Consistency of the graphical user interface and the core

Because the GUI converts its widgets into a dictionary and calls the same functions as the command-line interface (Section 4.4), its results can differ from those of the core only if this conversion loses or alters information. The tests of `test_gui.py` address this risk:

1. **Every example through the editor.** Each of the five examples was loaded into the input editor and collected again with `to_dict`. The resulting case had the same block positions and strengths, the same weight items (name, weight, LCG and extent, including the segments of the weight curve) and the same presence or absence of a spacing design as the case read directly from the file.
2. **Identical results.** The example vessel collected from the editor, the case saved from the main window and analysed with the command-line loader, and the case analysed through `AnalysisThread` with a Qt event loop all gave reactions identical to the direct analysis to a relative tolerance of $10^{-12}$. The survey case, opened with its 25 survey points from the CSV file in the table, did the same.
3. **Weights tab.** For the weight-curve example the editor showed the 18 curve points and the summary "Σ weight = 4643.0 t, LCG = 56.639 m" with "curve 4343.0 t"; the items kept their moments and no LCG was added. Replacing the LCG of the main engine of the example vessel by its moment of $380 \times 22$ t·m left the summary at 4500.0 t and 50.469 m and the item's LCG at 22.0 m; a disagreeing LCG was reported.
4. **Blocks tab.** Increasing the count of the spacing example to 46 added rows that copied the last spacing and strength; reducing it to 10, setting all spacings to 3.0 m and all strengths to 250 t, editing one spacing to 1.6 m, one gap to 0.5 mm and clearing one strength gave exactly these blocks, with the cleared strength equal to the type capacity of 200 t. The computed x column matched the case and was read-only, and an unreadable spacing was reported as "B04 spacing".
5. **Spacing design.** The settings of a loaded spacing design were written back unchanged on saving, F5 did not run the design, and *Auto-calculate spacing* wrote the design into the plan (Section 5.6.8).
6. **Sample welding data.** The example vessel, which has no welding data, was opened and the sample welding data were added through the editor, as the *Add sample welding data* action does; the welding chart was then shown, and the reactions agreed with a direct analysis of the welded-hull example. For a hull from $-5$ to 115 m the sample joints were mapped to 15, 35, 55, 75 and 95 m.
7. **Structure of the interface.** The input tabs (Ship, Weights, Hull sections, Blocks, Structure & limits, Welding deformation, Settings) and result tabs (Summary, Block forces, Plots, Block reactions, Checks, Welding data, Spacing design, Warnings) are exactly those described in Section 4.13, and every example can be opened and analysed from the menu with as many rows in the block table as blocks in the case.

Further tests concern input handling: an unreadable cell is reported with its row, a weight item with only one end of its extent is rejected, a case that refers to an unresolved survey file is reported, and an invalid input in the background thread (an LCG outside the middle third) is shown as an error message instead of crashing the application.

## 5.9 Performance and scalability

The complexity analysis of Section 3.13 predicts how the run time grows with the size of the problem. This section checks the prediction for one analysis by measurement and reports the cost of the spacing design, both in seconds and in the machine-independent number of analyses.

### 5.9.1 Measurement set-up

All measurements were made on one workstation (Intel64 Family 6 Model 158 processor, 8 logical CPUs, Windows 11, Python 3.14.6, NumPy 2.4.6, SciPy 1.18.0). The wall-clock time of each operation was measured with `time.perf_counter`, and the median of repeated runs is reported (21 runs for the contact solution, 11 for the strength post-processing and the complete analysis of the mesh study, 3 to 5 for the example cases, a single run for the welded-hull designs). The measurements are summarised by fitting a power law

$$ t = c\,x^{\alpha} \tag{5.12} $$

to the run time $t$ as a function of the problem size $x$ by least squares on a log–log scale.

### 5.9.2 Mesh refinement

The example vessel was analysed with the elastic hull only and maximum element lengths of 2.0, 1.0, 0.5, 0.25, 0.1 and 0.05 m, which gives between 96 and 2400 elements with the 44 blocks unchanged (Table 5.8, Figure 5.1). The contact solution grew from 1.3 ms to 6.3 ms while the number of elements grew 25-fold. The fitted exponent is 0.50 over the whole range and 0.67 for 504 or more elements. Both are below the linear complexity of a banded sparse solve (Section 3.13), because a fixed overhead of about 1 ms — creating the arrays, converting the sparse matrix and calling the solver from Python — dominates for small meshes. The strength post-processing grew with an exponent of 0.94 (0.99 for 504 or more elements), in close agreement with its linear complexity; it is the most expensive part of the analysis for fine meshes, because the element-wise integration of shear force and bending moment is a loop over the nodes in Python. A complete analysis took 5.4 ms with the coarsest and 68 ms with the finest mesh (exponent 0.92 for 504 or more elements). The default element length of 0.5 m, with 10 ms per complete analysis, is far from any performance limit.

| Max. element length [m] | $n_e$ | Unknowns $n$ | Contact solution [ms] | Strength [ms] | Complete analysis [ms] |
|---|---:|---:|---:|---:|---:|
| 2.0 | 96 | 194 | 1.32 | 2.84 | 5.36 |
| 1.0 | 147 | 296 | 1.34 | 4.03 | 6.81 |
| 0.5 | 255 | 512 | 1.68 | 7.63 | 10.02 |
| 0.25 | 504 | 1010 | 2.22 | 12.68 | 16.13 |
| 0.1 | 1200 | 2402 | 3.85 | 29.26 | 34.56 |
| 0.05 | 2400 | 4802 | 6.32 | 60.03 | 67.54 |

Table: Table 5.8 – Run times for mesh refinement of the example vessel (elastic hull only, one contact iteration in every case).

![Figure 5.1: Measured run times of the contact solution, the strength post-processing and the complete analysis as functions of the number of elements for the example vessel (log–log scale). The dotted line shows slope 1 for reference.](figures/fig_scaling.png)

### 5.9.3 Run time and analysis count of the spacing design

The cost of the spacing design is the number of analyses times the cost of one analysis (Section 4.10.6). The number of analyses depends only on the algorithm and the case, not on the computer, and is therefore the primary measure; Table 5.9 gives it per phase for the four designs of Chapter 6, together with the measured run time of `run_case`.

| Case (44 blocks, FEM) | Phase 1 | Phase 2 (accepted steps, Jacobians) | Phase 3 | Total analyses | $\max R_i/S_i$ | `run_case` [s] | per analysis [ms] |
|---|---:|---:|---:|---:|---:|---:|---:|
| Spacing example | 17 | 171 (4, 5) | 5 | 193 | 1.340 → 0.761 | 0.76 | 3.9 |
| Weight-curve example | 17 | 275 (6, 7) | 209 | 501 | 1.064 → 0.568 | 1.97 | 3.9 |
| Welded hull, half deformation | 17 | 81 (1, 2) | 6 | 104 | 1.547 → 0.917 | 0.52 | 5.0 |
| Welded hull, full deformation | 17 | 53 (1, 2) | 8 | 78 | 1.971 → 1.350 | 0.63 | 7.8 |

Table: Table 5.9 – Analyses per phase, result and run time of the spacing design (phase 1 includes the analyses of the given plan and of the starting layout; phase 3 includes the rounding).

Phase 1 always costs $2 + 15 = 17$ analyses: the given plan, the starting layout and fifteen equidistribution steps. It brings most of the improvement — for the spacing example from 1.340 to 0.788 — at a small fraction of the cost. Phase 2 is the most expensive phase, because each Jacobian costs up to $N$ analyses; the measured Jacobians cost between 22 and 40 analyses, fewer than 44, because blocks that are held between two neighbours at the minimum spacing cannot move in either direction and need no difference quotient. The weight-curve example needed six accepted steps and the longest local search (209 analyses, within the budget of $25N = 1100$). On the welded hull, an analysis costs more (5.0 and 7.8 ms instead of 3.9 ms) because the contact iteration needs up to three iterations when blocks are open, and phase 2 ended after one accepted step; with open blocks the utilisation is only piecewise smooth in the positions, which limits the range in which the linearisation predicts improvements.

All designs completed within two seconds, in line with requirement N6; the analysis of the given plan alone, with both methods, took 0.016–0.021 s for these cases. Because the number of analyses per phase is bounded — 17 in phase 1, at most $I_S(N + 1)$ in phase 2 with the iteration limit $I_S = 20$, and at most $1 + 25N$ in phase 3 — the cost of the design is bounded by $O(N^2)$ analyses for any case, each of cost $O(n)$, which agrees with the complexity analysis of Section 3.13. In the examples the actual counts were far below these bounds.

### 5.9.4 Limits of scale

For the problem sizes of dry docking — a few tens to about a hundred keel blocks and a few hundred elements — one analysis takes milliseconds and the spacing design a few seconds. The analysis scales linearly with the mesh size and could handle meshes far finer than necessary. The design scales with the number of analyses, which grows at most quadratically with the number of blocks through the finite-difference Jacobian. Should larger problems arise, the analyses of one Jacobian and of one local-search neighbourhood are independent and could be evaluated in parallel, and the Jacobian could be computed from one factorisation of the system matrix instead of by finite differences while the contact states do not change; Section 7.5 lists these as future work.

## 5.10 Limitations of the verification

The verification described in this chapter has the following limitations:

- **No validation.** All comparisons were made with analytical solutions of the same mathematical model, between the program's methods, or between the program and a substitute model. Whether the model predicts real block loads could not be assessed, given the idealisation of the hull as a beam, the block law, the kinematic welding model and the neglect of transverse effects. No measured block loads or keel surveys of a real hull were available. Validation against measurements, for example from load cells under instrumented blocks combined with a laser survey of the keel, would be the most valuable continuation of this work.
- **Simple benchmark geometries.** The benchmarks use uniform or simply varying loads and a few supports. The example cases of Chapter 6 are far more complex; for them, the program can only be checked through the self-checks of Section 5.7, the re-analysis of the designed layouts and engineering plausibility.
- **Contact iteration.** Convergence of the state iteration, with the gap continuation as a fallback, was demonstrated for all cases tested, including lift-off, crushing, welding gaps of several centimetres and every candidate layout of the spacing designs. No formal proof of termination was given for arbitrary combinations of gaps and crushing. The program therefore enforces an iteration limit and a minimum continuation step, and reports non-convergence as an error.
- **Optimality of the spacing design.** The minimax problem is non-convex and, with the floors, discrete. Sequential linear programming finds a local optimum, and the local search on the floors a local optimum with respect to its neighbourhood; neither is guaranteed to be global. The closed-form cases confirm that the method reaches the global optimum where it is known, but not in general. A status *does not fit* is therefore a statement about the best layout found, not a proof that no layout fits; the diagnosis of the reason (total strength, end block at the limit, gaps) makes this statement checkable by the engineer. Every designed layout, in contrast, is verified: it is the result of a full analysis and is checked like any given plan.
- **Design criterion.** The design minimises the largest ratio of reaction to block strength. The hull girder checks and the block bearing pressure are applied to the designed layout afterwards and are not constraints of the design; in the examples they were all satisfied (Chapter 6).
- **Synthetic survey data.** The identification was verified with synthetic surveys generated from the model itself. Real surveys contain systematic errors, such as datum and temperature effects, and deformation that the kink model does not represent. Their effect could not be assessed.
- **Graphical user interface.** The GUI tests cover the data path between the widgets and the core and the main interactive sequences, but not the visual layout. This was checked by inspection of screenshots.
- **Input data.** No program can detect input that is valid but wrong, such as a mistyped weight or a block strength from the wrong certificate. The cross-check of the summed weights against the declared weight and LCG reduces, but does not remove, this risk.

## 5.11 Summary

The program was verified by 131 automated tests.

- **Load, support and strength models.** The analytical benchmarks covered trapezoidal load integration, the weight curve with jumps, weight items given by their moment, determinate and indeterminate block arrangements, Euler–Bernoulli and Timoshenko beam deflection, continuous beams, shear and moment integration, block crushing and gap equivalence. They were reproduced to relative errors between zero and about $10^{-12}$; the one exception, at about $10^{-9}$, is explained physically by the finite support stiffness. The example's weight items expressed as one weight curve gave the same reactions to $10^{-14}$. The elastic-hull solution converged to the rigid-hull solution at the theoretically expected rate of $1/EI$, and the block reactions were independent of the mesh density.
- **Contact for large gaps.** Solutions with welding gaps of several centimetres satisfied the complementarity conditions and equilibrium to round-off, including a case in which only the gap continuation led to convergence.
- **Welding deformation and identification.** The kink and shrinkage shape functions matched their closed forms, the gaps and reactions were invariant to a straight-line tilt of the keel, and the end blocks were correctly unloaded by an end lift. Joint angles were recovered exactly from a noise-free survey, also with known zone shrinkage subtracted, and to within the accuracy the noise allows from a noisy one (0.070 mrad for the sample survey).
- **Spacing design.** With substitute analyses, the design reached the closed-form optimum for a uniform load to $2 \times 10^{-4}$, gave weak blocks one third of the spacing of three times stronger blocks, and respected the spacing limits, the extent, the grid and fixed ends; an insufficient total strength was diagnosed with a suggested number of blocks. With the full model, the spacing example was brought from 1.340 to 0.761 with every block under a floor and all checks passed; initial gaps moved with their blocks; on the welded hull the welding gaps were recomputed exactly at the designed positions, and a re-analysis from scratch reproduced the design; a moderate welding deformation was absorbed by the spacing and a severe one was diagnosed as gap-governed. The GUI produced the same design as the command-line program and wrote it into the block plan.
- **Run-time checks and interfaces.** Every analysis reports its equilibrium residuals and closure errors, which were of order $10^{-17}$ to $10^{-11}$ for all example cases, including the designed spacings. The GUI produced results identical to the command-line program for all examples.
- **Performance.** One analysis scaled linearly with the mesh size and took about 10 ms with the default mesh. The spacing designs of the examples needed between 78 and 501 analyses and took between 0.5 and 2 s.

The implementation is therefore considered verified with respect to the mathematical model of Chapter 3. Validation of that model against measured block loads and keel surveys remains to be done.
