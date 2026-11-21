# 6 Results and Discussion

This chapter applies the methods of Chapter 3 and the software of Chapter 4 to an illustrative docking case and its variants. Section 6.1 describes the vessel and the case files. Sections 6.2 to 6.5 report the analysis of a given block plan: the weight distribution, built from weight items and from a weight-distribution curve; the block reactions of the rigid-hull method (Method A) and of the elastic finite element model (Method B); the hull girder loads; and the block checks against the strength of every block. Section 6.6 quantifies, by five parametric studies, how the reactions depend on the block spacing and number, the extent of the block line, the hull and block stiffnesses and initial gaps. Section 6.7 treats the welded hull as built and the identification of the welding deformation from a keel survey. Section 6.8 presents the final result of the thesis, the automatic design of the block spacing from the block strengths: for the example block plan, for variants of the plan, for the weight-curve case and for the welded hull, where the welding deformation acts as a constraint. Section 6.9 discusses the results with respect to the research questions RQ1–RQ5 and states their limitations, and Section 6.10 summarises the chapter.

All numerical values in this chapter were produced with the current version of the `shipsupport` package, either by running the case files in `examples/` or by small scripts that build variants of these files (a changed block plan, stiffness, gap or design setting) and pass them to `case_from_dict` and `run_case`. Results are given in tonnes (t), tonne-metres (t·m), megapascals (MPa) and millimetres (mm). The sign conventions are those of Section 3.2: $x$ is measured from the aft perpendicular (AP), positive forward; the hull displacement $v$ is positive downward; a block reaction $R$ is positive upward on the hull; the shear force $V(x)$ is the sum of upward forces minus downward loads aft of $x$; the bending moment $M(x) = \int V\,dx$ is positive in sagging; and stresses are positive in tension. Blocks are labelled B01, B02, … from aft to forward. Run times refer to the development machine described in Section 5.9.

## 6.1 Case description

### 6.1.1 Status of the data

The vessel is a fictitious general-purpose ship of 120 m length. Its weights, section properties, block properties, strengths, permissible limits and welding parameters were chosen to be mutually consistent and of realistic magnitude, but they **are illustrative placeholders and do not describe any real ship or dock** (Section 1.4). The results therefore demonstrate the behaviour of the methods and the relative importance of the parameters; they are not design values. A real docking analysis requires the ship's docking plan and loading manual, the dock's block data and the limits of the governing standard [@milstd1625; @nstm997].

### 6.1.2 Principal particulars and block plan

Table 6.1 lists the principal particulars of the example vessel and its given block plan. The block plan is the only definition of the blocks in the program (Section 3.4.4): a block type, the position of the first block, the number of blocks, and the spacing, strength and initial gap of every block.

| Quantity | Value |
|---|---|
| Length between perpendiculars $L_{pp}$ | 120.0 m |
| Hull extent (AP to forward end) | 0.0 – 120.0 m |
| Docking weight $W$ (declared / computed) | 4500.0 t / 4500.0 t |
| Longitudinal centre of gravity (declared / computed) | 50.470 m / 50.469 m |
| Block plan (given) | 44 keel blocks from $x$ = 6.4 m, spacing 2.4 m (line 6.4 – 109.6 m), no initial gaps |
| Block strength | 200 t (type capacity) in `example_case.yaml`; individual strengths in `spacing_case.yaml` (Section 6.1.6) |
| Transverse structure (floors) | every 0.8 m from 0.0 to 120.0 m |

Table: Table 6.1 – Principal particulars and given block plan of the example vessel (illustrative data).

### 6.1.3 Weight items

The docking weight is described by nine weight items (Table 6.2). Eight are spread as trapezoids over their longitudinal extent; the rudder is a point load. Each trapezoid preserves the weight and the longitudinal centre of gravity of its item (Section 3.3.2), which is possible because every LCG lies in the middle third of its extent (Section 3.3.3).

| Item | Weight [t] | LCG [m] | Extent [m] |
|---|---:|---:|---|
| Hull steel | 2600 | 58.5 | 0.0 – 120.0 |
| Main engine | 380 | 22.0 | 16.0 – 28.0 |
| Auxiliary machinery | 220 | 26.0 | 14.0 – 38.0 |
| Superstructure | 420 | 14.0 | 4.0 – 24.0 |
| Outfit | 520 | 62.0 | 4.0 – 118.0 |
| Propeller and shaft | 60 | 6.0 | 0.0 – 12.0 |
| Rudder | 25 | 2.0 | point load |
| Anchors and chain | 75 | 112.0 | 104.0 – 120.0 |
| Residual fuel | 200 | 70.0 | 60.0 – 80.0 |
| **Total** | **4500** | **50.469** | |

Table: Table 6.2 – Weight items of the example vessel.

### 6.1.4 Hull girder properties

The hull girder is represented by three prismatic regions (Table 6.3). The bending stiffness $EI$ follows from the second moment of area with $E = 206$ GPa, and the shear stiffness $\kappa G A$ from the effective shear area with $G = 79.3$ GPa. The midship region is 2.5 times stiffer in bending than the aft region and 3.3 times stiffer than the forward region.

| Region [m] | $I$ [m⁴] | $EI$ [N·m²] | $A_s$ [m²] | $\kappa GA$ [N] | $Z_{deck}$ [m³] | $Z_{keel}$ [m³] |
|---|---:|---:|---:|---:|---:|---:|
| 0 – 20 | 8.0 | 1.648 × 10¹² | 0.5 | 3.965 × 10¹⁰ | 1.4 | 1.8 |
| 20 – 100 | 20.0 | 4.120 × 10¹² | 0.9 | 7.137 × 10¹⁰ | 3.6 | 4.4 |
| 100 – 120 | 6.0 | 1.236 × 10¹² | 0.4 | 3.172 × 10¹⁰ | 1.1 | 1.4 |

Table: Table 6.3 – Hull girder section properties.

The permissible still-water bending moment rises from 3000 t·m at the ends to 20 000 t·m at 30 m and 90 m and 25 000 t·m at 60 m (equal for sagging and hogging). The permissible shear force rises from 600 t at the ends to 1500 t between 30 m and 90 m. The allowable hull stresses are 175 MPa in bending and 110 MPa in shear.

### 6.1.5 Block type

All blocks are of one type (Table 6.4). The stiffness follows from three layers acting as springs in series over the contact area (Section 3.4.1); the softwood cap contributes 52 % of the compliance. The block follows the bilinear law of Section 3.4.2: linear up to the crush force, then 5 % of the initial stiffness.

| Property | Value |
|---|---|
| Layers (thickness, modulus) | concrete 0.90 m, 30 000 MPa; hardwood 0.20 m, 1 000 MPa; softwood cap 0.10 m, 400 MPa |
| Contact area | 0.8 m × 1.2 m = 0.96 m² |
| Initial stiffness $k$ | 2000 kN/mm |
| Capacity (default strength of every block) | 200 t |
| Allowable bearing pressure | 2.0 MPa (196 t over the contact area) |
| Crush stress / crush force | 4.0 MPa / 391.6 t |
| Compression at the onset of crushing | 1.92 mm |
| Post-crush stiffness ratio | 0.05 |

Table: Table 6.4 – Keel block type.

### 6.1.6 Case files

Table 6.5 lists the case files used in this chapter. They share the hull, the block type and the floors; they differ in the weight description, the block strengths, the welding deformation and whether the spacing design is enabled.

| Case file | Weights | Block plan and strengths | Welding | Spacing design | Sections |
|---|---|---|---|---|---|
| `example_case.yaml` | 9 items (LCG) | 44 × 200 t, 2.4 m | – | off | 6.2–6.6 |
| `weight_curve_case.yaml` | curve (18 points) + 3 items (moment) | 44 × 200 t, 2.4 m | – | on | 6.2, 6.8.3 |
| `spacing_case.yaml` | 9 items (LCG) | 44 blocks, 2.4 m; 160 t for B04–B12, 190 t otherwise (Σ S = 8090 t) | – | on | 6.5, 6.8.1–6.8.2 |
| `welding_case.yaml` | 9 items (LCG) | 44 × 200 t, 2.4 m | prescribed joints and zones | off (switched on in Section 6.8.4) | 6.7, 6.8.4 |
| `welding_survey_case.yaml` | 9 items (LCG) | 44 × 200 t, 2.4 m | survey `welding_survey.csv`, joint angles identified | off | 6.7.4 |

Table: Table 6.5 – Case files of the example vessel.

The block plan of `spacing_case.yaml` represents a realistic situation in a dock: the blocks available are not all equal, and the nine older blocks B04–B12, which in the given plan lie under the engine room, have a lower strength of 160 t. Both strengths are below the 196 t at which the bearing pressure reaches its allowable value, so that the block strength, and not the bearing pressure, is the binding block limit. The settings of the spacing design are those of Section 4.7.4: target utilisation 1.0, minimum spacing $s_{\min}$ = 0.8 m (the floor spacing), maximum spacing $s_{\max}$ = 4.8 m (twice the given spacing), block line within 6.4–109.6 m, final positions at the floors, at most 20 sequential-LP iterations.

## 6.2 Weight distribution

### 6.2.1 Weight items

Figure 6.1 shows the longitudinal weight distribution assembled from the nine items of Table 6.2. It is strongly non-uniform. Its maximum of 89.06 t/m lies just forward of $x$ = 16 m, where the main engine, auxiliary machinery, superstructure, outfit and hull steel overlap; forward of about 38 m the intensity is roughly 25–36 t/m. The mesh of the example has 256 nodes. Its node set includes every item boundary, section boundary and block position (Section 3.3.5), so every element carries an exactly linear load and the integrals of Section 3.3.6 are exact.

![Figure 6.1: Longitudinal weight distribution of the example vessel, built from the weight items of Table 6.2; the rudder is a point load at x = 2.0 m.](figures/fig_weight_distribution.png)

The distributed load integrates to $W$ = 4500.0 t with LCG = 50.469 m, which reproduces the sum of the items to machine precision, as the run-time weight check of the program requires (Section 5.7). The declared LCG of 50.470 m differs from the computed value by 0.001 m, well within the tolerance, so no weight warning is raised.

Two features of the distribution decide the block reactions. First, the LCG at 50.47 m lies aft of the centre of the block line at 58.0 m, so the aft blocks must carry more than the forward blocks. Second, a considerable part of the weight lies outside the block line: the hull aft of B01 ($x$ = 6.4 m) carries 266.3 t, and the hull forward of B44 ($x$ = 109.6 m) carries 298.8 t. These overhangs are cantilevered from the end blocks, and Section 6.3 shows that they control the end-block reactions.

### 6.2.2 Weight-distribution curve and items given by moment

In the weight-curve case the lightship is given as a piecewise-linear weight-distribution curve of 18 points with jumps at $x$ = 4, 16, 28, 38 and 104 m, and three deadweight items are given by their longitudinal moment about the AP instead of their LCG (Section 3.3.1). Each of the 12 segments of non-zero length of the curve is converted into an exact trapezoidal item (Section 3.3.4), so that the program treats curve and items in the same way. Table 6.6 lists the totals and Figure 6.2 shows the curve, the items and the total distribution.

| Component | Weight [t] | Moment about AP [t·m] | LCG [m] |
|---|---:|---:|---:|
| Lightship curve (18 points, 12 trapezoids) | 4343.0 | 241 775 | 55.670 |
| Residual fuel (60–80 m) | 200.0 | 14 000 | 70.000 |
| Fresh water (94–106 m) | 60.0 | 6 000 | 100.000 |
| Stores (point load) | 40.0 | 1 200 | 30.000 |
| Deadweight items together | 300.0 | 21 200 | 70.667 |
| **Total** | **4643.0** | **262 975** | **56.639** |

Table: Table 6.6 – Weight-curve case: lightship curve and deadweight items given by their moments. The moment of the curve is $W \cdot$ LCG from the program; the item LCGs follow from LCG = moment / weight.

![Figure 6.2: Weight distribution of the weight-curve case: lightship curve with its points and jumps, deadweight items given by their moments (dashed; the stores as a point load) and the total distribution w(x).](figures/fig_weight_curve.png)

The computed total of 4643.0 t at LCG 56.639 m agrees with the declared values, and the curve alone integrates to 4343.0 t at 55.670 m, as the label of Figure 6.2 shows. The largest intensity is the 64 t/m of the engine room (16–28 m); the jumps appear as vertical steps in the total. Compared with the item description of the example, this distribution is smoother and its LCG lies 6.2 m further forward, close to the centre of the block line. The overhangs carry 121.7 t aft of B01 and 243.3 t forward of B44.

## 6.3 Block reactions: rigid versus elastic hull

### 6.3.1 Overview

Figure 6.3 shows the force on each of the 44 blocks of `example_case.yaml` predicted by the elastic model (bars), with the rigid-hull results as markers, and Figure 6.4 shows the weight distribution and the reactions together with the settlement, shear force and bending moment of both methods. Both solutions satisfy global equilibrium: the reactions sum to 4500.0 t, and the relative force and moment residuals of the elastic solution were below $10^{-14}$. The contact iteration converged in one iteration in both methods, because no block lifted off or crushed.

![Figure 6.3: Weight distribution (top) and force on every block of the example vessel. Bars: elastic hull (Method B), coloured by utilisation of the block limit; diamonds: rigid hull (Method A); dashed line: block capacity of 200 t; dotted line: average load W/N = 102.3 t.](figures/fig_block_forces_example.png)

![Figure 6.4: Example vessel. From top: weight distribution and block reactions; hull settlement; shear force; bending moment. Both methods are shown with the permissible envelopes.](figures/fig_results_example.png)

### 6.3.2 Selected block reactions

Table 6.7 compares the two methods for the five aftmost blocks, a midship block and the five foremost blocks.

| Block | $x$ [m] | Rigid $R$ [t] | Elastic $R$ [t] | Elastic / rigid | Elastic $\delta$ [mm] | Elastic $p$ [MPa] |
|---|---:|---:|---:|---:|---:|---:|
| B01 | 6.4 | 145.06 | 254.63 | 1.76 | 1.249 | 2.601 |
| B02 | 8.8 | 143.07 | 211.75 | 1.48 | 1.038 | 2.163 |
| B03 | 11.2 | 141.08 | 185.17 | 1.31 | 0.908 | 1.892 |
| B04 | 13.6 | 139.09 | 171.95 | 1.24 | 0.843 | 1.756 |
| B05 | 16.0 | 137.10 | 169.52 | 1.24 | 0.831 | 1.732 |
| B21 | 54.4 | 105.26 | 66.61 | 0.63 | 0.327 | 0.680 |
| B40 | 100.0 | 67.44 | 70.28 | 1.04 | 0.345 | 0.718 |
| B41 | 102.4 | 65.45 | 82.34 | 1.26 | 0.404 | 0.841 |
| B42 | 104.8 | 63.46 | 112.31 | 1.77 | 0.551 | 1.147 |
| B43 | 107.2 | 61.47 | 164.47 | 2.68 | 0.806 | 1.680 |
| B44 | 109.6 | 59.48 | 246.26 | 4.14 | 1.207 | 2.516 |

Table: Table 6.7 – Reactions, compressions and bearing pressures of selected blocks (`example_case.yaml`).

The rigid method gives a linear distribution, from 145.06 t at B01 to 59.48 t at B44. The elastic model gives a markedly different one. The largest reaction is 254.63 t at B01, 76 % above the rigid value; the reaction at B44 is 246.26 t, 4.1 times the rigid value. Between about 30 m and 100 m the elastic reactions are lower than the rigid ones and follow the local weight intensity; the smallest is 54.39 t at B36 ($x$ = 90.4 m).

### 6.3.3 Physical explanation

The difference has two causes.

**Overhang cantilevers.** In the elastic model the 266.3 t of hull aft of B01 and the 298.8 t forward of B44 act as cantilevers. Their weight enters the hull girder near the end blocks and, because the hull can bend, it is carried mainly by the end blocks and their neighbours. In the rigid model the same weight enters only through the global force and moment balance and is spread over all 44 blocks by the planar displacement $v = a + bx$, which cannot represent a local concentration of load.

**Local load following.** A flexible hull settles more where it is heavily loaded. The elastic reactions therefore track the weight distribution: about 166–172 t under the machinery space (13.6–20.8 m) and about 65–68 t under the light midbody (45–50 m). A rigid hull forces all blocks into one plane, so its reactions are linear whatever the weight distribution.

The distance over which the hull spreads a concentrated load is governed by the characteristic parameter of a beam on an elastic foundation [@hetenyi1946] (Section 2.6):

$$ \lambda = \left(\frac{k_f}{4\,EI}\right)^{1/4}, \qquad k_f = \frac{k}{s} \tag{6.1} $$

where $k_f$ is the foundation modulus obtained by smearing the block stiffness $k$ over the spacing $s$. For the midship region ($EI = 4.12 \times 10^{12}$ N·m², $k = 2.0 \times 10^9$ N/m, $s$ = 2.4 m), Equation (6.1) gives $\lambda = 0.0843$ m⁻¹, a characteristic length $1/\lambda$ = 11.86 m and a half-wavelength $\pi/\lambda$ = 37.3 m. A load is therefore spread over about 15 blocks on each side, much less than the 103.2 m of the block line, and the hull behaves as a long beam on an elastic foundation rather than as a rigid body. The end regions are more flexible still, since their $EI$ is 2.5–3.3 times smaller.

The settlement curves of Figure 6.4 confirm this. The rigid solution is a straight line from 0.74 mm at the stern to 0.25 mm at the bow. The elastic hull settles by 1.96 mm at the aft tip and 3.60 mm at the forward tip: both overhangs sag onto the end blocks.

### 6.3.4 Weight-curve case

For the smoother weight distribution of the weight-curve case (Section 6.2.2), the rigid method gives reactions between 113.50 t (B01) and 97.54 t (B44), and every block passes. The elastic model again concentrates the load at the ends: 128.14 t at B01 and 212.80 t at B44, which exceeds the block capacity of 200 t by 6.4 % and fails the block checks (utilisation 1.064 against strength, 1.087 in bearing pressure). The forward end governs here because the forward overhang (243.3 t) is twice the aft one (121.7 t). The same qualitative difference between the two methods appears with a different weight description, which shows that it is a property of the support problem and not of a particular load case.

## 6.4 Hull girder loads and stresses

### 6.4.1 Shear force and bending moment

The shear force and bending moment are obtained by exact integration of the linear element loads, point loads and reactions (Section 3.8). In both methods the curves close to zero at the forward end within round-off, which confirms global equilibrium.

The two methods give fundamentally different bending moments (Figure 6.4, bottom). With the rigid hull, the moment is hogging over almost the whole length and reaches −12 298.9 t·m near $x$ = 54 m. With the elastic hull, the largest moment in magnitude is −1 723.4 t·m at $x$ = 107.2 m, in the forward overhang region, and the largest sagging moment is +143.1 t·m at $x$ = 69.9 m. The largest shear forces are 525.9 t (rigid) and 298.8 t (elastic).

### 6.4.2 Governing checks

Table 6.8 lists the governing hull girder checks of both methods. The governing location is where the ratio of demand to capacity is largest, which is not necessarily where the demand is largest, because the permissible envelopes vary along the length. The utilisation of a check is

$$ \eta = \frac{|D|}{C} \tag{6.2} $$

where $D$ is the demand and $C$ the capacity at the same location; for the bending moment, $C$ is the sagging limit where $M \ge 0$ and the hogging limit where $M < 0$.

| Check | Rigid: demand | Rigid: limit | Rigid: $\eta$ | Elastic: demand | Elastic: limit | Elastic: $\eta$ |
|---|---:|---:|---:|---:|---:|---:|
| Bending moment [t·m] | −12 191.8 (49.6 m) | 23 266.7 | 0.524 | −1 511.4 (109.6 m) | 8 893.3 | 0.170 |
| Shear force [t] | −525.9 (28.0 m) | 1 440.0 | 0.365 | −266.3 (6.4 m) | 792.0 | 0.336 |
| Deck stress [MPa] | 37.66 (100.06 m) | 175 | 0.215 | 15.36 (107.2 m) | 175 | 0.088 |
| Keel stress [MPa] | −29.59 (100.06 m) | 175 | 0.169 | −12.07 (107.2 m) | 175 | 0.069 |
| Shear stress [MPa] | 7.89 (102.4 m) | 110 | 0.072 | 7.33 (109.6 m) | 110 | 0.067 |

Table: Table 6.8 – Governing hull girder checks of `example_case.yaml` (location of the governing value in brackets).

All hull girder checks pass with both methods; the highest utilisation is 0.524, for the bending moment of the rigid method. The deck stresses are tensile because the hull is in hogging.

### 6.4.3 Why the rigid method overstates hogging

The rigid method overstates the largest bending moment by a factor of 7.1 (12 298.9 / 1 723.4). This follows from equilibrium. The bending moment at a section is the moment of the net load $w - R$ on one side of it. In the rigid solution the reactions are linear in $x$, while the weight is concentrated aft (machinery, superstructure) and in the overhangs; integrating this mismatch over half the length produces a large hogging moment amidships. In the elastic solution the reactions follow the weight locally (Section 6.3.3), so the net load on every segment is small and so is the moment. The remaining elastic moments lie in the overhangs, where they are statically determined by the cantilevered weight: the governing value of −1 511.4 t·m at $x$ = 109.6 m is the moment of the forward overhang about the last block.

For the hull girder, the rigid method is therefore conservative in this example. For the blocks it is not: it underestimates the most critical block load by 43 % (145.06 t instead of 254.63 t). A method that is conservative for one quantity is not necessarily conservative for another; this is the central practical finding of the comparison, and it is the reason why the spacing design of Section 6.8 uses the elastic model.

## 6.5 Block checks against block strength

The block checks compare the bearing pressure $p = R/A$ with the allowable pressure and the reaction with the strength $S_i$ of each block (Section 3.8). Table 6.9 lists, for the elastic solution of both block plans, every block whose utilisation exceeds 0.9.

| Case | Block | $x$ [m] | $R$ [t] | $S_i$ [t] | $R/S_i$ | $p$ [MPa] | $p/p_{allow}$ |
|---|---|---:|---:|---:|---:|---:|---:|
| `example_case` | B01 | 6.4 | 254.63 | 200 | 1.273 | 2.601 | 1.301 |
| | B02 | 8.8 | 211.75 | 200 | 1.059 | 2.163 | 1.082 |
| | B03 | 11.2 | 185.17 | 200 | 0.926 | 1.892 | 0.946 |
| | B44 | 109.6 | 246.26 | 200 | 1.231 | 2.516 | 1.258 |
| `spacing_case` | B01 | 6.4 | 254.63 | 190 | 1.340 | 2.601 | 1.301 |
| | B02 | 8.8 | 211.75 | 190 | 1.114 | 2.163 | 1.082 |
| | B03 | 11.2 | 185.17 | 190 | 0.975 | 1.892 | 0.946 |
| | B04 – B07 | 13.6 – 20.8 | 165.7 – 171.9 | 160 | 1.036 – 1.075 | 1.69 – 1.76 | 0.85 – 0.88 |
| | B08 | 23.2 | 155.2 | 160 | 0.970 | 1.585 | 0.793 |
| | B43 | 107.2 | 164.47 | 190 | 0.866 | 1.680 | 0.840 |
| | B44 | 109.6 | 246.26 | 190 | 1.296 | 2.516 | 1.258 |

Table: Table 6.9 – Highly utilised blocks in the elastic solution of the given block plans (2.4 m spacing). B43 of `spacing_case` is listed for comparison.

With the uniform capacity of 200 t (`example_case.yaml`), three blocks fail: B01, B02 and B44, with a largest utilisation of 1.301 (bearing pressure at B01). No block reaches the crush force of 391.6 t. With the rigid method the largest utilisations are 0.741 (pressure) and 0.725 (strength), both at B01, and every check passes. **The rigid method thus fails to identify three overloaded blocks**, and the overall verdict of the case changes from "passed" to "failed" when the hull is modelled as elastic.

With the individual strengths of `spacing_case.yaml`, seven blocks are loaded beyond their strength: the end blocks B01 and B44, B02, and four of the 160 t blocks under the engine room (B04–B07), where the elastic hull follows the heavy machinery weight. The largest ratio of reaction to strength is 1.340 at B01. The rigid method gives a largest ratio of 0.869, at B04 (139.1 t on a 160 t block), and again passes every check. The given block plan is therefore not acceptable, although the total strength of its blocks (8090 t) is 1.8 times the docking weight: the problem is not the amount of strength but its distribution along the hull. This is the situation that the spacing design of Section 6.8 resolves.

All blocks of the given plans lie at floors, because 6.4 m and the spacing of 2.4 m are multiples of the floor spacing of 0.8 m, so no off-structure warning is raised.

## 6.6 Parametric studies

The parametric studies quantify how the reactions of the given block plan depend on the main design and modelling parameters. They use `example_case.yaml` (44 blocks of 200 t, no initial gaps, both methods) and change one group of parameters at a time. Its baseline elastic solution has a largest reaction of 254.63 t at B01, 246.26 t at B44 and a smallest reaction of 54.39 t; the rigid baseline has a largest reaction of 145.06 t and a largest bending moment of 12 298.9 t·m. In all tables "failing blocks" is the number of blocks with a pressure or strength utilisation above 1.0 in the elastic solution. The studies motivate the spacing design: they show which changes of the block plan relieve the governing blocks, and by how much.

### 6.6.1 Block spacing and number of blocks

In the first study the block line starts at 6.4 m, the spacing is varied from 1.6 m to 4.8 m, and the number of blocks is the largest that keeps the last block at or before 109.6 m (Table 6.10 and Figure 6.5).

| Spacing [m] | Blocks | Last block [m] | Max. $R$ [t] | Aft end $R$ [t] | Fwd end $R$ [t] | Min. $R$ [t] | Rigid max. $R$ [t] | Failing blocks |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.6 | 65 | 108.8 | 216.54 | 195.37 | 216.54 | 33.06 | 97.28 | 1 |
| 2.0 | 52 | 108.4 | 248.70 | 226.03 | 248.70 | 42.45 | 120.49 | 2 |
| 2.4 | 44 | 109.6 | 254.63 | 254.63 | 246.26 | 54.39 | 145.06 | 3 |
| 2.8 | 37 | 107.2 | 315.41 | 281.42 | 315.41 | 61.67 | 165.04 | 7 |
| 3.2 | 33 | 108.8 | 306.90 | 306.90 | 299.64 | 74.76 | 189.99 | 7 |
| 3.6 | 29 | 107.2 | 352.17 | 330.93 | 352.17 | 83.94 | 209.75 | 8 |
| 4.0 | 26 | 106.4 | 386.32 | 354.11 | 386.32 | 93.28 | 230.11 | 8 |
| 4.8 | 22 | 107.2 | 392.23 | 392.20 | 392.23 | 118.02 | 274.93 | 8 |

Table: Table 6.10 – Effect of a uniform block spacing (block line from 6.4 m, ending at or before 109.6 m).

![Figure 6.5: Effect of a uniform block spacing on the largest and end-block reactions of the elastic model and on the largest reaction of the rigid model; N is the number of blocks.](figures/fig_spacing_study.png)

The largest reaction grows with the spacing from 216.5 t to 392.2 t. Even the closest spacing of 1.6 m (65 blocks) leaves one block above its capacity. The growth is not monotonic: 2.8 m gives a higher maximum than 3.2 m, because the last block then lies at 107.2 m instead of 108.8 m and the forward overhang is 1.6 m longer. At 4.8 m both end blocks crush. The ratio of the largest elastic to the largest rigid reaction falls from 2.23 at 1.6 m to 1.43 at 4.8 m, but the absolute discrepancy stays large. Spacings of 2.0, 2.8 and 3.6 m are not multiples of the floor spacing and place 26, 18 and 14 blocks off the floors.

Table 6.11 prescribes the number of blocks on the fixed line 6.4–109.6 m. Reducing the count from 60 to 20 nearly doubles the largest reaction.

| Blocks | Spacing [m] | Max. $R$ [t] | Min. $R$ [t] | Rigid max. $R$ [t] | Crushed | Failing blocks | Off floors |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 5.432 | 397.18 | 138.54 | 314.13 | 2 | 8 | 12 |
| 30 | 3.559 | 328.46 | 86.22 | 211.44 | 0 | 7 | 18 |
| 44 | 2.400 | 254.63 | 54.39 | 145.06 | 0 | 3 | 0 |
| 60 | 1.749 | 211.40 | 37.72 | 106.76 | 0 | 2 | 36 |

Table: Table 6.11 – Effect of the number of blocks on the fixed block line 6.4 – 109.6 m (uniform spacing).

Both tables show the limitation of a *uniform* spacing: more blocks reduce the largest reaction only slowly, because the end blocks always carry the overhangs, while the least loaded blocks, amidships, use only 17–27 % of their capacity when 44 or more blocks are used. A uniform plan cannot fit the reactions to the strengths; a designed, non-uniform spacing can (Section 6.8).

### 6.6.2 Extent of the block line

The second study varies the extent of the block line at a constant spacing of 2.4 m (Table 6.12), which directly changes the overhangs.

| Block line [m] | Blocks | Max. $R$ [t] (block) | Aft end $R$ [t] | Fwd end $R$ [t] | Max. \|M\| elastic [t·m] | Max. \|M\| rigid [t·m] | Failing blocks | All checks |
|---|---:|---|---:|---:|---:|---:|---:|---|
| 6.4 – 109.6 (baseline) | 44 | 254.63 (B01) | 254.63 | 246.26 | 1 723.4 | 12 298.9 | 3 | fail |
| 4.0 – 109.6 | 45 | 246.26 (B45) | 174.24 | 246.26 | 1 723.4 | 10 322.5 | 1 | fail |
| 6.4 – 114.4 | 46 | 254.63 (B01) | 254.63 | 146.05 | 943.0 | 11 112.1 | 2 | fail |
| 4.0 – 114.4 | 47 | 174.25 (B01) | 174.25 | 146.05 | 539.2 | 8 844.8 | 0 | **pass** |
| 6.4 – 100.0 | 40 | 393.49 (B40) | 254.66 | 393.49 | 6 207.5 | 16 060.4 | 5 | fail |
| 10.0 – 108.4 | 42 | 391.78 (B01) | 391.78 | 270.87 | 2 379.5 | 15 911.9 | 5 | fail |

Table: Table 6.12 – Effect of the extent of the block line (uniform spacing 2.4 m).

The extent has the strongest influence of all parameters studied. Extending the line by 2.4 m aft and 4.8 m forward (4.0–114.4 m, three more blocks) reduces the largest reaction by 32 %, from 254.63 t to 174.25 t, and the largest elastic bending moment by 69 %, to 539.2 t·m. It is the only configuration of this study in which every check passes, with largest utilisations of 0.890 (bearing pressure) and 0.871 (strength) at the new aft end block. Each extension acts locally: extending only aft relieves the aft end block but leaves the forward end unchanged, and vice versa. Shortening the line forward to 100.0 m raises the forward end-block reaction to 393.49 t, crushes that block and raises the elastic bending moment to 6 207.5 t·m because of the 20 m cantilever. Starting the line at 10.0 m has the same effect at the aft end (391.78 t) and, since 10.0 m is not a multiple of 0.8 m, places all 42 blocks off the floors.

### 6.6.3 Hull bending stiffness

In the third study the second moments of area of all regions are scaled by a common factor, with unchanged shear areas (Table 6.13 and Figure 6.6). The rigid solution does not depend on $EI$ and stays at 145.06 t (aft end) and 59.48 t (forward end).

| $EI$ factor | Max. $R$ [t] (block) | Aft end $R$ [t] | Fwd end $R$ [t] | Min. $R$ [t] | Max. \|M\| [t·m] | Failing blocks |
|---:|---|---:|---:|---:|---:|---:|
| 0.1 | 394.68 (B44) | 322.38 | 394.68 | 14.65 | 1 511.4 | 3 |
| 0.25 | 336.55 (B44) | 287.75 | 336.55 | 38.06 | 1 511.4 | 3 |
| 0.5 | 286.44 (B44) | 269.04 | 286.44 | 47.82 | 1 627.0 | 3 |
| 1 | 254.63 (B01) | 254.63 | 246.26 | 54.39 | 1 723.4 | 3 |
| 2 | 242.96 (B01) | 242.96 | 216.16 | 59.37 | 1 879.2 | 3 |
| 5 | 230.04 (B01) | 230.04 | 188.94 | 63.95 | 2 104.2 | 2 |
| 10 | 221.41 (B01) | 221.41 | 174.78 | 66.29 | 2 294.1 | 2 |
| 100 | 203.73 (B01) | 203.73 | 152.07 | 69.52 | 3 470.7 | 1 |
| 1000 | 200.34 (B01) | 200.34 | 148.24 | 69.96 | 3 829.9 | 1 |

Table: Table 6.13 – Effect of the hull bending stiffness (second moments of area of all regions scaled; shear areas unchanged).

![Figure 6.6: Effect of the hull bending stiffness on the end-block and largest reactions (left) and on the largest bending moment (right); dotted lines show the rigid-hull end-block reactions.](figures/fig_ei_study.png)

A stiffer hull spreads the overhang loads over more blocks: the largest reaction falls from 394.68 t at a tenth of the baseline stiffness to 200.34 t at a thousandfold stiffness, and the bending moment rises from 1 511.4 t·m to 3 829.9 t·m. Even at a thousandfold bending stiffness the aft end-block reaction is 38 % above the rigid value and the moment is less than a third of the rigid one. The reason is the shear flexibility, which was not scaled: the rigid method is the limit of the elastic model when both $EI$ and $\kappa GA$ tend to infinity (verified in Section 5.4 for a hull without shear deformation), and real hulls are far from that limit. For the most flexible variants the largest moment equals the statically determined moment of the forward overhang about the last block (1 511.4 t·m), because the interior moments are smaller still.

### 6.6.4 Block stiffness

In the fourth study the block stiffness is prescribed directly between 500 and 8000 kN/mm; in a second variant only the modulus of the softwood cap is changed (Table 6.14).

| $k$ [kN/mm] | Max. $R$ [t] (block) | Aft end $R$ [t] | Fwd end $R$ [t] | Min. $R$ [t] | Max. \|M\| [t·m] | Failing blocks |
|---:|---|---:|---:|---:|---:|---:|
| 500 | 229.25 (B01) | 229.25 | 176.64 | 63.84 | 2 177.0 | 2 |
| 1000 | 240.46 (B01) | 240.46 | 207.81 | 59.64 | 1 917.1 | 3 |
| 2000 | 254.63 (B01) | 254.63 | 246.26 | 54.39 | 1 723.4 | 3 |
| 4000 | 290.97 (B44) | 274.14 | 290.97 | 50.07 | 1 616.1 | 3 |
| 8000 | 339.68 (B44) | 299.72 | 339.68 | 47.79 | 1 511.4 | 3 |
| 1315 (cap 200 MPa) | 245.57 (B01) | 245.57 | 222.16 | 57.67 | 1 838.9 | 3 |
| 2704 (cap 800 MPa) | 265.05 (B44) | 262.36 | 265.05 | 52.19 | 1 678.3 | 3 |
| 3282 (cap 1600 MPa) | 277.67 (B44) | 267.93 | 277.67 | 50.98 | 1 648.0 | 3 |

Table: Table 6.14 – Effect of the block stiffness (upper part: prescribed stiffness; lower part: modulus of the softwood cap varied, stiffness from the layers).

Stiffer blocks raise the end-block reactions; softer blocks act like a stiffer hull, because both reduce $\lambda$ in Equation (6.1). The equivalence holds quantitatively: quadrupling the block stiffness (8000 kN/mm) gives a largest reaction of 339.68 t, within 1 % of the 336.55 t obtained by quartering $EI$; the small difference arises because the shear stiffness is not scaled. Doubling the cap modulus raises the block stiffness by only 35 % (from 2000 to 2704 kN/mm), because the other layers then dominate the compliance. Soft caps are a means of spreading load, but over the tested range they do not remove the overload of the end blocks.

### 6.6.5 Initial gaps

The fifth study gives a single block an initial gap through the `gap_mm` entry of the block plan (Section 3.4.4), which moves with its block (Table 6.15).

| Gapped block | Gap [mm] | $R$ of gapped block [t] | Neighbour reactions [t] | Max. $R$ [t] (block) | State |
|---|---:|---:|---|---|---|
| B21 (54.4 m) | 0.0 | 66.61 | 65.44 / 68.56 | 254.63 (B01) | elastic |
| B21 (54.4 m) | 0.3 | 14.85 | 72.94 / 76.06 | 254.50 (B01) | elastic |
| B21 (54.4 m) | 0.5, 1.0 or 2.0 | 0.00 | 75.09 / 78.21 | 254.46 (B01) | open |
| B08 (23.2 m) | 1.0 | 0.00 | 190.98 / 165.55 | 254.50 (B01) | open |
| B11 (30.4 m) | 1.0 | 0.00 | 143.13 / 116.91 | 253.01 (B01) | open |
| B01 (6.4 m) | 0.5 | 198.16 | 241.41 (B02) | 246.26 (B44) | elastic |
| B01 (6.4 m) | 1.0 | 141.70 | 271.06 (B02) | 271.06 (B02) | elastic |
| B44 (109.6 m) | 1.0 | 135.58 | 222.78 (B43) | 254.63 (B01) | elastic |

Table: Table 6.15 – Effect of a single initial gap.

A gap in the midbody reduces the load of the gapped block, and once it exceeds the local elastic compression (0.33 mm at B21) the block carries nothing. Its load moves mainly to the two neighbours, which gain about 10 t each; the effect stays within roughly one characteristic length, and the largest reaction changes by less than 0.1 %. Near the ends the behaviour differs. A 1.0 mm gap at B01 reduces its reaction from 254.63 t to 141.70 t but shifts the load to B02, which then carries 271.06 t, more than B01 carried without the gap: a gap at the most loaded block moves the peak inwards and can increase it. A 1.0 mm gap at B44 shifts load to B43 and raises the largest bending moment from 1 723.4 to 2 103.8 t·m, because the effective forward overhang becomes longer. A 1.0 mm gap at B08, under the machinery, raises its aft neighbour to 190.98 t, close to the block capacity. Height deviations of the order of one millimetre are thus significant compared with the elastic compressions of 0.3–1.2 mm (Table 6.7). This sensitivity is the key to the welded hull of Section 6.7, whose keel deviates by tens of millimetres.

### 6.6.6 Ranking of the parameters

Table 6.16 ranks the parameters by their influence on the largest block reaction over the ranges studied.

| Rank | Parameter | Range studied | Max. $R$ range [t] | Character of the effect |
|---:|---|---|---|---|
| 1 | Extent of the block line (overhangs) | ends at 4.0 / 6.4 / 10.0 m and 100.0 / 108.4 – 109.6 / 114.4 m | 174.3 – 393.5 | dominant, local to each end |
| 2 | Uniform spacing / number of blocks | 1.6 – 5.4 m (65 – 20 blocks) | 211.4 – 397.2 | global, roughly proportional to the spacing |
| 3 | Block stiffness | 500 – 8000 kN/mm | 229.3 – 339.7 | through $\lambda$; stiffer blocks concentrate load |
| 4 | Hull bending stiffness | 0.25 – 100 × baseline | 203.7 – 336.6 | through $\lambda$; a stiffer hull spreads load |
| 5 | Single initial gap | 0 – 2 mm at various blocks | 246.3 – 271.1 | local; significant near the ends only |

Table: Table 6.16 – Ranking of the parameters by their influence on the largest block reaction (straight hull, uniform spacing).

The parameters that can be chosen in the docking plan — the extent of the block line and the positions of the blocks — are more effective than the structural parameters, which cannot be changed. Because $EI$ and $k$ enter mainly through $k/(s\,EI)$, an uncertainty in one of them can be expressed as an equivalent uncertainty in the other; quartering the hull stiffness or quadrupling the block stiffness raised the largest reaction by about one third. Welding deformation, a distributed gap, is treated in Section 6.7.

## 6.7 Welded hull: as-built reactions and identification from a survey

### 6.7.1 Case description

In `welding_case.yaml` the ship of Section 6.1 is assembled from six construction blocks of 20 m length, with erection joints at $x$ = 20, 40, 60, 80 and 100 m. Two welding-deformation mechanisms are prescribed (Section 3.9): angular distortion at the joints of 0.9, 0.5, 0.4, 0.5 and 0.6 mrad (aft to forward), positive values making the keel concave upwards; and a uniform shrinkage curvature in every construction block from a deck shrinkage of 0.24 mm/m, a bottom shrinkage of 0.20 mm/m and a depth of 10 m, i.e. $\kappa = 4 \times 10^{-6}$ m⁻¹. These values are placeholders of plausible magnitude, not measurements. The block plan is that of `example_case.yaml` (44 blocks of 200 t, 2.4 m). The keel lift is referred to the least-squares line through the block positions and shifted so that the lowest block point has zero gap (Section 3.9). The resulting gap is 41.4 mm at B01 and 39.0 mm at B44 and vanishes at B23 ($x$ = 59.2 m).

### 6.7.2 As-built block reactions

Table 6.17 lists the reactions of selected blocks.

| Block | $x$ [m] | Gap [mm] | Rigid $R$ [t] | Elastic $R$ [t] |
|---|---:|---:|---:|---:|
| B01 | 6.4 | 41.42 | 0.00 | 0.00 |
| B02 | 8.8 | 37.26 | 0.00 | 0.00 |
| B03 | 11.2 | 33.12 | 0.00 | 0.00 |
| B04 | 13.6 | 29.00 | 0.00 | 9.20 |
| B06 | 18.4 | 20.83 | 0.00 | 338.80 |
| B07 | 20.8 | 17.50 | 0.00 | 394.25 (crushed) |
| B08 | 23.2 | 15.63 | 0.00 | 317.06 |
| B09 | 25.6 | 13.79 | 0.00 | 219.40 |
| B18 | 47.2 | 1.83 | 401.69 (crushed) | 16.62 |
| B19 | 49.6 | 1.42 | 401.89 (crushed) | 0.00 |
| B20 | 52.0 | 1.03 | 401.85 (crushed) | 0.00 |
| B21 | 54.4 | 0.66 | 401.57 (crushed) | 21.93 |
| B41 | 102.4 | 27.76 | 0.00 | 32.80 |
| B42 | 104.8 | 31.48 | 0.00 | 0.00 |
| B43 | 107.2 | 35.22 | 0.00 | 0.00 |
| B44 | 109.6 | 38.98 | 0.00 | 0.00 |

Table: Table 6.17 – Selected block reactions of the welded hull as built (`welding_case.yaml`, given block plan).

With the elastic model, the welding deformation unloads eight blocks completely: the three aftmost (B01–B03), the three foremost (B42–B44) and two midship blocks (B19 and B20), where the kinks at 40 m and 60 m leave the keel locally high relative to its neighbours. The load moves to the quarter regions. B07 ($x$ = 20.8 m), just forward of the joint at 20 m and under the machinery, carries 394.25 t, exceeds its crush force and is the only crushed block; B06 and B08 carry 338.8 t and 317.1 t. The largest ratio of reaction to strength is 1.971 and the largest bearing-pressure utilisation 2.014, against 1.273 and 1.301 for the straight hull. The contact iteration needed three iterations, and the equilibrium residuals stayed below $10^{-13}$.

The rigid-hull model gives a qualitatively different and physically meaningless picture. A rigid hull cannot bend to follow a concave keel, so it rests on the lowest part of the line only: 31 blocks are unloaded, the remaining 13 blocks between 35.2 m and 64.0 m carry the whole ship, 10 of them crush and the largest reaction is 401.9 t at B19, a block that the elastic model finds unloaded. The corresponding hogging moment of −54 211.8 t·m (bending-moment utilisation 2.314) is an artefact of the rigid assumption; the elastic model gives −13 119.9 t·m at $x$ = 40.0 m, with a governing bending-moment utilisation of 0.611 and a shear-force utilisation of 0.715. For a keel deviation of several centimetres, the support problem can only be analysed with a flexible hull girder: the elastic settlement of the overhangs under their own weight, of the same order as the welding lift, decides which blocks touch.

### 6.7.3 Identification of the joint angles from a synthetic survey

The welding deformation is normally measured rather than predicted. To test the identification of Section 3.10, a synthetic survey was generated from the joint kinks of Section 6.7.1 with stations every 5 m from 0 to 120 m (25 points), an arbitrary datum and tilt and normally distributed noise with a standard deviation of 1 mm; the five joint angles and the straight-line terms were fitted by linear least squares. The survey is synthetic and does not originate from a real ship. Figure 6.7 shows the fit and the identified angles, and Table 6.18 summarises variants of the survey. Each row of the table uses the same random seed as Figure 6.7; the standard deviations were obtained from 200 further noise realisations per row.

![Figure 6.7: Identification of the joint angles from a synthetic keel survey (stations every 5 m, 1 mm noise). Left: survey points with the fitted straight line removed and the identified kink model; right: true and identified joint angles.](figures/fig_identification.png)

| Stations [m] | Points | Noise [mm] | Regularisation | Identified angles [mrad] (true 0.9 / 0.5 / 0.4 / 0.5 / 0.6) | Max. error [mrad] | RMS residual [mm] | Std. dev. (200 runs) [mrad] |
|---:|---:|---:|---:|---|---:|---:|---|
| 5.0 | 25 | 1.0 | 0 | 1.148 / 0.301 / 0.468 / 0.512 / 0.598 | 0.248 | 0.756 | 0.087 – 0.098 |
| 2.5 | 49 | 1.0 | 0 | 0.843 / 0.552 / 0.456 / 0.380 / 0.704 | 0.120 | 0.900 | 0.062 – 0.076 |
| 1.0 | 121 | 1.0 | 0 | 0.972 / 0.477 / 0.420 / 0.454 / 0.637 | 0.072 | 0.895 | 0.039 – 0.049 |
| 5.0 | 25 | 0.5 | 0 | 1.024 / 0.400 / 0.434 / 0.506 / 0.599 | 0.124 | 0.378 | 0.044 – 0.049 |
| 5.0 | 25 | 1.0 | 10⁻³ | 1.032 / 0.374 / 0.439 / 0.533 / 0.558 | 0.132 | 0.799 | 0.069 – 0.081 |
| 5.0 | 25 | 1.0 | 10⁻² | 0.640 / 0.518 / 0.463 / 0.529 / 0.405 | 0.260 | 1.486 | 0.030 – 0.037 |

Table: Table 6.18 – Identification of the joint angles from synthetic keel surveys (`identify_joint_angles`; one realisation shown per row, standard deviations over 200 realisations).

The fit reproduces the survey with an RMS residual of 0.76 mm, at the level of the noise, but the individual angles of this realisation deviate by up to 0.25 mrad, most at the joints at 20 m and 40 m. The realisation of Figure 6.7 is an unfavourable one: over 200 realisations the standard deviation is about 0.09 mrad per angle for 5 m stations, falling to 0.06–0.08 mrad for 2.5 m and 0.04–0.05 mrad for 1 m stations, and the mean error (bias) is below 0.01 mrad. Halving the noise halves the scatter, as expected for a linear estimator. The problem is moderately ill-conditioned [@hansen1998]: noise near a joint can be explained either by that joint or by its neighbour. The sum of the angles, which governs the end lift, is better determined (3.03 mrad identified against 2.90 mrad true). Tikhonov regularisation trades scatter for bias: $10^{-3}$ reduces the scatter slightly and halves the largest error of this realisation, while $10^{-2}$ biases the outer angles towards zero (by −0.34 mrad at 20 m) and doubles the residual, the usual trade-off between noise suppression and bias [@hansen1998]. For the support problem the identified model is used only to compute the gaps at the blocks, whose error is of the order of the survey noise.

### 6.7.4 Survey-based welding case

`welding_survey_case.yaml` demonstrates the complete data path from a survey file to the block reactions. The survey `welding_survey.csv` contains 25 keel heights every 5 m. It was generated from the welding model of Section 6.7.1, shifted by an arbitrary instrument datum and tilt, disturbed by noise with a standard deviation of 0.5 mm and rounded to 0.1 mm. The case gives the joint positions and the zone shrinkage and requests the identification of the joint angles; the known shrinkage part is subtracted from the survey before the fit (Section 3.10).

The identified angles are 0.938, 0.477, 0.384, 0.497 and 0.670 mrad, against 0.9, 0.5, 0.4, 0.5 and 0.6 mrad used to generate the data, i.e. within 0.07 mrad. The survey is reproduced with an RMS residual of 0.295 mm and a largest residual of 0.698 mm, consistent with the noise. The gaps at the end blocks become 41.04 mm aft and 38.89 mm forward, against 41.42 mm and 38.98 mm for the prescribed model. The elastic as-built state is practically the same as in Section 6.7.2: nine blocks carry no load (B01–B04, B19, B20, B42–B44), B07 is crushed with 395.10 t, and the largest ratio of reaction to strength is 1.975 (1.971 with the prescribed model).

Figure 6.8 shows the welding data figure that the program produces for every case with welding deformation (file `welding_data.png`, GUI result tab "Welding data"). The top panel shows the keel lift used for the gaps, its components — the joint kinks and the zone shrinkage, each with its straight-line part removed — the 25 survey points and the block positions. The middle panel shows the identified joint angles (0.94, 0.48, 0.38, 0.50 and 0.67 mrad) and the uniform zone curvature of $4 \times 10^{-6}$ m⁻¹. The bottom panel shows the gap from welding at every block, with the nine blocks without contact as built highlighted. The figure makes the cause of the lift-off directly visible: the joint kinks dominate the keel shape, and the end blocks lie 28–41 mm below the keel.

![Figure 6.8: Welding deformation data of the survey-based case (synthetic survey). Top: keel lift with its components, the survey points and the block positions; middle: identified joint angles and zone shrinkage curvature; bottom: gap from welding at every block, blocks without contact as built highlighted.](figures/fig_welding_data.png)

## 6.8 Automatic spacing design from block strengths

This section presents the final result of the thesis. For a block plan with a given number of blocks and a given strength of every block, the program computes the spacing of the blocks at which every reaction is within the strength of its block, by the three-phase method of Section 3.12: strength-weighted equidistribution (phase 1) [@huang2011], sequential linear programming on the minimax problem with a trust region (phase 2) [@griffith1961; @conn2000], and rounding to the floors with a local search (phase 3). Every candidate layout is analysed with the full support model, including, for a welded hull, the welding gaps at the candidate block positions (Section 3.12.2). The quantity minimised is the largest utilisation

$$ U(\mathbf{x}) = \max_i \frac{R_i(\mathbf{x})}{S_i} \tag{6.3} $$

Before the results, a simple bound is useful for judging them. Since every reaction satisfies $R_i \le U S_i$ and the reactions sum to the weight, $W = \sum_i R_i \le U \sum_i S_i$, so that for *any* spacing

$$ U(\mathbf{x}) \;\ge\; \bar u \equiv \frac{W}{\sum_{i=1}^{N} S_i}, \tag{6.4} $$

with equality only if every block is used to the same fraction of its strength. The bound $\bar u$ is also the mean utilisation $\sum R_i / \sum S_i$. A plan can fit only if $\bar u \le u_{target}$; for $N$ blocks of equal strength $S$ this gives the smallest possible number of blocks

$$ N_{\min} = \left\lceil \frac{W}{u_{target}\, S} \right\rceil . \tag{6.5} $$

The bound is not attainable in general, because the overhangs, the extent of the block line and the spacing limits prevent an equal utilisation; the closed-form case of Section 3.12.7, in which it is attained, was used for verification (Section 5.6). The ratio $U/\bar u$ is used below as a measure of how evenly a design uses its blocks.

### 6.8.1 Spacing design of the example block plan

The block plan of `spacing_case.yaml` has 44 blocks with $\sum S_i$ = 8090 t, so $\bar u$ = 4500/8090 = 0.556. On the given uniform spacing of 2.4 m the elastic model loads B01 to 1.340 times its strength and seven blocks beyond their strength (Section 6.5). Figure 6.9 shows the design, and Table 6.19 its progress by phase.

![Figure 6.9: Spacing design for the block plan of `spacing_case.yaml` (44 blocks, strengths 190 t and 160 t, elastic hull). Top: weight distribution with the given (light) and designed (dark) block positions; second: given and designed spacing to the next block along the hull; third: block reactions of the designed spacing (bars) against the block strengths (dashed), with the reactions at the given spacing as markers; bottom: largest ratio of reaction to strength over the number of analyses in the three phases.](figures/fig_spacing_design.png)

| Phase | Step | Analyses (cumulative) | Largest $R_i/S_i$ |
|---|---|---:|---:|
| start | given plan | 1 | 1.340 |
| 1 – equidistribution | iteration 1 | 3 | 1.041 |
| | iteration 2 | 4 | 0.881 |
| | iteration 3 | 5 | 0.807 |
| | iteration 4 | 6 | 0.788 |
| | best (iteration 5); iterations 6–15 stationary at 0.788 | 7 (17) | 0.787 |
| 2 – sequential LP | accepted step 1 | 50 | 0.778 |
| | accepted step 2 | 87 | 0.770 |
| | accepted step 3 | 122 | 0.762 |
| | accepted step 5 (step 4 rejected) | 156 | 0.761 |
| 3 – floors and local search | final | 193 | 0.761 |

Table: Table 6.19 – Convergence of the spacing design of `spacing_case.yaml` by phase.

**Convergence.** Phase 1 does most of the work: four fixed-point steps reduce the largest ratio from 1.340 to 0.788 — below the target already after the second step (analysis 4) — and the iteration then becomes stationary. Its fixed point is not the minimax optimum, because equidistribution balances an estimated load density between neighbours and not the reactions of the elastic hull themselves (Section 3.12.3). Phase 2 removes the remaining 3.5 % in four accepted steps; each sequential-LP iteration costs up to one finite-difference analysis per block plus one trial analysis, 33–37 analyses here, which is why phase 2 accounts for most of the 193 analyses. Phase 3 rounds the positions to the floors and searches the neighbouring floor positions of the governing blocks; here it confirmed the continuous result (0.761) with 37 analyses. The whole run, including the rigid and elastic analyses of the given plan and the post-processing, took about 0.7 s.

**Designed spacing.** The designed spacing follows the load and the strengths (Figure 6.9, top). Aft of 14.4 m, under the stern overhang, the superstructure and the start of the engine room, nine blocks are set at spacings of 0.8 and 1.6 m; the nine 160 t blocks B04–B12 now lie between 9.6 m and 20.0 m. Between 14.4 m and 40.8 m, under the engine room, the spacing is 1.6–2.4 m. Amidships it opens to 3.2 m (40.8–63.2 m), then 4.0 m and 4.8 m up to 105.6 m, where the elastic hull spreads the load most evenly; and the six foremost blocks are set at 0.8 m to share the bow overhang. The pattern alternates between 0.8 and 1.6 m in the aft part because the positions must lie at floors 0.8 m apart; without rounding (Table 6.21) the same design consists of smooth runs of 1.045, 1.84, 3.185 and 4.8 m with the same largest ratio, the effect of the smoothing term in the objective of phase 2 (Section 3.12.4).

**Reactions versus strength.** The designed reactions lie between 73.8 t and 144.6 t (Figure 6.9, middle). Five blocks are within 2.5 % of the largest ratio: B44 (0.761), B01 (0.760), B32 at 72.0 m (0.750), B12 at 20.0 m (0.745, a 160 t block) and B33 (0.744). Several active peaks of nearly equal height are the characteristic signature of a minimax optimum: moving any block to relieve one of them loads another. The least used block is B22 at 38.4 m (0.388). The end blocks are among the governing ones because they sit at the limits of the allowed line (6.4 and 109.6 m) and must carry the overhangs; this is why the design stays at $U/\bar u$ = 1.37 and why a longer line helps (Section 6.8.2).

**All checks.** Table 6.20 compares the checks of the elastic model for the given and the designed spacing.

| Check | Given spacing: $\eta$ | Designed spacing: $\eta$ |
|---|---:|---:|
| Hull bending moment | 0.170 | 0.176 |
| Hull shear force | 0.336 | 0.336 |
| Deck bending stress | 0.088 | 0.087 |
| Keel bending stress | 0.069 | 0.068 |
| Hull shear stress | 0.067 | 0.067 |
| Block bearing pressure | 1.301 | 0.738 |
| Block load vs strength | 1.340 | 0.761 |
| All checks | failed | **passed** |

Table: Table 6.20 – Checks of the elastic model for the given and the designed spacing of `spacing_case.yaml`.

Every block of the designed plan is within its strength, the largest bearing pressure is 1.477 MPa, no block is open or crushed, and all hull girder checks pass. The hull girder loads hardly change (largest moment 1 700.8 t·m instead of 1 723.4 t·m), because they are governed by the overhangs, which the design does not change. Every block lies at a floor, so no off-structure warning is raised. The designed spacings are written to `spacing.csv` and, in the GUI, back into the block plan.

### 6.8.2 Variants of the block plan and the settings

Table 6.21 shows how the result depends on the settings and on the block plan. All variants were computed with the current program from `spacing_case.yaml` with the changes listed; the given plans have a uniform spacing of 2.4 m from $x$ = 6.4 m, so that the "given" ratio of the variants with fewer blocks refers to a block line that ends well before 109.6 m.

| Variant | Blocks × strength | Extent [m] | $s_{\max}$ [m] | Analyses | $U$ given → designed | Mean $= \bar u$ / min. $R_i/S_i$ | Designed spacing [m] | Result |
|---|---|---|---:|---:|---|---|---|---|
| example (floors) | 44, 190 / 160 | 6.4–109.6 | 4.8 | 193 | 1.340 → 0.761 | 0.558 / 0.388 | 0.8–4.8 | fits |
| without rounding to floors | 44, 190 / 160 | 6.4–109.6 | 4.8 | 188 | 1.340 → 0.761 | 0.558 / 0.382 | 0.8–4.8 | fits |
| smaller $s_{\max}$ | 44, 190 / 160 | 6.4–109.6 | 3.2 | 251 | 1.340 → 0.760 | 0.562 / 0.296 | 0.8–3.2 | fits |
| $s_{\max}$ = given spacing | 44, 190 / 160 | 6.4–109.6 | 2.4 | 29 | 1.340 → 1.310 | 0.567 / 0.274 | 0.8–2.4 | does not fit |
| longer line | 44, 190 / 160 | 4.0–114.4 | 4.8 | 502 | 1.340 → 0.575 | 0.557 / 0.493 | 0.8–4.8 | fits |
| rigid hull model | 44, 190 / 160 | 6.4–109.6 | 4.8 | 242 | 0.869 → 0.600 | 0.558 / 0.497 | 0.8–4.0 | fits |
| fewer blocks | 40 × 190 | 6.4–109.6 | 4.8 | 407 | 2.071 → 0.749 | 0.592 / 0.365 | 0.8–4.8 | fits |
| fewer blocks, longer line | 32 × 190 | 4.0–114.4 | 6.4 | 402 | 2.585 → 0.751 | 0.740 / 0.721 | 0.8–6.4 | fits |
| minimum number | 29 × 160 | 4.0–114.4 | 6.4 | 220 | 3.547 → 0.984 | 0.970 / 0.933 | 0.8–6.4 | fits |
| one block fewer | 28 × 160 | 4.0–114.4 | 6.4 | 421 | 3.725 → 1.019 | 1.004 / 0.968 | 0.8–6.4 | does not fit (Σ S = 4480 t < W) |
| three blocks fewer | 26 × 160 | 4.0–114.4 | 6.4 | 198 | 4.115 → 1.114 | 1.082 / 1.007 | 1.6–6.4 | does not fit (Σ S = 4160 t < W) |

Table: Table 6.21 – Spacing design for variants of the block plan of the example vessel (elastic hull unless stated; $s_{\min}$ = 0.8 m; positions at the floors unless stated). "Mean" is the mean ratio of reaction to strength of the designed plan, which equals the bound $\bar u$ of Equation (6.4) up to round-off.

The variants support five conclusions.

1. **Rounding to the floors costs nothing here.** With and without rounding the largest ratio is 0.761; the floors restrict where the blocks may stand but not the quality of the design, because the floor spacing (0.8 m) is small compared with the designed spacings.

2. **The maximum spacing must leave freedom, but not much.** With $s_{\max}$ = 3.2 m the design reaches the same largest ratio (0.760) by a different pattern: it cannot open the midbody as far, so it places more blocks there (25 spacings of 3.2 m) and uses the least loaded block less (0.296). With $s_{\max}$ equal to the given spacing of 2.4 m, the 43 spacings of at most 2.4 m exactly span the line, no redistribution is possible, and phases 1 and 2 cannot move; only the local search shifts the end blocks inwards by two floors, which reduces the ratio to 1.310. The program reports that the plan does not fit, that the aft end block at the limit of the line governs, and suggests a longer line or about 58 blocks.

3. **The extent matters as much as the spacing.** With the line extended to 4.0–114.4 m, as in Section 6.6.2, the same 44 blocks reach 0.575, only 3.4 % above the bound $\bar u$ = 0.556, and the least used block carries 0.493 of its strength: the blocks are used almost evenly. The remaining distance to the bound in the example plan (0.761 against 0.556) is therefore due to the extent of the block line, a physical constraint, and not to the algorithm. This run needed 502 analyses because the local search of phase 3 found several improving moves after rounding (0.568 continuous, 0.575 at the floors).

4. **The number of blocks can be reduced to the minimum of Equation (6.5).** 40 blocks of 190 t fit on the original line (0.749), and 32 blocks of 190 t on the longer line with $s_{\max}$ = 6.4 m (0.751, only 1.5 % above $\bar u$ = 0.740). For blocks of 160 t, Equation (6.5) gives $N_{\min} = \lceil 4500/160 \rceil = 29$. The design found a spacing for exactly 29 blocks with a largest ratio of 0.984, 1.5 % above the bound $\bar u$ = 4500/4640 = 0.970, and with every block used to between 93 % and 98 % of its strength: a nearly uniform utilisation, the ideal of a support fitted to its blocks. With 28 blocks the total strength of 4480 t is less than the load, no spacing can fit, and the program says so, gives the reason (Σ S < W) and suggests 29 blocks; the same holds for 26 blocks. The design thus brackets the minimum: 29 fits and 28 provably cannot.

5. **The model matters.** The rigid-hull model, which spreads the load more evenly (Section 6.3), gives a lower ratio for the given plan (0.869) and a different design (0.600), with a dense spacing of 0.8 m over the aft 10 m of the line. Analysed with the elastic model, this rigid-designed layout gives a largest ratio of 0.967 at B44 instead of the 0.600 predicted: it still fits, by a margin of 3 %, but far less than the elastic design's 0.761. Since the end-block loads are the main effect of hull flexibility, design decisions should be based on the elastic model.

### 6.8.3 Weight-curve case

The spacing design of `weight_curve_case.yaml` uses 44 blocks of 200 t ($\sum S_i$ = 8800 t, $\bar u$ = 4643/8800 = 0.528) with the same settings as the example. On the given plan the elastic model loads B44 to 1.064 of its strength (Section 6.3.4). The design brings the largest ratio to 0.568 in 501 analyses (about 1.9 s for the complete run): phase 1 reaches 0.629 in 15 iterations, phase 2 0.567 in six accepted steps, and rounding to the floors with the local search ends at 0.568. All checks pass (bending moment 0.118, shear force 0.267, bearing pressure 0.580, strength 0.568).

The designed spacing is 0.8–4.8 m. Over 34 blocks from 6.4 m to 92.0 m, where the weight curve is smooth and nearly level, the utilisation is almost constant (0.543–0.568) at spacings of mostly 1.6–3.2 m, with closer spacing (0.8–2.4 m) under the engine room (16–28 m), where the curve reaches 64 t/m. The eight foremost blocks are packed at the minimum spacing of 0.8 m between 104.0 m and 109.6 m to carry the forward overhang of 243.3 t. Their utilisation falls from 0.568 at B44 to about 0.33 four to five blocks aft of it, because the end block at the line limit carries most of the overhang and its neighbours cannot be moved closer than $s_{\min}$; the forward end, constrained by the extent and the minimum spacing, therefore keeps the design 7.6 % above $\bar u$. The case confirms that the spacing design works unchanged with a weight curve and items given by moment, because both are converted into the same load model (Section 3.3).

### 6.8.4 Welding deformation as a constraint

For a welded hull, every candidate layout of the spacing design is analysed with the welding gaps recomputed at its own block positions (Section 3.12.2). The deformation is therefore a constraint of the design: a block moved along the hull also moves along the curved keel and changes its gap. To find how much deformation the spacing alone can compensate, the welding deformation of `welding_case.yaml` was scaled by a factor (all joint angles and the difference of deck and bottom shrinkage), and the spacing design was run for each factor with 44 blocks of 200 t, given spacing 2.4 m, $s_{\min}$ = 0.8 m, $s_{\max}$ = 4.8 m and the line 6.4–109.6 m. Table 6.22 and Figure 6.10 show the result. The column "continuous" gives the result of the same design without rounding to the floors.

| Factor | End lift [mm] | $U$ given | $U$ designed (floors) | $U$ designed (continuous) | Result (floors) |
|---:|---:|---:|---:|---:|---|
| 0.00 | 0.0 | 1.273 | 0.701 | 0.699 | fits |
| 0.10 | 4.1 | 1.081 | 0.611 | 0.602 | fits |
| 0.25 | 10.4 | 1.188 | 0.710 | 0.657 | fits |
| 0.40 | 16.6 | 1.404 | 0.892 | 0.816 | fits |
| 0.50 | 20.7 | 1.547 | 0.917 | 0.915 | fits |
| 0.55 | 22.8 | 1.619 | 1.045 | 0.950 | does not fit |
| 0.60 | 24.9 | 1.691 | 1.111 | 0.997 | does not fit |
| 0.75 | 31.1 | 1.904 | 1.266 | 1.155 | does not fit (1 block open) |
| 1.00 | 41.4 | 1.971 | 1.350 | 1.306 | does not fit (5 blocks open) |

Table: Table 6.22 – Spacing design under scaled welding deformation (44 blocks of 200 t, elastic hull).

![Figure 6.10: Largest ratio of reaction to strength at the given spacing (2.4 m) and at the designed spacing (at the floors) as a function of the welding end lift of the keel (44 blocks of 200 t, elastic hull).](figures/fig_welding_spacing_study.png)

A small end lift helps: at 4.1 mm the given plan improves from 1.273 to 1.081 and the design reaches 0.611, better than for the straight hull (0.701), because a slightly concave keel counteracts the load concentration at the overhangs found in Section 6.3. Beyond about 5 mm the lift works against the design, and the designed ratio rises steadily. **With the positions at the floors, the spacing design fits up to an end lift of about 21 mm (factor 0.5, ratio 0.917) and not beyond** (22.8 mm: 1.045). Without the floor constraint the threshold lies at about 25 mm (factor 0.6, 0.997). The rounding, which cost nothing for the straight hull, costs up to 0.11 here (0.950 against 1.045 at factor 0.55): with large gaps the reactions depend on the positions in a strongly non-smooth way, so that a shift of one floor can bring a block into or out of contact.

**The full deformation.** Figure 6.11 shows the design for the full deformation of `welding_case.yaml` (end lift 41.4 mm). The run took 78 analyses and 0.5 s. Phase 1 reduced the ratio from 1.971 to 1.309 within six analyses, but the fixed-point iteration then diverged to about 1.59 (Figure 6.11, bottom): with gaps the reactions are not a smooth function of the positions, and the load-density estimate on which equidistribution relies no longer describes the redistribution; the best layout found was kept, as the algorithm prescribes. Phase 2 reached 1.306 in one accepted step, after five trial steps had been rejected and the trust region reduced, and rounding to the floors gave 1.350. In the designed layout five blocks carry no load (B01 at 6.4 m, B02 at 11.2 m, B26 at 50.4 m, B43 at 104.8 m and B44 at 109.6 m), and six blocks exceed their strength, all of them at or next to an erection joint: 269.9 t at 20.8 m (the governing block, 1.350), 241.0 t at 19.2 m, 225.3 t at 21.6 m, 210.7 t at 60.0 m, 216.7 t at 80.0 m and 236.9 t at 100.0 m. The design has moved blocks towards the joints — dense groups at 0.8 m around 18.4–29.6 m and 79.2–83.2 m — but the governing block stays at the joint at 20 m.

![Figure 6.11: Spacing design of the welded hull with the full welding deformation (`welding_case.yaml` with the spacing design enabled; end lift 41.4 mm). Top: weight distribution with the given and designed block positions; second: given and designed spacing; third: reactions of the designed spacing against the strength of 200 t (red: above strength), with the reactions at the given spacing as markers; bottom: largest ratio of reaction to strength over the analyses.](figures/fig_spacing_welding.png)

**More blocks do not help.** Repeating the design with more blocks on the same line gives 1.350 for 44 blocks, 1.382 for 52, 1.385 for 60, 1.340 for 68 and 1.351 for 76 blocks (Table 6.23). In every case the governing block lies at the 20 m joint, and between 5 and 14 blocks are out of contact. Adding blocks does not reduce the governing reaction, unlike for the straight hull (Table 6.11).

| Blocks | 44 | 52 | 60 | 68 | 76 |
|---|---:|---:|---:|---:|---:|
| $U$ designed | 1.350 | 1.382 | 1.385 | 1.340 | 1.351 |

Table: Table 6.23 – Spacing design for the full welding deformation with more blocks on the line 6.4–109.6 m (governing block always at the erection joint at 20 m; 5–14 blocks open).

**The program's diagnosis.** The program reports the case as "does not fit" and explains why (Section 3.12.6): the best spacing gives a largest ratio of 1.350, five blocks carry no load because of gaps under the keel, "so the gaps rather than the spacing govern and more blocks may not help", and it recommends adjusting the block heights to the keel line or using stronger blocks where the keel bears. The suggested number of blocks (about 60, from $\lceil N U / u_{target} \rceil$) is explicitly qualified as valid only "if the load spread", which the study of Table 6.23 confirms it does not.

**Physical interpretation.** A concave-up kink makes the keel at an erection joint a local low point: on both sides the keel rises linearly away from the joint. The block under the joint is therefore compressed first, and a neighbouring block at a distance $d$ comes into contact only when the hull has bent, or the joint block has been compressed, by the additional gap of order $\theta d$, where $\theta$ is the local change of slope. A block reaches its strength at the compression

$$ \delta_S = \frac{S}{k} = \frac{200 \times 9.80665\ \text{kN}}{2000\ \text{kN/mm}} \approx 0.98\ \text{mm}, \tag{6.6} $$

while a joint angle of 0.9 mrad opens an extra gap of $\theta d$ ≈ 0.7 mm over a single floor spacing ($d$ = 0.8 m) — of the same order. Only blocks within one or two floor spacings of the joint can therefore share the load before the joint block reaches its strength, and how much load collects at the joint is determined by the bending of the hull between the joints, not by the number of blocks. The joint acts as a *hard point*: blocks placed beside it stand in a gap and stay idle, and blocks placed on it overload. Spacing changes the positions of the supports, but it cannot change the shape of the keel relative to the block tops. Beyond an end lift of about 21 mm, therefore, the spacing design must be combined with an adjustment of the block heights (shims, wedges or adjustable blocks) to the surveyed keel line; this combination is outside the scope of the thesis and is recommended as future work (Section 7.5).

## 6.9 Discussion

This section discusses the results with respect to the research questions of Section 1.5 and states their limitations.

### 6.9.1 RQ1 – Rigid versus elastic reactions

For the straight hull, hull flexibility changes the results qualitatively. The rigid-hull method underestimated the most heavily loaded block by 43 % (145.06 t instead of 254.63 t) and the forward end block by a factor of 4.1; it passed all block checks where the elastic model found three overloaded blocks (seven with the individual strengths of the spacing case); and it overestimated the largest bending moment by a factor of 7.1. The same pattern appeared for the weight-curve case (113.50 t against 212.80 t at the forward end block). The characteristic length of 11.9 m (Equation 6.1) is short compared with the 103 m block line, which explains why the hull behaves as a beam on an elastic foundation [@hetenyi1946] and why the linear distribution of the rigid method [@nstm997] cannot capture the load concentration at the overhangs. The difference occurs at the ends of the block line and, to a lesser extent, under concentrated weights such as the machinery. Even with a hundredfold bending stiffness the end-block reaction remained 40 % above the rigid value, because the shear flexibility remains. The rigid method remains useful as an independent equilibrium check and quick estimate, but not for verifying block loads where overhangs exist — and, consequently, not as the model of the spacing design (Section 6.8.2, item 5).

### 6.9.2 RQ2 – Effect of welding deformation; identification

The welding deformation of a block-built hull acts as a smooth, concave-up distribution of initial gaps of tens of millimetres, an order of magnitude larger than the elastic block compressions (about 1 mm). Small end lifts (about 4–10 mm) are absorbed by the settlement of the overhangs and even relieve the end blocks. Larger lifts unload the end blocks and then further blocks, and concentrate the load at the erection joints in the quarter regions: at 41 mm end lift, eight blocks carried no load, one crushed, and the largest ratio of reaction to strength rose from 1.273 to 1.971. Which blocks lose contact depends on the keel shape, the weight distribution and the hull stiffness together; two midship blocks were unloaded because of the local shape of the kinked keel, which no simple rule would predict. The rigid model is meaningless for such a hull. The joint angles can be identified from a keel survey by linear least squares: with 1 mm noise the residual equalled the noise and the angles had a standard deviation of 0.04–0.09 mrad depending on the station spacing; for the sample survey with 0.5 mm noise the angles were recovered within 0.07 mrad and the as-built results matched those of the prescribed model (largest ratio 1.975 against 1.971).

### 6.9.3 RQ3 – Automatic spacing design: feasibility, optimality, efficiency, impossibility

**Feasibility.** Yes, the spacing can be computed automatically such that every reaction is within its block strength whenever such a spacing exists within the limits. For the example plan the design turned a layout that overloaded seven blocks (largest ratio 1.340) into one in which every block uses at most 76.1 % of its strength, by changing only the positions, with all blocks at floors and all hull girder and block checks passed. The same held for the weight-curve case (1.064 → 0.568) and for all variants of Table 6.21 whose total strength exceeds the load and whose spacing limits leave freedom.

**Closeness to the optimum.** The design problem is non-convex and non-smooth, so the method cannot prove global optimality (Section 1.4). Three observations indicate that it comes close. First, where the extent of the line does not bind, the result lies within 1.5–3.4 % of the lower bound $\bar u$ of Equation (6.4), which no spacing can beat: 0.575 against 0.556 for the longer line, 0.751 against 0.740 for 32 blocks, 0.984 against 0.970 for 29 blocks. Second, the designed layouts show several active peaks of nearly equal height, the signature of a minimax optimum. Third, the method reproduced the closed-form optimum of Section 3.12.7 in the verification (Section 5.6). Where the result lies far above the bound (0.761 against 0.556), the cause is a binding physical constraint — the end blocks at the limits of the line carrying the overhangs — as the longer-line variant demonstrates.

**Efficiency.** A design needed 193–502 analyses for the example vessel, and between 0.7 s and about 2 s, which is interactive speed. The cost is dominated by the finite-difference sensitivities of phase 2 (one analysis per block and iteration), each analysis being a sparse banded solution linear in the number of elements (Section 3.13). Phase 1 is cheap — one analysis per iteration — and brought the largest ratio below the target within four analyses in the example. Phase 3 is bounded by its budget of $25N$ analyses.

**Impossibility.** The design fails, and says why, in three situations, all of which occurred: the total strength is less than the load (28 and 26 blocks of 160 t; Equation 6.4 shows that no spacing can then fit); the spacing limits leave no freedom and an end block at the limit of the line governs ($s_{\max}$ = 2.4 m; the remedy is a longer line or more blocks); and blocks are out of contact because of gaps (the welded hull; more blocks do not help). In each case the program's message names the cause and a remedy, which turns a failed design into information for the docking engineer.

### 6.9.4 RQ4 – Welding deformation as a constraint of the spacing design

Because every candidate is analysed with the welding gaps at its own positions, the designed layout is valid for the welded hull and not only for an idealised straight one. For the example, the spacing alone compensated the welding deformation up to an end lift of about 21 mm (about half of the modelled deformation) with the blocks at the floors, and about 25 mm without the floor constraint. Beyond that limit no spacing fits, and more blocks do not help (1.340–1.385 for 44–76 blocks), because the erection joints act as hard points: the keel rises away from a kink faster than a block can be compressed before reaching its strength (Equation 6.6), so that load collects at the joint whatever the block count. The spacing design recognises this situation and reports that the gaps rather than the spacing govern. Two further effects were observed under large gaps: the fixed-point iteration of phase 1 can diverge (the best iterate is kept), and rounding to the floors matters more, because the reactions are no longer smooth functions of the positions. Practically, the welding limit means that a hull whose surveyed keel deviates by more than about 20 mm at the ends (for this ship and block type) needs an adjustment of the block heights in addition to the designed spacing.

### 6.9.5 RQ5 – Trustworthiness through verification and architecture

The results of this chapter are trustworthy to the extent that the program solves its equations correctly, and this was established in Chapter 5 by 131 automated tests against closed-form solutions, limiting cases and independent computations, including the closed-form optimum of the spacing design. The results also carry their own evidence: every analysis reports its equilibrium residuals (below $10^{-13}$ throughout this chapter) and the closure of shear force and bending moment; the weight check confirmed the load model to machine precision (Section 6.2); and the designed spacing is not accepted from the optimiser's linear model but verified by a full analysis of the final layout, whose checks are reported like those of the given plan (Table 6.20). The architecture contributes in two ways. Because the spacing design receives the analysis as an injected `evaluate()` function (Section 4.10), the layouts it designs are analysed by exactly the code that analyses a given plan, including the welding gaps; and because the command-line program and the GUI use the same core, the numbers in this chapter are reproduced by either interface. All results were obtained from plain-text case files and short scripts, so they can be reproduced exactly. Verification does not establish validity, however [@oberkampf2010].

### 6.9.6 Limitations

- **Illustrative data.** All input data are placeholders, including the strengths and the welding parameters, and the keel surveys are synthetic. The magnitudes are realistic, but statements about specific blocks, and the welding limit of about 21 mm, apply only to the example. The trends — end-block concentration, the effect of the extent, near-uniform utilisation of a designed spacing, the hard-point effect of joint kinks — are expected to be general for ships with overhangs and block-built hulls.
- **Two-dimensional model.** The hull girder is a beam in the vertical-longitudinal plane; transverse load sharing with side blocks, local bottom deformation and heel are not represented.
- **Simplified welding model.** The deformation is described by joint kinks and uniform zone curvatures and treated as a stress-free initial shape; thermo-mechanical effects, residual stresses and the erection sequence are not modelled.
- **Fully docked condition; fixed block heights.** The landing phase and any adjustment of block heights are outside the scope (Section 1.4). The spacing design therefore assumes block tops at the nominal line plus the given gaps.
- **Deterministic design.** The design assumes exact gaps, stiffnesses and weights. Survey errors and block height tolerances of the order found in Section 6.6.5 would change the reactions; a robust or probabilistic design was not attempted. The margin of a design ($u_{target} - U$) is the only protection against such scatter.
- **Local optimum and cost.** The method finds a good local optimum; Equation (6.4) bounds the remaining gap but does not close it. The finite-difference sensitivities cost one analysis per block and iteration, acceptable for keel block lines but limiting for much larger problems.
- **Linear-elastic hull, simplified block law, rigid dock floor.** Hull plasticity, timber creep, load history and the flexibility of a floating dock are not modelled.
- **No validation.** The software was verified (Chapter 5), but no comparison with measured block loads was possible [@oberkampf2010].

## 6.10 Summary

- The weight distribution integrated exactly to 4500.0 t at LCG 50.469 m (items) and to 4643.0 t at 56.639 m (lightship curve of 4343.0 t and deadweight items of 300.0 t given by moment). In the example, 266.3 t lie aft and 298.8 t forward of the block line.
- The elastic model predicted end-block reactions of 254.63 t (B01) and 246.26 t (B44), against 145.06 t and 59.48 t with the rigid method; three blocks of 200 t, or seven blocks of the individual-strength plan (largest ratio 1.340), were overloaded, whereas the rigid method found none. The rigid method overestimated the largest bending moment by a factor of 7.1. The characteristic length of 11.9 m explains the difference.
- The extent of the block line was the most influential parameter (4.0–114.4 m: −32 % largest reaction, all checks passed), followed by the uniform spacing, the block stiffness and the hull stiffness; single gaps of about 1 mm matter near the ends. A uniform spacing cannot fit the reactions to the strengths.
- A welding end lift of 41.4 mm unloaded eight blocks, crushed one (394.25 t) and raised the largest ratio of reaction to strength to 1.971; the rigid model was meaningless for this case. Joint angles were identified from a synthetic survey within 0.07 mrad (0.5 mm noise), with a scatter of 0.04–0.09 mrad for 1 mm noise.
- **Final result.** For the 44-block plan (190 t, nine blocks of 160 t), the spacing design reduced the largest ratio of reaction to strength from 1.340 to 0.761 in 193 analyses (about 0.7 s), with a spacing of 0.8–4.8 m at the floors and all checks passed; most of the reduction came from phase 1, the last 3.5 % from phase 2.
- On a line of 4.0–114.4 m the same blocks reached 0.575, 3.4 % above the bound $W/\sum S_i$ = 0.556. With blocks of 160 t, the minimum number $\lceil 4500/160 \rceil$ = 29 fitted with every block used to 93–98 % of its strength (0.984 against the bound 0.970), and 28 and 26 blocks were correctly reported as impossible (Σ S < W). For the weight-curve case the design gave 1.064 → 0.568.
- Under welding deformation as a constraint, the spacing design fitted up to an end lift of about 21 mm (25 mm without the floor constraint) and not beyond. For the full 41 mm the best spacing gave 1.350 (1.306 before rounding), more blocks (44–76) did not help, the governing block stayed at the erection joint at 20 m, and the program correctly diagnosed that the gaps rather than the spacing govern. Beyond this limit the block heights must be adjusted to the keel line.

These results form the basis for the conclusions and recommendations in Chapter 7.
