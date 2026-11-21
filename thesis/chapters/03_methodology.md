# 3 Methodology

This chapter sets out the mathematical models and numerical procedures developed in this thesis. Their purpose is a single result: the spacing of a given set of keel blocks, each with its own strength, such that no block reaction exceeds the strength of its block while the ship rests in dry dock. To obtain this result the chapter first builds the load model of the ship, the model of the blocks and two models for the block reactions, together with the hull girder checks and the welding deformation of the assembled hull. It then formulates the spacing design as an optimisation problem and solves it by a three-phase algorithm. Each model is derived from first principles to the level of detail needed for implementation. The software realisation is described in Chapter 4 and the verification in Chapter 5.

## 3.1 Overview of the approach

A ship supported on $N$ blocks is a statically indeterminate system whenever $N > 2$. Vertical force and moment equilibrium give only two equations for $N$ unknown reactions. The remaining $N-2$ conditions follow from compatibility: the deflected keel line and the compressed block tops must coincide at every block that is in contact. The block reactions therefore depend on:

- the longitudinal weight distribution;
- the block positions, that is, the spacing;
- the block stiffness;
- the initial fit between the block tops and the keel, which for a hull assembled from welded construction blocks is governed by the welding deformation of the hull;
- the bending and shear stiffness of the hull girder.

Of these quantities, the spacing is the one that the docking engineer can choose freely within practical limits. The methodology is therefore organised as a pipeline that ends in a design step for the spacing, shown in Figure 3.1.

![Figure 3.1: Overall calculation procedure: input data, validation, load model and welding gaps, the two support models, the contact iteration, the block reactions and the spacing design, the hull girder and block checks, the verdict and the report.](figures/fig_method_flow.png)

The boxes of Figure 3.1 correspond to the following steps:

1. **Input (Section 3.4.4 and Chapter 4).** Weight items with an LCG or a longitudinal moment and a weight-distribution curve, hull sections with bending stiffness $EI$, shear stiffness $\kappa G A$ and section moduli $Z$, block types, the block plan (spacing, strength and gap of every block), the floors, the permissible bending moment and shear force, and the welding deformation.
2. **Validation and unit conversion (Section 3.2.3).** The input is checked and converted to SI units by the loader and the data model of the program.
3. **Load model and welding gaps (Sections 3.3 and 3.9–3.10).** Weight items and the weight curve are converted into a piecewise-linear load intensity $q(x)$ on a mesh that makes all load integrals exact. In parallel, the welding deformation of the assembled hull (joint kinks, shrinkage curvature and, optionally, a survey or identified joint angles) is converted into initial block gaps $g_i$.
4. **Support models (Sections 3.4–3.6).** Each block is a compression-only, possibly crushable spring. Method A treats the hull as rigid, so that its displacement is the straight line $v(x) = a + b\,x$ and the problem reduces to a $2 \times 2$ system. Method B models the hull as a Timoshenko beam on the block springs, discretised by finite elements.
5. **Contact iteration (Sections 3.7 and 3.11).** For both models the contact state of every block (open, elastic or crushed) is found by an active-set iteration on the affine block law $R = k_t\,\delta + c$, protected against cycling and embedded in a gap continuation when the gaps are large.
6. **Block reactions and spacing design (Section 3.12).** The converged reactions $R_i$ are compared with the block strengths $S_i$. The spacing design then moves the blocks so that $R_i \le S_i$ for every block, by a strength-weighted equidistribution, a sequential linear programme for the minimax problem (solved with HiGHS) and a final rounding to the floors with a local search (labelled 1-opt in the figure). Every candidate layout is re-analysed with the support model, including the welding gaps at its own block positions.
7. **Checks (Section 3.8).** For the given and for the designed layout, the shear force, bending moment and stresses of the hull girder are compared with the permissible envelopes, and every block is checked for bearing pressure, load against its strength and position relative to the floors.
8. **Verdict and report (Chapter 4).** The hull and block checks are combined with the outcome of the spacing design (the designed spacing fits or does not fit, with the reason), and the results are written as a summary, CSV tables, plots and the result tabs of the graphical user interface.

The analysis assumes that the dock is dry, that the ship is supported only by the keel block line, and that the problem is two-dimensional in the longitudinal–vertical plane.

## 3.2 Coordinate system, sign conventions and units

### 3.2.1 Coordinates and kinematics

The longitudinal coordinate $x$ is measured from the aft perpendicular (AP), positive forward. The hull extends from $x_{\text{aft}}$ to $x_{\text{fwd}}$; these normally coincide with $0$ and $L_{pp}$, but either end may lie outside the perpendiculars. The vertical displacement of the hull girder is $v(x)$, taken **positive downward**, that is, into the blocks. The rotation is $\theta = \mathrm{d}v/\mathrm{d}x$.

The block at position $x_i$ has an initial gap $g_i$, defined as the vertical distance between the undeformed keel line and the top of the block. A positive gap means the block top lies below the keel, so the hull must settle by $g_i$ before contact is made. The compression of block $i$ is

$$ \delta_i = v_i - g_i , \qquad v_i = v(x_i) \tag{3.1} $$

A block is in contact only when $\delta_i > 0$.

### 3.2.2 Forces and stress resultants

The load intensity $q(x)$ (weight per unit length) is positive downward. The block reaction $R_i$ is positive upward on the hull, which corresponds to compression in the block. The shear force at a section $x$ is defined as the sum of upward forces minus the sum of downward forces acting on the part of the hull aft of $x$:

$$ V(x) = \sum_{x_i < x} R_i - \int_{x_{\text{aft}}}^{x} q(\xi)\,\mathrm{d}\xi - \sum_{x_j < x} P_j \tag{3.2} $$

where $P_j$ are concentrated (point) weights. The bending moment is the integral of the shear force:

$$ M(x) = \int_{x_{\text{aft}}}^{x} V(\xi)\,\mathrm{d}\xi \tag{3.3} $$

With these definitions a simply supported beam under a uniform downward load has a positive mid-span moment. Hence **sagging is positive** and hogging negative. Normal stresses are taken positive in tension, so a sagging moment produces compression in the deck and tension in the keel.

The conventions are collected in Table 3.1.

| Quantity | Symbol | Positive sense |
|---|---|---|
| Longitudinal coordinate | $x$ | forward, origin at AP |
| Vertical displacement | $v$ | downward |
| Rotation | $\theta = \mathrm{d}v/\mathrm{d}x$ | bow moves down relative to stern |
| Load intensity, point weight | $q$, $P$ | downward |
| Longitudinal moment of a weight | $M_j = W_j\,\text{LCG}_j$ | weight forward of AP |
| Block reaction | $R_i$ | upward on hull (block in compression) |
| Block gap | $g_i$ | block top below keel line |
| Block compression | $\delta_i = v_i - g_i$ | block compressed |
| Keel lift from welding | $z$ | upward, away from the blocks |
| Shear force | $V$ | net upward force aft of the section |
| Bending moment | $M$ | sagging |
| Normal and shear stress | $\sigma$, $\tau$ | tension; same sign as $V$ |

Table: Table 3.1 – Sign conventions used throughout the methodology.

### 3.2.3 Units

All internal calculations use SI units (N, m, Pa). Input and output use the units customary in docking practice, converted at the boundary of the program:

- weights and block strengths in tonnes, with $1\ \text{t} = 9.80665\ \text{kN}$;
- longitudinal moments in t·m and weight-curve intensities in t/m;
- stresses and pressures in MPa, layer moduli in MPa and block stiffnesses in kN/mm;
- block gaps and keel lifts in millimetres, joint angles in mrad and shrinkage strains in mm/m.

## 3.3 Longitudinal weight distribution

The load model has to reproduce the total weight and its longitudinal centre exactly, because the block reactions are in equilibrium with them, and it has to keep the local distribution, because the elastic hull transfers local weights mainly to nearby blocks. Two forms of input are supported, which can be combined: discrete weight items and a weight-distribution curve.

### 3.3.1 Weight items with LCG or moment

The ship's docking weight is described by a list of weight items. Each item $j$ has a weight $W_j$, a longitudinal position of its centre of gravity and an extent $[x_{a,j},\,x_{b,j}]$. Hull structure, machinery, superstructure, outfit, tank contents and similar groups are entered in this way. This representation follows the usual practice of building the longitudinal weight curve from grouped weights for still-water strength calculations [@hughes2010; @lewis1988].

Weight tables in practice often list the longitudinal moment of each group about the aft perpendicular instead of its LCG, because moments can be added directly. The program therefore accepts the position either as $\text{LCG}_j$ or as the moment $M_j$ [t·m]; in the latter case

$$ \text{LCG}_j = \frac{M_j}{W_j} \tag{3.4} $$

If both are given, they must agree. The check uses a combined relative and absolute tolerance,

$$ \left| \frac{M_j}{W_j} - \text{LCG}_j \right| \le \max\!\left( 10^{-6} \max\!\left( \left|\tfrac{M_j}{W_j}\right|, |\text{LCG}_j| \right),\ 10^{-6}\ \text{m} \right) \tag{3.5} $$

and an item that violates it is rejected with an error that names both values. Accepting redundant but inconsistent data silently would hide input errors that shift the centre of gravity and thus every reaction. An item without an extent, or with a single position $x$, is a point weight at its LCG. The weight $W_j$ must be non-zero; negative values are allowed for deductions.

### 3.3.2 Trapezoidal spreading

Over its extent, each distributed item is replaced by a linearly varying load intensity, that is, a trapezoid. With $L = x_b - x_a$, local coordinate $s = x - x_a$ and end intensities $w_a$ and $w_b$, the intensity is

$$ w(s) = w_a + (w_b - w_a)\,\frac{s}{L}, \qquad 0 \le s \le L \tag{3.6} $$

The two end intensities are fixed by requiring the trapezoid to reproduce the item's weight and its centre of gravity. The resultant gives

$$ W_j = \int_0^L w(s)\,\mathrm{d}s = \frac{L}{2}\,(w_a + w_b) \tag{3.7} $$

and the first moment about the aft end of the item, with $\bar{x} = \text{LCG}_j - x_a$, gives

$$ W_j\,\bar{x} = \int_0^L w(s)\,s\,\mathrm{d}s = \frac{L^2}{6}\,(w_a + 2w_b) \tag{3.8} $$

Solving Equations (3.7) and (3.8) for the end intensities yields

$$ w_a = \frac{2W_j}{L}\left(2 - \frac{3\bar{x}}{L}\right), \qquad w_b = \frac{2W_j}{L}\left(\frac{3\bar{x}}{L} - 1\right) \tag{3.9} $$

The linear system of Equations (3.7)–(3.8) in the unknowns $(w_a, w_b)$ has the determinant $\tfrac{L}{2}\cdot\tfrac{2L^2}{6} - \tfrac{L}{2}\cdot\tfrac{L^2}{6} = \tfrac{L^3}{12} > 0$, so the trapezoid with a prescribed weight and centre of gravity is unique. For an item with its LCG at mid-length, $\bar{x} = L/2$, Equation (3.9) gives $w_a = w_b = W_j/L$, a uniform load as expected.

### 3.3.3 Middle-third condition

A physically meaningful weight distribution must not contain negative intensities (for a positive item weight). Requiring $w_a \ge 0$ and $w_b \ge 0$ in Equation (3.9) gives

$$ \frac{L}{3} \le \bar{x} \le \frac{2L}{3} \tag{3.10} $$

so the LCG of each item must lie in the middle third of its extent. At the limits $\bar{x} = L/3$ and $\bar{x} = 2L/3$ the trapezoid degenerates into a triangle with its peak at the aft or forward end respectively. Figure 3.2 shows the three cases for an item of 100 t spread over 12 m.

If the condition is violated, the program rejects the item. A balancing distribution with a negative end would be mathematically valid but physically meaningless; the item must instead be split, or its extent changed.

![Figure 3.2: Trapezoidal spreading of a 100 t item over 12 m with the LCG at L/3, inside the middle third and at 2L/3. The grey band marks the admissible middle third.](figures/fig_trapezoid.png)

### 3.3.4 Weight-distribution curve as exact trapezoids

When a detailed weight distribution is available, for example the lightship weight curve of the loading manual, it is entered directly as a curve: a sequence of points $(x_k, q_k)$, $k = 0, \dots, n$, with $x_0 \le x_1 \le \dots \le x_n$ and $q_k \ge 0$ [t/m]. The curve defines the piecewise-linear intensity

$$ q_c(x) = q_k + (q_{k+1} - q_k)\,\frac{x - x_k}{x_{k+1} - x_k}, \qquad x_k \le x \le x_{k+1}, \qquad q_c(x) = 0 \ \text{outside}\ [x_0, x_n] \tag{3.11} $$

Two consecutive points with the same $x$ define a jump of the intensity, for example at the end of a superstructure or of a cargo hold. The points must be given from aft to forward, every intensity must be finite and non-negative, and the curve must carry a positive weight; otherwise it is rejected.

Instead of adding a second kind of load to the program, every segment $k$ of length $L_k = x_{k+1} - x_k > 0$ with $q_k + q_{k+1} > 0$ is converted into one trapezoidal weight item with the extent $[x_k, x_{k+1}]$, the weight and the LCG

$$ W_k = \frac{(q_k + q_{k+1})\,L_k}{2}, \qquad \text{LCG}_k = x_k + L_k\,\frac{q_k + 2\,q_{k+1}}{3\,(q_k + q_{k+1})} \tag{3.12} $$

Segments of zero length (jumps) and unloaded segments ($q_k = q_{k+1} = 0$) carry no weight and are skipped. The conversion has three properties, which together show that it is exact.

*The item satisfies the middle-third condition automatically.* With $\bar{x}_k = \text{LCG}_k - x_k$, Equation (3.12) gives

$$ \frac{\bar{x}_k}{L_k} - \frac{1}{3} = \frac{q_{k+1}}{3\,(q_k + q_{k+1})} \ge 0, \qquad \frac{2}{3} - \frac{\bar{x}_k}{L_k} = \frac{q_k}{3\,(q_k + q_{k+1})} \ge 0 \tag{3.13} $$

because $q_k, q_{k+1} \ge 0$, so Equation (3.10) holds and the item is always accepted.

*The trapezoid reproduces the segment exactly.* Substituting Equation (3.12) into Equation (3.9), with $2W_k/L_k = q_k + q_{k+1}$ and $3\bar{x}_k/L_k = (q_k + 2q_{k+1})/(q_k + q_{k+1})$, gives

$$ w_a = (q_k + q_{k+1}) \cdot 2 - (q_k + 2\,q_{k+1}) = q_k, \qquad w_b = (q_k + 2\,q_{k+1}) - (q_k + q_{k+1}) = q_{k+1} \tag{3.14} $$

A linear function on an interval is determined by its two end values, so the spread item coincides with $q_c$ on $[x_k, x_{k+1}]$ pointwise, not only in its weight and centre. Equivalently, since the trapezoid with a given weight and LCG is unique (Section 3.3.2) and the segment itself is such a trapezoid, the conversion cannot produce anything else.

*Weight and moment are preserved.* Summing over the loaded segments $\mathcal{K}$,

$$ W_c = \sum_{k \in \mathcal{K}} \frac{(q_k + q_{k+1})\,L_k}{2}, \qquad M_c = \sum_{k \in \mathcal{K}} W_k\,\text{LCG}_k = \sum_{k \in \mathcal{K}} \frac{L_k}{6}\Big[q_k\,(2x_k + x_{k+1}) + q_{k+1}\,(x_k + 2x_{k+1})\Big] \tag{3.15} $$

which are the exact integrals $\int q_c\,\mathrm{d}x$ and $\int q_c\,x\,\mathrm{d}x$ of Equation (3.11); the second form follows from $W_k\,\text{LCG}_k = W_k\,x_k + L_k^2\,(q_k + 2q_{k+1})/6$. The centre of the curve is $\text{LCG}_c = M_c / W_c$.

Because the segments become ordinary items, the mesh, the load model, the exact integration of Section 3.3.6 and all later steps apply unchanged; the curve needs no special treatment anywhere else in the program. Its segment ends become key points of the mesh (Section 3.3.5), so a jump is represented exactly as a discontinuity between two elements. Figure 3.3 shows the weight-curve example used in Chapter 6: a lightship curve of 18 points with five jumps, converted into 12 trapezoidal items of 4343.0 t in total at LCG 55.67 m, and three deadweight items entered with their moments (200 t with 14 000 t·m, 60 t with 6000 t·m, and a point weight of 40 t with 1200 t·m), which give 4643.0 t at LCG 56.639 m.

![Figure 3.3: Weight distribution of the weight-curve example: lightship weight curve with jumps (line and points), deadweight items entered with their longitudinal moments (dashed, the point weight as an arrow) and the total intensity (shaded).](figures/fig_weight_curve.png)

### 3.3.5 Superposition and mesh construction

The load intensities of all items, including the items obtained from the curve, are superposed. The resulting $q(x)$ is piecewise linear, with discontinuities wherever an item starts or ends. To represent it without approximation, the one-dimensional mesh used by all subsequent calculations is built from a set of *key points*:

- the hull ends;
- every item boundary $x_{a,j}$ and $x_{b,j}$;
- every boundary of a hull girder section;
- every block position $x_i$.

The key points are sorted, and points closer than $10^{-6}$ m are merged. Each interval of length $\Delta x$ between consecutive key points is then divided uniformly into $\lceil \Delta x / h_{\max} \rceil$ elements, where $h_{\max}$ is a user-defined maximum element length (0.5 m by default).

Because no element straddles an item boundary, the load on every element $e$ between nodes $x_1$ and $x_2$ is exactly linear. It is stored by its end values $q_1^{(e)}$ and $q_2^{(e)}$. Point weights are stored as nodal loads $P_n$ at the node coinciding with their position. Because the block positions are key points, the mesh depends on the layout; the spacing design of Section 3.12 therefore rebuilds the mesh and the load model for every candidate layout.

Figure 3.4 shows the distribution obtained for the example vessel used in Chapter 6. Its data are placeholder values chosen for illustration.

![Figure 3.4: Longitudinal weight distribution of the example vessel built from nine weight items (placeholder data). The rudder is a point weight.](figures/fig_weight_distribution.png)

### 3.3.6 Exact integration and weight checks

For a linear load on an element of length $L_e = x_2 - x_1$, the resultant is

$$ F^{(e)} = \frac{L_e}{2}\left(q_1^{(e)} + q_2^{(e)}\right) \tag{3.16} $$

The first moment about $x = 0$ is obtained by writing the load in terms of the linear interpolation functions $(x_2 - x)/L_e$ and $(x - x_1)/L_e$ and integrating:

$$ \int_{x_1}^{x_2} q(x)\,x\,\mathrm{d}x = \frac{L_e}{6}\Big[q_1^{(e)}\,(2x_1 + x_2) + q_2^{(e)}\,(x_1 + 2x_2)\Big] \tag{3.17} $$

The total weight and its LCG then follow from

$$ W = \sum_e F^{(e)} + \sum_n P_n, \qquad \text{LCG} = \frac{1}{W}\left(\sum_e \int_{x_1}^{x_2} q\,x\,\mathrm{d}x + \sum_n P_n x_n\right) \tag{3.18} $$

Because the integration is exact, $W$ and LCG from Equation (3.18) agree with the sums over the input items, and with Equation (3.15) for the curve, to machine precision. The program uses this agreement as an internal consistency check. A separate, user-controlled tolerance (0.5 % by default) compares the computed weight and LCG with the declared docking weight and LCG from the docking plan, and a warning is issued if they differ.

For the example vessel the mesh contains 256 nodes. The computed weight is 4500.0 t at LCG = 50.469 m, against a declared 4500.0 t at 50.47 m.

## 3.4 Block model and block plan

### 3.4.1 Stiffness of a layered block

A keel block typically consists of a concrete or steel base, one or more hardwood layers and a softwood cap that conforms to the hull. The layers carry the same force and their compressions add, so they act as springs in series. For layer $m$ with thickness $t_m$ and modulus $E_m$ in the loading direction (perpendicular to the grain for timber), loaded over the contact area $A$, the block stiffness is

$$ \frac{1}{k_i} = \sum_m \frac{t_m}{E_m\,A} \tag{3.19} $$

The contact area is the product of the transverse contact width $b$ (keel or flat-bottom width resting on the block) and the longitudinal block length $\ell$:

$$ A = b\,\ell \tag{3.20} $$

The softest layer, usually the softwood cap, dominates the sum in Equation (3.19), so the block stiffness is sensitive to the cap material and thickness. The program also accepts a measured or prescribed stiffness, which overrides Equation (3.19).

For the example keel block the layers are 0.90 m of concrete ($E$ = 30 GPa), 0.20 m of hardwood ($E$ = 1.0 GPa) and 0.10 m of softwood ($E$ = 0.4 GPa), with $b$ = 0.8 m and $\ell$ = 1.2 m. These give $A = 0.96\ \text{m}^2$ and $k = 2.0 \times 10^{9}$ N/m (2000 kN/mm). These values are placeholders, not recommended design data.

### 3.4.2 Compression-only bilinear law

A block can push on the hull but cannot pull it down. Beyond a certain stress, the timber cap crushes, that is, it deforms permanently with a much reduced tangent stiffness. Both effects are captured by the piecewise-linear law

$$ R_i(\delta_i) = \begin{cases} 0, & \delta_i \le 0 \quad \text{(open)} \\ k_i\,\delta_i, & 0 < \delta_i \le \delta_{y,i} \quad \text{(elastic)} \\ k_{2,i}\,\delta_i + (k_i - k_{2,i})\,\delta_{y,i}, & \delta_i > \delta_{y,i} \quad \text{(crushed)} \end{cases} \tag{3.21} $$

with the crush force, the compression at which crushing starts, and the post-crush stiffness

$$ F_{y,i} = \sigma_{c,i}\,A_i, \qquad \delta_{y,i} = \frac{F_{y,i}}{k_i}, \qquad k_{2,i} = \rho_i\,k_i \tag{3.22} $$

Here $\sigma_{c}$ is the crush stress of the block and $\rho \in [0,1]$ is the post-crush stiffness ratio (default 0.05). If no crush stress is specified, $\delta_{y,i} = \infty$ and the law reduces to a linear compression-only spring. The law is continuous and monotonically non-decreasing in $\delta_i$, a property used in Sections 3.7 and 3.11.

Figure 3.5 shows the law for the example block with $\sigma_c$ = 4.0 MPa. The crush force is then $F_y$ = 3.84 MN (391.6 t), reached at $\delta_y$ = 1.92 mm, and the post-crush stiffness is $k_2 = 1.0 \times 10^8$ N/m.

![Figure 3.5: Bilinear compression-only force law of the example keel block (placeholder properties).](figures/fig_block_spring_law.png)

### 3.4.3 State-wise linearisation

Within each branch of Equation (3.21) the force is an affine function of the compression. Introducing a block state $s_i \in \{\text{open}, \text{elastic}, \text{crushed}\}$, the law can be written uniformly as

$$ R_i = k_{t,i}(s_i)\,\delta_i + c_i(s_i) \tag{3.23} $$

with the tangent stiffness $k_{t,i}$ and offset $c_i$ listed in Table 3.2. For a fixed set of states, every support model in this thesis becomes linear in the unknowns. Section 3.7 describes how the correct states are found.

| State $s_i$ | Range of $\delta_i$ | $k_{t,i}$ | $c_i$ |
|---|---|---|---|
| open | $\delta_i \le 0$ | $0$ | $0$ |
| elastic | $0 < \delta_i \le \delta_{y,i}$ | $k_i$ | $0$ |
| crushed | $\delta_i > \delta_{y,i}$ | $k_{2,i}$ | $(k_i - k_{2,i})\,\delta_{y,i}$ |

Table: Table 3.2 – Linearisation of the block law for each contact state.

Substituting Equation (3.1) into Equation (3.23) expresses the reaction in terms of the hull displacement:

$$ R_i = k_{t,i}\,v_i - \left(k_{t,i}\,g_i - c_i\right) = k_{t,i}\,v_i - r_i, \qquad r_i = k_{t,i}\,g_i - c_i \tag{3.24} $$

The term $r_i$ collects the effect of the gap and the crushing offset. In both support models it appears as a known load.

### 3.4.4 Block plan: number, spacing, strength and gap of every block

The blocks are defined by a single *block plan*, which is the object of the spacing design. A plan names one block type (Section 3.4.1–3.4.2), the position $x_1$ of the first (aftmost) block, the number of blocks $N \ge 2$, and three per-block quantities, each given either as one value for all blocks or as a list:

- the spacings $s_1, \dots, s_{N-1} > 0$, where $s_i$ is the distance from block $i$ to block $i+1$;
- the strengths $S_1, \dots, S_N > 0$ [t], the force that each block may carry;
- the initial gaps $g_1^{0}, \dots, g_N^{0}$ [mm], for example from uneven block heights.

If the number is omitted it is inferred from the length of a list. The block positions follow by accumulating the spacings:

$$ x_1 = x_{\text{start}}, \qquad x_{i+1} = x_i + s_i, \qquad i = 1, \dots, N-1 \tag{3.25} $$

so the positions are strictly increasing by construction. The blocks are labelled $B01, B02, \ldots$ from aft to forward. The strength of block $i$ is

$$ S_i = \begin{cases} S_i^{\text{plan}}, & \text{if given in the plan} \\ R_{\text{cap}}, & \text{otherwise (capacity of the block type)} \end{cases} \tag{3.26} $$

The strength and the initial gap belong to the block, not to a position: when the spacing design moves block $i$, its strength $S_i$ and its gap $g_i^{0}$ move with it. The design keeps the number and the order of the blocks, so every block keeps its strength. This reflects the practical situation in the dock, where a known set of blocks, possibly of different condition or rating, is available and only their positions are to be chosen.

The strength is a limit, not a stiffness: it does not enter the support models. The reactions therefore do not depend on the strengths, which enter only the block checks (Section 3.8.5) and the objective of the spacing design (Section 3.12).

## 3.5 Method A: rigid hull on elastic blocks

### 3.5.1 Formulation

If the hull girder is much stiffer than the block line, its deflection can be neglected relative to the block compressions, and the hull moves as a rigid body. To improve the numerical conditioning, positions are measured from the centroid of the block positions, $x_0 = \frac{1}{N}\sum_i x_i$, so that $x_{c,i} = x_i - x_0$. The displacement is then linear:

$$ v(x) = a + b\,(x - x_0), \qquad v_i = a + b\,x_{c,i} \tag{3.27} $$

where $a$ is the settlement at $x_0$ and $b$ the rotation. The reactions follow from Equation (3.24) as $R_i = k_{t,i}(a + b\,x_{c,i}) - r_i$. Vertical force equilibrium and moment equilibrium about $x_0$ require

$$ \sum_i R_i = W, \qquad \sum_i R_i\,x_{c,i} = M_0 \tag{3.28} $$

where $M_0 = M_W - W x_0$ is the first moment of the load about $x_0$ and $M_W = W \cdot \text{LCG}$ is the first moment about the origin from Equation (3.18). Substitution gives the $2 \times 2$ linear system

$$ \begin{bmatrix} \sum k_{t,i} & \sum k_{t,i}\,x_{c,i} \\ \sum k_{t,i}\,x_{c,i} & \sum k_{t,i}\,x_{c,i}^2 \end{bmatrix} \begin{bmatrix} a \\ b \end{bmatrix} = \begin{bmatrix} W + \sum r_i \\ M_0 + \sum r_i\,x_{c,i} \end{bmatrix} \tag{3.29} $$

The coefficient matrix is symmetric. It is non-singular if and only if at least two blocks with $k_{t,i} > 0$ are in contact at distinct positions, which is the same support condition that applies to the elastic model (Section 3.6.6).

Centring does not change the solution, but it reduces the condition number of the matrix. Without centring, the off-diagonal terms grow with the distance of the block line from the origin, while with centring they vanish for a symmetric block arrangement of equal stiffness.

### 3.5.2 Closed form for equal stiffness

For $N$ blocks of equal stiffness $k$, no gaps and no crushing, $\sum x_{c,i} = 0$ and Equation (3.29) decouples. The solution is the familiar "P/A ± Mc/I" form:

$$ R_i = \frac{W}{N} + \frac{W\,e\,(x_i - \bar{x})}{\sum_m (x_m - \bar{x})^2}, \qquad e = \text{LCG} - \bar{x} \tag{3.30} $$

where $\bar{x} = x_0$ is the block centroid and $e$ the eccentricity of the weight relative to it. The reactions therefore vary linearly along the block line, increasing towards the side on which the centre of gravity lies. In particular, all reactions are equal to $W/N$ when the block centroid lies under the centre of gravity, a fact used in Section 3.12.7.

### 3.5.3 Continuous block line: trapezoidal and triangular distributions

When the blocks are closely and evenly spaced, the line can be idealised as a continuous elastic support of length $L_b$, and the reaction becomes a load per unit length $p(x)$. Applying Equation (3.30) in continuous form, the end values of the trapezoidal reaction distribution are

$$ p_{\text{end}} = \frac{W}{L_b}\left(1 \pm \frac{6e}{L_b}\right), \qquad |e| \le \frac{L_b}{6} \tag{3.31} $$

This is the classical trapezoidal method of docking calculations [@nstm997; @milstd1625]. If $|e| > L_b/6$, Equation (3.31) would give a negative end value, which a compression-only support cannot provide. Contact is then lost over part of the line, and the reaction becomes triangular. Its length $c$ and peak value follow from force and moment equilibrium of a triangle whose centroid lies at the LCG:

$$ c = 3\left(\frac{L_b}{2} - |e|\right), \qquad p_{\text{peak}} = \frac{2W}{c} \tag{3.32} $$

Equations (3.31) and (3.32) are not used directly by the program. They are the continuous limit of the discrete formulation, and they explain the lift-off behaviour that the contact iteration reproduces for discrete blocks.

### 3.5.4 Limitations

Method A satisfies equilibrium exactly but ignores the deformation of the hull. It therefore cannot represent:

- the concentration of reactions under heavy local weights;
- the additional load attracted by the end blocks from unsupported overhangs;
- the redistribution caused by gaps between the keel and individual blocks.

It does, however, provide an important reference. It is the limit that the elastic model must approach as the hull stiffness tends to infinity (Chapter 5), and it is the basis of established docking practice.

## 3.6 Method B: elastic hull girder FEM

### 3.6.1 Governing equations

The hull girder is idealised as a beam with bending stiffness $EI(x)$ and, optionally, effective shear stiffness $\kappa G A(x)$ [@hughes2010]. In Euler–Bernoulli theory [@gere2013], plane sections remain plane and normal to the deflected axis, and the deflection under the net load $q_{\text{net}}(x)$ satisfies

$$ \frac{\mathrm{d}^2}{\mathrm{d}x^2}\left(EI\,\frac{\mathrm{d}^2 v}{\mathrm{d}x^2}\right) = q_{\text{net}}(x) \tag{3.33} $$

where $q_{\text{net}}$ is the weight intensity minus the block reactions, the latter acting as concentrated forces.

The Timoshenko theory [@timoshenko1921] relaxes the normality assumption. It introduces the section rotation $\psi$ as an independent field, so that the shear strain is $\gamma = \mathrm{d}v/\mathrm{d}x - \psi$, and the equations become

$$ \frac{\mathrm{d}}{\mathrm{d}x}\left(EI\,\frac{\mathrm{d}\psi}{\mathrm{d}x}\right) + \kappa G A\left(\frac{\mathrm{d}v}{\mathrm{d}x} - \psi\right) = 0, \qquad \frac{\mathrm{d}}{\mathrm{d}x}\left[\kappa G A\left(\frac{\mathrm{d}v}{\mathrm{d}x} - \psi\right)\right] = -q_{\text{net}} \tag{3.34} $$

In the limit $\kappa G A \to \infty$ the shear strain vanishes, $\psi \to \mathrm{d}v/\mathrm{d}x$, and eliminating the shear force between the two equations recovers Equation (3.33). Shear deformation is significant for short, deep girders and for the short spans between closely spaced blocks. The relative importance of shear over an element of length $L$ is measured by the dimensionless parameter

$$ \Phi = \frac{12\,EI}{\kappa G A\,L^2} \tag{3.35} $$

For the example hull ($EI = 4.12 \times 10^{12}$ N·m², $\kappa G A = 7.14 \times 10^{10}$ N amidships), $\Phi$ is 0.048 for the full length of 120 m but 120 for a block spacing of 2.4 m. At the scale of individual block spans, shear deformation therefore dominates the local flexibility.

### 3.6.2 Finite element discretisation

The beam is divided into two-node elements on the mesh of Section 3.3.5, with two degrees of freedom per node, $v$ and $\theta$. For node $n$ the global degrees of freedom are numbered $2n$ (displacement) and $2n+1$ (rotation). Within an element of length $L$, with local coordinate $\xi = s/L \in [0,1]$, the Euler–Bernoulli displacement is interpolated by the Hermite cubic shape functions [@przemieniecki1968; @cook2002]:

$$ N_1 = 1 - 3\xi^2 + 2\xi^3, \quad N_2 = L\,(\xi - 2\xi^2 + \xi^3), \quad N_3 = 3\xi^2 - 2\xi^3, \quad N_4 = L\,(-\xi^2 + \xi^3) \tag{3.36} $$

so that $v(\xi) = N_1 v_1 + N_2 \theta_1 + N_3 v_2 + N_4 \theta_2$. The functions are plotted in Figure 3.6. Each takes the value one for its own degree of freedom and zero for the other three, in both displacement and slope.

![Figure 3.6: Hermite cubic shape functions of the two-node beam element; N₂ and N₄ are shown divided by the element length.](figures/fig_hermite_shapes.png)

### 3.6.3 Element stiffness matrix

For the Euler–Bernoulli element, the stiffness matrix follows from the strain energy $U = \frac12 \int_0^L EI\,(v'')^2\,\mathrm{d}s$:

$$ \mathbf{K}_e = \int_0^L EI\,\mathbf{B}^{\mathsf T}\mathbf{B}\,\mathrm{d}s, \qquad \mathbf{B} = \frac{\mathrm{d}^2}{\mathrm{d}s^2}\begin{bmatrix} N_1 & N_2 & N_3 & N_4 \end{bmatrix} \tag{3.37} $$

Evaluating Equation (3.37) for constant $EI$ gives the classical matrix. If the Timoshenko element is formulated with interdependent interpolation of $v$ and $\psi$, which reproduces the exact homogeneous solution of Equation (3.34) [@przemieniecki1968; @bathe2014], the stiffness matrix becomes

$$ \mathbf{K}_e = \frac{EI}{L^3\,(1+\Phi)} \begin{bmatrix} 12 & 6L & -12 & 6L \\ 6L & (4+\Phi)L^2 & -6L & (2-\Phi)L^2 \\ -12 & -6L & 12 & -6L \\ 6L & (2-\Phi)L^2 & -6L & (4+\Phi)L^2 \end{bmatrix} \tag{3.38} $$

with $\Phi$ from Equation (3.35). Setting $\Phi = 0$ recovers the Euler–Bernoulli matrix exactly. The program therefore uses Equation (3.38) for every element and sets $\Phi = 0$ wherever no shear stiffness is given, or where shear deformation is switched off by the user.

### 3.6.4 Consistent load vector

The distributed load is transferred to the nodes as work-equivalent (consistent) nodal forces and moments, $\mathbf{f}_e = \int_0^L \mathbf{N}^{\mathsf T} q(s)\,\mathrm{d}s$. For the linear load $q(\xi) = q_1 (1-\xi) + q_2\,\xi$, integration of the products of the Hermite functions with the two linear functions gives

$$ \mathbf{f}_e = \begin{bmatrix} \dfrac{L\,(7q_1 + 3q_2)}{20} & \dfrac{L^2\,(3q_1 + 2q_2)}{60} & \dfrac{L\,(3q_1 + 7q_2)}{20} & -\dfrac{L^2\,(2q_1 + 3q_2)}{60} \end{bmatrix}^{\mathsf T} \tag{3.39} $$

For a uniform load $q_1 = q_2 = q$ this reduces to the familiar fixed-end values $qL/2$ and $\pm qL^2/12$. Point weights $P_n$ are added directly to the displacement degree of freedom of their node.

### 3.6.5 Assembly and support springs

The global stiffness matrix $\mathbf{K}_0$ and load vector $\mathbf{F}_0$ are assembled by adding each element contribution to the rows and columns of its four degrees of freedom $(2e, 2e+1, 2e+2, 2e+3)$:

$$ \mathbf{K}_0 = \sum_e \mathbf{A}_e^{\mathsf T}\,\mathbf{K}_e\,\mathbf{A}_e, \qquad \mathbf{F}_0 = \sum_e \mathbf{A}_e^{\mathsf T}\,\mathbf{f}_e + \sum_n P_n\,\mathbf{e}_{2n} \tag{3.40} $$

where $\mathbf{A}_e$ is the Boolean gather matrix of element $e$ and $\mathbf{e}_{2n}$ is the unit vector of the displacement degree of freedom of node $n$. In the implementation the assembly is performed in coordinate (COO) sparse format, and duplicate entries are summed.

The blocks enter through Equation (3.24). Moving the stiffness term to the left-hand side, each block in state $s_i$, located at node $n_i$, contributes

$$ K_{2n_i,\,2n_i} \mathrel{+}= k_{t,i}, \qquad F_{2n_i} \mathrel{+}= r_i = k_{t,i}\,g_i - c_i \tag{3.41} $$

and the linear system for a given set of block states is

$$ \left(\mathbf{K}_0 + \mathbf{K}_s\right)\mathbf{u} = \mathbf{F}_0 + \mathbf{F}_s \tag{3.42} $$

where $\mathbf{u}$ contains all nodal displacements and rotations, $\mathbf{K}_s$ is the diagonal spring matrix and $\mathbf{F}_s$ the vector of spring loads. After solution, the block compressions follow from Equation (3.1) and the reactions from Equation (3.23).

### 3.6.6 Free–free support condition

No displacement boundary conditions are imposed: the ship is held only by the blocks. The unsupported stiffness matrix $\mathbf{K}_0$ is therefore singular. Its null space is spanned by the two rigid-body modes, vertical translation ($v = 1$, $\theta = 0$ at every node) and rotation ($v = x$, $\theta = 1$). The spring matrix $\mathbf{K}_s$ removes both modes only if at least two blocks with $k_{t,i} > 0$ are in contact at distinct longitudinal positions.

The program checks this condition before every solve. If it is not met, the support is reported as unstable instead of attempting to solve a singular system. The usual physical cause is a centre of gravity outside the block line, for which no compression-only support can provide moment equilibrium. A preliminary check rejects such cases before the analysis starts.

### 3.6.7 Solution of the linear system

The matrix in Equation (3.42) is symmetric, sparse and banded (each row has at most six non-zeros), and it is positive definite whenever the support condition of Section 3.6.6 holds. It is solved with a sparse direct method [@davis2006], using the SuperLU-based solver [@li2005superlu] provided by SciPy [@virtanen2020scipy]. For the meshes of interest, with a few hundred to a few thousand degrees of freedom, a solve takes milliseconds. This matters for the spacing design, which analyses hundreds of candidate layouts (Section 3.13).

### 3.6.8 Accuracy of the discretisation

For the Euler–Bernoulli beam with Hermite cubic elements and consistent loads, the nodal displacements are exact for any load that is polynomial within each element. This is because the element Green's functions lie in the space of cubic polynomials [@cook2002]. Since the mesh of Section 3.3.5 makes the load exactly linear within each element, and blocks act only at nodes, the computed reactions do not depend on the element size. Refining the mesh only increases the resolution of the plotted displacement, shear force and moment curves. The Timoshenko element of Equation (3.38) is likewise exact for uniform loads. Chapter 5 confirms both properties numerically. For the spacing design, this property has a further consequence: the reactions of a candidate layout do not contain a discretisation error that changes with the mesh, which would otherwise appear as noise in the finite-difference sensitivities of Section 3.12.4.

### 3.6.9 Beam on elastic foundation: characteristic length

Insight into the elastic model is obtained by smearing closely spaced blocks of stiffness $k$ at spacing $s$ into a continuous foundation of modulus $k_f = k/s$. The governing equation becomes that of a beam on an elastic foundation [@hetenyi1946]:

$$ EI\,\frac{\mathrm{d}^4 v}{\mathrm{d}x^4} + k_f\,v = q(x) \tag{3.43} $$

Its homogeneous solutions decay with the characteristic parameter

$$ \lambda = \left(\frac{k_f}{4\,EI}\right)^{1/4} \tag{3.44} $$

A concentrated load is effectively spread over a length of the order of $\pi/\lambda$ on either side of its point of application. If this length is large compared with the block line, the hull behaves almost rigidly and Method A is adequate. If it is short, reactions are dominated by the local weight distribution and overhangs.

For the example vessel, $k_f = 2.0 \times 10^9 / 2.4 = 8.33 \times 10^8$ N/m² and $EI = 4.12 \times 10^{12}$ N·m². These give $\lambda = 0.0843$ m⁻¹, $1/\lambda = 11.86$ m and $\pi/\lambda = 37.3$ m. This is substantially shorter than the 103.2 m block line, so significant differences between Methods A and B are to be expected (Section 6.3). Equation (3.44) also shows why the spacing is an effective design variable for the elastic hull: a local reduction of the spacing raises the local foundation modulus $k_f$ and draws load to that region, while a larger spacing sheds load to the neighbours.

## 3.7 Contact and crushing iteration

### 3.7.1 Complementarity formulation

Without crushing, the compression-only block law of Equation (3.21) is a set of complementarity conditions [@cottle1992; @wriggers2006]. For each block,

$$ R_i \ge 0, \qquad k_i\,\delta_i - R_i \le 0, \qquad R_i\,\left(k_i\,\delta_i - R_i\right) = 0 \tag{3.45} $$

In words: the reaction is non-negative, it is not smaller than the elastic spring force $k_i\,\delta_i$ (which is negative for an open gap), and at least one of the two inequalities holds with equality. Either the block is open with $R_i = 0$ and $\delta_i \le 0$, or it is in contact with $R_i = k_i\,\delta_i > 0$. Together with the linear equilibrium equations (3.29) or (3.42), this forms a linear complementarity problem. Crushing adds a third linear branch, so the problem becomes a piecewise-linear system with three possible states per block.

### 3.7.2 State iteration

Because each branch is linear, the problem is solved by iterating on the block states. This is an active-set strategy. Each iteration solves one linear system (Equation (3.29) or (3.42)) for the assumed states and updates the states from the computed compressions. The algorithm used for both support models is given as Algorithm 3.1.

**Algorithm 3.1 – Block-state iteration for compression-only, crushing supports.**

```text
input : block data (x_i, k_i, g_i, δy_i, k2_i), linear solver S(kt, c) -> δ
output: compressions δ_i, states s_i, reactions R_i

 1  s_i <- elastic for all blocks                  (assume full contact)
 2  seen <- {} ; one_at_a_time <- false
 3  for iteration = 1, 2, ..., max_iterations
 4      (kt_i, c_i) <- linearisation of s_i       (Table 3.2)
 5      check that >= 2 blocks with kt_i > 0 at distinct x   (else: unstable)
 6      δ <- S(kt, c)                              (solve Eq. 3.29 or 3.42)
 7      s'_i <- state implied by δ_i              (Eq. 3.21 ranges)
 8      if s' = s then return δ, s, R_i = kt_i δ_i + c_i
 9      add s to seen
10      if one_at_a_time or s' in seen or s' leaves < 2 supports then
11          one_at_a_time <- true
12          j <- block with largest violation η_j among changed blocks
13          s'' <- s ; s''_j <- s'_j ; s' <- s''
14      s <- s'
15  error: no convergence within max_iterations
```

### 3.7.3 Cycle detection and fallback

In the standard mode (line 14), all inconsistent states are updated simultaneously. This usually converges in very few iterations: for the example vessel a single solve suffices with the initial all-elastic assumption, and for the welded hull of Chapter 6 three (elastic hull) and five (rigid hull) iterations were needed. For some configurations, however, simultaneous updates can cycle between sets of states. To prevent such cycling, every visited state set is recorded. If a state set is about to be revisited — or if the simultaneous update would leave fewer than two separated blocks in contact, which can happen with large initial gaps — the algorithm switches permanently to a one-at-a-time mode, in which only the block with the largest violation of its assumed state is changed per iteration. The violation measure is

$$ \eta_i = \begin{cases} \max(\delta_i, 0), & s_i = \text{open} \\ \max(-\delta_i, 0) + \max(\delta_i - \delta_{y,i}, 0), & s_i = \text{elastic} \\ \max(\delta_{y,i} - \delta_i, 0), & s_i = \text{crushed} \end{cases} \tag{3.46} $$

that is, the distance of the computed compression from the admissible range of the assumed state. Changing a single block at a time corresponds to a classical pivoting step and avoids the oscillation of simultaneous updates. Because the force law is monotone, this procedure terminates in practice. A maximum iteration count (200 by default) is nevertheless retained as a safeguard, and Section 3.11 provides a robust fallback when it is reached.

### 3.7.4 Support stability check

At line 5 of Algorithm 3.1 the algorithm verifies the free–free support condition of Section 3.6.6 for the current state set: at least two blocks with positive tangent stiffness at distinct positions. Crushed blocks with $\rho = 0$ carry a constant force but no stiffness, so they do not count towards this condition. If it fails, the iteration stops with an instability error. Before any iteration, the analysis also verifies that the LCG lies within the extent of the block line. This is a necessary condition for any compression-only equilibrium.

### 3.7.5 Equilibrium residuals

After convergence, the global equilibrium residuals

$$ \varepsilon_F = \frac{\sum_i R_i - W}{W}, \qquad \varepsilon_M = \frac{\sum_i R_i\,x_i - W\cdot\text{LCG}}{W\,(x_{\text{fwd}} - x_{\text{aft}})} \tag{3.47} $$

are evaluated and reported. For the example vessel they were $1.7 \times 10^{-16}$ and $0$ for Method A and $0$ and $5.1 \times 10^{-15}$ for Method B, that is, at the level of round-off. Force equilibrium, $\sum_i R_i = W$, is also the basis of the lower bound on the achievable utilisation in Section 3.12.6.

## 3.8 Hull girder loads and checks; block strength and utilisation

### 3.8.1 Exact piecewise integration

Once the reactions are known, the shear force and bending moment follow from Equations (3.2) and (3.3). Because the load is linear within each element and all concentrated forces (point weights and block reactions) act at nodes, both curves can be integrated exactly.

Let $V_j^{+}$ and $M_j$ denote the shear force just forward of node $j$ and the bending moment at node $j$. Within element $j$, at distance $s \in [0, L]$ from its aft node and with $\Delta q = q_2 - q_1$,

$$ V(s) = V_j^{+} - \left(q_1\,s + \frac{\Delta q\,s^2}{2L}\right) \tag{3.48} $$

$$ M(s) = M_j + V_j^{+}\,s - \left(\frac{q_1\,s^2}{2} + \frac{\Delta q\,s^3}{6L}\right) \tag{3.49} $$

Thus $V$ is quadratic and $M$ cubic within each element. At node $j+1$ the shear force jumps by the net concentrated force applied there, while the bending moment is continuous:

$$ V_{j+1}^{+} = V(L) + \sum_{i \in j+1} R_i - P_{j+1}, \qquad M_{j+1} = M(L) \tag{3.50} $$

Starting from $V = M = 0$ at the aft end of the hull, Equations (3.48)–(3.50) are applied element by element. Several interior sample points per element are retained for plotting and for locating maxima. Each node is stored twice, with the shear values on either side of the jump, so that the stepped shape of the shear force curve is represented exactly.

### 3.8.2 Closure check

At the forward end of the hull, both $V$ and $M$ must return to zero, because the reactions are in equilibrium with the weight. The normalised closure errors

$$ \varepsilon_V = \frac{V(x_{\text{fwd}})}{\max |V|}, \qquad \varepsilon_M = \frac{M(x_{\text{fwd}})}{\max |M|} \tag{3.51} $$

provide an independent check on the combined weight, reaction and integration procedures. Values above $10^{-6}$ are reported as warnings. For the example vessel the largest of them was $1.6 \times 10^{-12}$.

### 3.8.3 Hull girder stresses

With section moduli $Z_{\text{deck}}$ and $Z_{\text{keel}}$ and the tension-positive convention, the longitudinal bending stresses at deck and keel are

$$ \sigma_{\text{deck}} = -\frac{M}{Z_{\text{deck}}}, \qquad \sigma_{\text{keel}} = \frac{M}{Z_{\text{keel}}} \tag{3.52} $$

so a sagging moment compresses the deck and stretches the keel. The mean shear stress over the effective shear area $A_s$ is

$$ \tau = \frac{V}{A_s} \tag{3.53} $$

Equation (3.53) is a simplification of the shear flow distribution in thin-walled sections [@hughes2010]. It is adequate for a screening check of the docked condition. Section properties are taken piecewise constant along the hull, as defined in the input.

### 3.8.4 Permissible envelopes and utilisation

The loading manual of a ship defines permissible still-water bending moments, with separate sagging and hogging limits, and permissible shear forces as functions of $x$ [@iacss11]. These are interpolated linearly to every sample point. The utilisation of the bending moment selects the limit corresponding to the sign of the moment:

$$ u_M(x) = \begin{cases} M(x) / M_{\text{sag}}(x), & M(x) \ge 0 \\ -M(x) / M_{\text{hog}}(x), & M(x) < 0 \end{cases}, \qquad u_V(x) = \frac{|V(x)|}{V_{\text{perm}}(x)} \tag{3.54} $$

The stresses of Equations (3.52) and (3.53) are compared in the same way with allowable bending and shear stresses. For every check the program reports the *governing* value: the location of maximum utilisation, together with the demand and limit there. A check passes when $u \le 1$. If no limit is defined for a quantity, its maximum absolute value is reported for information only. The hull girder checks are not constraints of the spacing design; they are evaluated for the given and for the designed layout, so that the effect of the new spacing on the hull girder is reported (Chapter 6).

### 3.8.5 Block strength and utilisation

The central quantity of this thesis is the utilisation of a block with respect to its own strength $S_i$ from Equation (3.26),

$$ u_i = \frac{R_i}{S_i}, \qquad U = \max_{i = 1, \dots, N} u_i \tag{3.55} $$

A layout is acceptable when $U \le u_{\text{target}}$, with a target utilisation $u_{\text{target}}$ that is normally 1; a smaller target provides an additional margin. $U$ is the objective of the spacing design (Section 3.12), and the same $u_i$ is reported as the block check "block load vs capacity". Two further block checks are evaluated. The bearing pressure on the contact area and its utilisation against the allowable pressure of the block type are

$$ p_i = \frac{R_i}{A_i}, \qquad u_{p,i} = \frac{p_i}{p_{\text{allow},i}} \tag{3.56} $$

and blocks in the crushed state of Equation (3.21) are reported separately, because they indicate permanent deformation even where the strength check passes. For every category the governing block, with the largest utilisation, is reported.

### 3.8.6 Position relative to the floors

Blocks should be located under transverse structure (floors, bulkheads, web frames), so that the reaction is carried into the hull without local overload of the bottom plating. Given the positions $x_k^{s}$ of such structure, the offset of block $i$ is

$$ d_i = \min_k \left|x_i - x_k^{s}\right| \tag{3.57} $$

A warning is issued if $d_i$ exceeds a user-defined tolerance (0.15 m by default). This is a positional check only; the local strength of the bottom structure is outside the scope of the present model. In the spacing design, the same structure positions serve as the admissible final block positions (Section 3.12.5), so that a designed layout satisfies $d_i = 0$.

## 3.9 Welding deformation of the assembled hull

### 3.9.1 Origin of the initial keel shape

In block construction the hull is assembled from prefabricated construction blocks that are erected on the building berth or in the dock and joined by erection butt welds [@storch1995; @okumoto2009]. The weld metal and the heat-affected zone shrink on cooling. Two effects of this shrinkage are important for the global keel shape [@ueda2012]:

- **angular distortion at the erection joints**, which makes the keel line kink at each joint;
- **differential longitudinal shrinkage** within a construction block, where the deck and upper structure, with their larger weld volume, shorten more than the bottom. This bends the block as a whole.

Both effects commonly make the keel concave upwards, so that the ends of the completed hull lift away from the block line. The first and last blocks can then lose contact, and their load is transferred to the remaining blocks. The model below describes this initial, stress-free keel shape kinematically. It does not attempt to compute welding distortion from the welding process; the distortion parameters are inputs, obtained from experience, from welding simulation [@ueda2012] or from a survey of the hull (Section 3.10). The implementation is in `src/shipsupport/deformation.py`.

### 3.9.2 Kinematic model

Let $z(x)$ denote the upward lift of the keel relative to a straight line. An angular distortion $\varphi_j$ at an erection joint at $x_j$ is represented by the symmetric kink function

$$ z_j(x) = \frac{\varphi_j}{2}\,\left| x - x_j \right| \tag{3.58} $$

whose slope jumps by exactly the joint angle:

$$ z_j'\!\left(x_j^{+}\right) - z_j'\!\left(x_j^{-}\right) = \varphi_j \tag{3.59} $$

The symmetric form lifts both hull ends equally; any other choice differs from it only by a straight line, which, as shown in Section 3.9.4, has no influence on the reactions.

A uniform shrinkage curvature $\kappa$ over a construction block $[x_a, x_b]$ is represented by $\kappa\,\psi(x)$, where $\psi$ has unit curvature inside the block and is straight outside it:

$$ \psi(x) = \begin{cases} 0, & x < x_a \\ \tfrac{1}{2}\,(x - x_a)^2, & x_a \le x \le x_b \\ (x_b - x_a)\left(x - \tfrac{1}{2}(x_a + x_b)\right), & x > x_b \end{cases} \tag{3.60} $$

$\psi$ and its slope are continuous at $x_a$ and $x_b$. If the deck and bottom shrinkage strains of the block are $\varepsilon_{\text{deck}}$ and $\varepsilon_{\text{bottom}}$ (shortening positive) and $D$ is the depth between them, plane sections give the curvature

$$ \kappa = \frac{\varepsilon_{\text{deck}} - \varepsilon_{\text{bottom}}}{D} \tag{3.61} $$

A positive $\varphi_j$ or $\kappa$ makes the keel concave upwards and lifts the hull ends. In the welding example of Chapter 6, a deck shrinkage of 0.24 mm/m, a bottom shrinkage of 0.20 mm/m and $D$ = 10 m give $\kappa = 4.0 \times 10^{-6}$ m⁻¹ (placeholder values).

When a survey of the keel line is available, its measured lifts $z_m$ at stations $x_m$ are interpolated to give a third contribution $z_s(x)$. Three interpolants are offered: shape-preserving piecewise cubic Hermite interpolation (PCHIP), which avoids overshoot between stations [@fritsch1980]; piecewise linear interpolation; and a natural cubic spline. All three pass exactly through the survey points and are extrapolated beyond the outermost stations. A survey measures the complete keel shape. When it is used directly, it therefore replaces the model terms: joints or zones given in addition are not added to it, and the program issues a warning. When it is used to identify the joint angles (Section 3.10), the zone curvatures are treated as known and only the joint angles are fitted. Formally, the keel lift is

$$ z(x) = \sum_j \frac{\varphi_j}{2}\,\left| x - x_j \right| + \sum_z \kappa_z\,\psi_z(x) \quad \text{(model)}, \qquad z(x) = z_s(x) \quad \text{(survey used directly)} \tag{3.62} $$

### 3.9.3 Conversion to block gaps

A lift of the keel above the block tops increases the gap of the block below it. Because the absolute datum and inclination of $z$ are arbitrary, the straight line that best fits the lift at the block positions is first removed:

$$ (a^{*}, b^{*}) = \arg\min_{a,\,b} \sum_{i=1}^{N} \left( z(x_i) - a - b\,x_i \right)^2 \tag{3.63} $$

The welding gap of block $i$ is then the remaining lift, shifted so that the lowest block point has zero gap:

$$ g_i^{w} = \tilde{z}_i - \min_k \tilde{z}_k, \qquad \tilde{z}_i = z(x_i) - a^{*} - b^{*} x_i \tag{3.64} $$

It is added to the initial gap $g_i^{0}$ of the block plan (Section 3.4.4), so the effective gap used by both support models is

$$ g_i = g_i^{0} + g_i^{w} \tag{3.65} $$

For the welding example, Equations (3.58)–(3.64) give gaps of 41.4 mm at the aft end block and 39.0 mm at the forward end block. Equations (3.63)–(3.64) depend on the block positions; Section 3.12.2 uses this dependence to make the welding deformation a constraint of the spacing design.

### 3.9.4 Invariance to a straight-line part

The normalisation of Equations (3.63)–(3.64) is admissible because a change of all gaps by a straight line, $\Delta g_i = \alpha + \beta\,x_i$, does not change any reaction or block state in either support model.

In the rigid-hull model (Section 3.5) the block compressions are $\delta_i = a + b\,x_{c,i} - g_i$. Replacing $g_i$ by $g_i + \alpha + \beta x_i$ and $(a, b)$ by $(a + \alpha + \beta x_0,\, b + \beta)$ leaves every $\delta_i$ unchanged. The equilibrium equations (3.28) contain only the reactions, so the shifted pair is the solution of the shifted problem.

In the elastic model (Section 3.6) the same argument applies to the rigid-body modes of the free–free beam. The displacement field $v(x) + \alpha + \beta x$, $\theta(x) + \beta$ produces no strain, so $\mathbf{K}_0$ of Equation (3.42) annihilates it, and the consistent load vector is unchanged:

$$ \mathbf{K}_0\,\mathbf{u}_{\text{rb}} = \mathbf{0}, \qquad \mathbf{u}_{\text{rb}} = \left[\ \alpha + \beta x_n,\ \beta\ \right]_{n} \tag{3.66} $$

Adding $\mathbf{u}_{\text{rb}}$ to the solution therefore satisfies the shifted equations with unchanged compressions and reactions. The automated tests of Chapter 5 confirm this numerically: the reactions with and without an added tilt and offset agreed to a relative tolerance of $10^{-8}$.

## 3.10 Identification of joint angles from a keel survey

A raw survey contains measurement noise and an unknown datum. When the erection joint positions are known, the survey can instead be used to estimate the joint angles, which gives a smooth, physically interpretable keel shape. The shrinkage curvature of the construction blocks is usually estimated separately, from the welding procedure, and is contained in the survey. Its contribution is therefore subtracted from the measured lifts before the fit; otherwise it would be counted twice, once in the identified angles and once in the zone term of Equation (3.62). The corrected survey lifts are modelled as

$$ z_m - \sum_z \kappa_z\,\psi_z(x_m) \approx c_0 + c_1\,\frac{x_m}{L_s} + \sum_{j=1}^{J} \varphi_j\,\frac{\left| x_m - x_j \right|}{2}, \qquad m = 1, \dots, p \tag{3.67} $$

where $c_0$ and $c_1$ absorb the survey datum and tilt and $L_s$ is the length of the surveyed range, which scales the second column. In matrix form $\mathbf{z} \approx \mathbf{A}\,\mathbf{p}$, where $\mathbf{z}$ now denotes the corrected lifts, with $\mathbf{p} = [c_0, c_1, \varphi_1, \dots, \varphi_J]^{\mathsf T}$. The parameters minimise the regularised least-squares functional [@hansen1998]

$$ \min_{\mathbf{p}}\ \left\| \mathbf{A}\,\mathbf{p} - \mathbf{z} \right\|_2^2 + \lambda_R\,L_s^2 \sum_{j=1}^{J} \varphi_j^2 \tag{3.68} $$

with Tikhonov weight $\lambda_R \ge 0$. Only the joint angles are penalised, not the datum terms. The minimiser satisfies the normal equations

$$ \left( \mathbf{A}^{\mathsf T}\mathbf{A} + \lambda_R\,L_s^2\,\mathbf{D} \right) \mathbf{p} = \mathbf{A}^{\mathsf T}\mathbf{z}, \qquad \mathbf{D} = \operatorname{diag}(0, 0, 1, \dots, 1) \tag{3.69} $$

In the implementation the problem is solved as an augmented least-squares problem by an orthogonal factorisation, which avoids forming $\mathbf{A}^{\mathsf T}\mathbf{A}$.

Between two consecutive joints each kink column of $\mathbf{A}$ is linear in $x$, so it is locally indistinguishable from the datum columns. A joint angle is therefore identifiable only through the change of slope across the joint, and at least one station is needed in every segment between joints. Sparse stations and survey noise make the problem ill-conditioned. For the welding example with stations every 5 m and Gaussian noise of 1 mm standard deviation, the condition number of $\mathbf{A}$ was about $1.2 \times 10^{3}$. The fit residual was 0.76 mm rms and the largest error of an identified angle was 0.25 mrad; halving the station spacing to 2.5 m reduced the largest angle error to 0.12 mrad. For the sample survey supplied with the program (Section 6.7: 25 stations, 0.5 mm noise, known zone shrinkage), the angles were recovered within 0.07 mrad with a residual of 0.295 mm rms. The regularisation term stabilises the estimate when stations are sparse or the joint spacing is small, at the cost of a bias towards zero angles. The identified angles define the keel lift $z(x)$ by Equation (3.62), which can then be evaluated at any block position, including the candidate positions of the spacing design.

## 3.11 Robust contact solution for large gaps

### 3.11.1 Failure of the direct state iteration

The state iteration of Algorithm 3.1 starts with all blocks elastic. For gaps of the order of the elastic settlement this start is close to the solution. Welding gaps are often an order of magnitude larger. The first linear solution then predicts tension in most blocks, and a simultaneous update of all states either leaves fewer than two supporting blocks or starts a cycle. The fallback of Section 3.7.3 then changes one state per iteration; on its own, however, it does not guarantee convergence within the iteration limit. For the example vessel with joint angles of 0.2 mrad only (end lifts of 13.6 and 13.7 mm), the direct iteration of the rigid-hull model did not converge in 200 iterations.

### 3.11.2 Gap continuation

The solution of the contact problem depends continuously on the gaps. For linear or bilinear monotone block laws it is the unique minimiser of the convex potential energy

$$ \Pi(\mathbf{u}) = \tfrac{1}{2}\,\mathbf{u}^{\mathsf T}\mathbf{K}_0\,\mathbf{u} - \mathbf{F}_0^{\mathsf T}\mathbf{u} + \sum_{i=1}^{N} \int_0^{\max(0,\,\delta_i(\mathbf{u}))} R_i(\delta)\,\mathrm{d}\delta \tag{3.70} $$

provided that at least two separated blocks are in contact [@cottle1992; @wriggers2006]. Small changes of the gaps therefore change only a few block states. This motivates a continuation (homotopy) in the gaps:

$$ \mathbf{g}(\lambda) = g_{\min}\,\mathbf{1} + \lambda \left( \mathbf{g} - g_{\min}\,\mathbf{1} \right), \qquad 0 \le \lambda \le 1 \tag{3.71} $$

At $\lambda = 0$ all gaps are equal to the smallest gap, which is a rigid shift of the block line; every block is in contact and the direct iteration converges as in the ordinary docking problem. The continuation parameter is then increased in steps, each started from the converged states of the previous step. A failed step is halved, a successful step is doubled up to 0.5. The procedure is given as Algorithm 3.2.

**Algorithm 3.2 – Gap continuation for the block-state iteration.**

```text
input : springs (k, crush law), gaps g, linear solver, iteration limit n_max
output: compressions δ, states S, total iterations
 1  try Algorithm 3.1 with gaps g, all states elastic           -> return if converged
 2  g_min <- min(g)
 3  (δ, S) <- Algorithm 3.1 with gaps g_min·1, all states elastic
 4  λ <- 0 ; Δλ <- 0.25
 5  while λ < 1:
 6      λ' <- min(1, λ + Δλ)
 7      try (δ', S') <- Algorithm 3.1 with gaps g(λ') of Eq. (3.71), starting from S
 8      if failed:   Δλ <- Δλ / 2 ;  stop with error if Δλ < 1/4096
 9      else:        λ <- λ' ; S <- S' ; Δλ <- min(2 Δλ, 0.5)
10  return δ', S, total iterations
```

For the 0.2 mrad case of Section 3.11.1, Algorithm 3.2 converged with one solve at $\lambda = 0$ and three continuation steps (11 state iterations in total). The final state had 28 open blocks, as expected for a rigid hull on a concave keel line. The continuation is used only when the direct iteration fails, so the results of ordinary cases are unaffected. For the spacing design the robustness of the contact solution is essential: every candidate layout on a welded hull has its own, large and uneven gaps, and a single failed analysis would stop the design.

## 3.12 Automatic spacing design from block strengths

The methods of Sections 3.3–3.11 analyse a given layout. This section turns them into a design method. The blocks are given by the block plan: their number $N$, their order and the strength $S_i$ of every block are fixed. The method computes the positions $x_1 < \dots < x_N$, and thus the spacing $x_{i+1} - x_i$ of every block, such that no block reaction exceeds the strength of its block. This spacing is the final result of the program. The implementation is in `src/shipsupport/spacing.py` (Section 4.10). The problem belongs to the class of support-position optimisation problems of structural optimisation [@mroz1975; @wang2004; @haftka1992], here with unilateral, crushable supports, a rigid or elastic hull, an individual strength for every support and an initial deformation of the supported structure.

### 3.12.1 Problem statement

Let $\mathbf{x} = (x_1, \dots, x_N)$ denote the block positions, and let $R_i(\mathbf{x})$ be the reaction of block $i$ computed by one of the support models (Method A or B, Sections 3.5–3.7 and 3.11) for the layout $\mathbf{x}$, with the gaps of Section 3.12.2. With the utilisations of Equation (3.55), the objective is

$$ U(\mathbf{x}) = \max_{i = 1, \dots, N} u_i(\mathbf{x}), \qquad u_i(\mathbf{x}) = \frac{\max\!\big(R_i(\mathbf{x}),\,0\big)}{S_i} \tag{3.72} $$

and the spacing design problem is the minimax problem

$$ \begin{aligned} \min_{\mathbf{x} \in \mathbb{R}^N} \quad & U(\mathbf{x}) \\ \text{s.t.} \quad & s_{\min} \le x_{i+1} - x_i \le s_{\max}, \qquad i = 1, \dots, N-1, \\ & x_{\text{lo}} \le x_1, \qquad x_N \le x_{\text{hi}}, \\ & x_i \in \mathcal{G} \quad \text{(admissible final positions, if given)}. \end{aligned} \tag{3.73} $$

The layout *fits* when its value satisfies $U \le u_{\text{target}}$. The constraints have the following meaning and defaults:

- $s_{\min}$ (default 1 m; in the examples the floor spacing of 0.8 m) is the smallest practicable distance between two blocks.
- $s_{\max}$ is the largest admissible span. By default it is twice the mean spacing of the given plan, $s_{\max} = 2\,(x_N^{0} - x_1^{0})/(N-1)$. A maximum spacing is necessary: without it, the minimax optimum leaves very long unsupported spans amidships, where the reactions are small, which is unacceptable in practice and would load the hull girder locally.
- $[x_{\text{lo}}, x_{\text{hi}}]$ is the allowed extent of the block line, by default the extent of the given plan, limited to the hull. With the option *fix end blocks*, $x_1$ and $x_N$ stay at their given positions and only the interior blocks move.
- $\mathcal{G}$ is the set of structure (floor) positions, or a regular grid, at which blocks may stand (Section 3.12.5).

The order constraint $x_1 < \dots < x_N$ is implied by $s_{\min} > 0$. The problem has feasible continuous solutions if and only if

$$ (N-1)\,s_{\min} \le x_{\text{hi}} - x_{\text{lo}} \quad \text{and, with fixed ends,} \quad x_{\text{hi}} - x_{\text{lo}} \le (N-1)\,s_{\max} \tag{3.74} $$

and the program rejects the problem with an explanatory error otherwise. With free ends, if $(N-1)\,s_{\max} < x_{\text{hi}} - x_{\text{lo}}$, the blocks cannot span the whole extent, and the extent is shortened to $x_{\text{hi}} = x_{\text{lo}} + (N-1)\,s_{\max}$.

Minimising $U$, rather than merely finding some layout with $U \le u_{\text{target}}$, gives the most even safety margin that the given blocks allow, and it gives a meaningful answer also when the plan does not fit: the smallest achievable utilisation, which quantifies how far the plan is from being feasible.

Problem (3.73) is difficult for three reasons. It is *non-smooth*: $U$ is the maximum of $N$ functions and is not differentiable where the maximising block changes; the reactions are only piecewise smooth, because block states change with the positions (Section 3.7) and the welding gaps have kinks (Section 3.12.2). It is *non-convex*. And every evaluation of $U$ is a complete nonlinear support analysis, with no analytical derivatives. The method therefore treats the analysis as a black box, `evaluate(blocks)`, which is passed to the design procedure (dependency injection, Section 4.10), and solves Problem (3.73) in three phases: a fast global heuristic (Section 3.12.3), a local minimax refinement (Section 3.12.4) and a discrete phase for the admissible positions (Section 3.12.5). Throughout, layouts are compared by the lexicographic score

$$ \operatorname{score}(\mathbf{x}) = \Big( \operatorname{round}\big(U(\mathbf{x}),\,10^{-9}\big),\ \sum_{i=1}^{N} u_i(\mathbf{x})^{8} \Big) \tag{3.75} $$

The first component is the objective. The second, a smooth measure dominated by the highest utilisations, breaks ties between layouts with practically the same maximum, so that a step that lowers the second-highest utilisations without changing the maximum still counts as progress. This matters for minimax problems, where several blocks are typically at the maximum at the optimum.

### 3.12.2 Welding deformation as a constraint

The welding deformation of the assembled hull (Section 3.9) is not a property of a layout but of the hull: the keel lift $z(x)$ is fixed, and the gap that a block sees depends on where it stands. Every candidate layout $\mathbf{x}$ is therefore analysed with the gaps recomputed at its own block positions, by Equations (3.63)–(3.65):

$$ g_i(\mathbf{x}) = g_i^{0} + g_i^{w}(\mathbf{x}), \qquad g_i^{w}(\mathbf{x}) = \tilde z_i(\mathbf{x}) - \min_k \tilde z_k(\mathbf{x}), \qquad \tilde z_i(\mathbf{x}) = z(x_i) - a^{*}(\mathbf{x}) - b^{*}(\mathbf{x})\,x_i \tag{3.76} $$

where $a^{*}(\mathbf{x})$ and $b^{*}(\mathbf{x})$ are the least-squares line of Equation (3.63) through the candidate positions and $g_i^{0}$ is the plan gap that moves with block $i$. The reactions of a candidate are thus

$$ R_i(\mathbf{x}) = \mathcal{R}_i\big(\mathbf{x},\ \mathbf{g}(\mathbf{x})\big) \tag{3.77} $$

where $\mathcal{R}$ denotes the solution of the contact problem (Algorithms 3.1 and 3.2) for given positions and gaps. In this sense the welding deformation is a constraint of the design: a layout is acceptable only if it satisfies $U \le u_{\text{target}}$ on the actual, deformed hull. Optimising the spacing for a straight keel and checking the welded hull afterwards would ignore that a block moved by the design also changes its gap.

The least-squares line in Equation (3.76) couples all blocks, but by the invariance of Section 3.9.4 it does not affect the reactions: $\mathbf{g}(\mathbf{x})$ and $\mathbf{g}^{0} + \big(z(x_1), \dots, z(x_N)\big)$ differ only by a straight line in $x_i$, so that

$$ \mathcal{R}\big(\mathbf{x},\ \mathbf{g}(\mathbf{x})\big) = \mathcal{R}\big(\mathbf{x},\ \mathbf{g}^{0} + z(\mathbf{x})\big) \tag{3.78} $$

The welding deformation therefore acts on block $i$ only through the value of the keel lift $z(x_i)$ at its own position. Moving block $i$ by $\Delta$ changes its gap by approximately $z'(x_i)\,\Delta$, relative to a straight line. This has two consequences.

First, the effect of the deformation on the reactions is not smooth in $\mathbf{x}$: by Equation (3.59), $z'$ jumps by $\varphi_j$ at every erection joint. When a block crosses a joint, the sensitivity of its gap to its position changes abruptly. This is one of the sources of non-smoothness of Problem (3.73) noted in Section 3.12.1.

Second, an erection joint acts as a *hard point* that spacing alone can relieve only partly. For a positive joint angle the kink of Equation (3.58) points towards the blocks, so near the joint the keel is lowest at the joint itself. Relative to the tangent of the smooth part of $z$ at the joint, a block at the distance $d$ from the joint sees the additional gap $\varphi_j\,d/2$. The hull can close this gap only by bending over the short distance $d$, which requires a large local bending moment, so the block nearest to the joint tends to carry a disproportionate share of the load. A scale comparison shows when this governs. The compression of a block that carries exactly its strength is

$$ \delta_{S,i} = \frac{S_i}{k_i} \tag{3.79} $$

For the welded example (blocks of 200 t at 2000 kN/mm) $\delta_S = 0.98$ mm, while the largest joint angle, $\varphi = 0.9$ mrad at $x = 20$ m, gives $\varphi\,d/2 = 1.08$ mm at the given spacing $d = 2.4$ m. A neighbour at the given spacing therefore starts with an additional gap larger than the compression at which the joint block reaches its strength. Reducing the distance to one floor spacing (0.8 m) reduces the additional gap to 0.36 mm, which explains why moderate deformations can be absorbed by the spacing design. The neighbours cannot, however, come closer than $s_{\min}$, and blocks that stand on the lifted parts of the keel near the hull ends remain out of contact whatever their spacing. Beyond a certain deformation, adding blocks or moving them no longer lowers $U$, and the gaps rather than the spacing govern. In the welded example of Chapter 6 the design fits for deformations up to half of the full welding deformation and does not fit for the full deformation, with the governing block at the 20 m joint for every number of blocks tried (Section 6.8). The diagnosis of Section 3.12.6 detects this situation.

### 3.12.3 Phase 1: strength-weighted equidistribution

The first phase computes a good global layout with few analyses. It adapts the equidistribution principle of adaptive mesh generation [@deboor1974; @huang2011] to the block line: the blocks play the role of mesh nodes, the load carried by the blocks plays the role of the monitor function, and the strengths play the role of the target cell weights.

*Start layout.* Unless the end blocks are fixed, the end blocks are placed at the limits of the extent, $x_1 = x_{\text{lo}}$ and $x_N = x_{\text{hi}}$, which shortens the overhangs as far as allowed, and the interior blocks are mapped affinely from the given layout onto $[x_{\text{lo}}, x_{\text{hi}}]$. The layout is then made admissible by the sweep of Equation (3.83) below.

*Load density.* From the reactions $R_i$ of the current layout, the tributary cell of block $i$ is $[e_{i-1}, e_i]$ with

$$ e_0 = x_1, \qquad e_i = \tfrac12\,(x_i + x_{i+1}), \quad i = 1, \dots, N-1, \qquad e_N = x_N \tag{3.80} $$

and the load $\ell_i = \max(R_i, 0)$, raised to at least $5\,\%$ of the mean of the $\ell_i$, is spread uniformly over the cell. The floor keeps the density positive where blocks are lifted off, for example by welding gaps, so that the cumulative load is strictly increasing and can be inverted. The cumulative load $Q(x)$ is the piecewise-linear function through the points

$$ Q(e_0) = 0, \qquad Q(e_i) = \sum_{k=1}^{i} \ell_k, \quad i = 1, \dots, N \tag{3.81} $$

*Equidistribution map.* The new positions divide the total load in proportion to the mean strength of neighbouring blocks:

$$ y_j = Q^{-1}\!\big( F_j\,Q(e_N) \big), \qquad F_j = \frac{\sum_{i < j} w_i}{\sum_{i=1}^{N-1} w_i}, \qquad w_i = \tfrac12\,(S_i + S_{i+1}), \qquad j = 1, \dots, N \tag{3.82} $$

with $y_1 = x_1$ and $y_N = x_N$ (since $F_1 = 0$, $F_N = 1$). Because $Q$ is strictly increasing and the $F_j$ are strictly increasing, the $y_j$ are strictly increasing. The load between two neighbouring new positions is $(F_{j+1} - F_j)\,Q(e_N) \propto w_j$, that is, proportional to the mean strength of the two blocks that share it. Where the hull presses hard or the blocks are weak, the spacing becomes small; where the load is small or the blocks are strong, it becomes large. For equal strengths, Equation (3.82) is the classical equidistribution of the load. At a fixed point with locally uniform spacing, the load between interior neighbours, $(R_j + R_{j+1})/2$, is proportional to $(S_j + S_{j+1})/2$, which is satisfied by $R_j = c\,S_j$, that is, by equal utilisations. The map is a heuristic, not an exact optimality condition: the end blocks also carry the overhang loads, which the tributary view assigns to the end intervals only approximately, and the elastic hull redistributes load when blocks move. Phase 2 corrects both effects.

*Admissible sweep.* The proposed positions are made admissible by one sweep from aft to forward with fixed end blocks. With $x_1$ and $x_N$ fixed, block $j = 2, \dots, N-1$ is kept in the window

$$ \max\!\big(x_{j-1} + s_{\min},\ x_N - (N-j)\,s_{\max}\big) \;\le\; x_j \;\le\; \min\!\big(x_{j-1} + s_{\max},\ x_N - (N-j)\,s_{\min}\big) \tag{3.83} $$

that is, it keeps the spacing limits to its (already placed) predecessor and leaves room for the $N - j$ remaining spacings. Without admissible positions, $x_j$ is the proposed position clipped to the window; with admissible positions (Section 3.12.5), it is the admissible position in the window nearest to the proposed one. Finally the last spacing $x_N - x_{N-1}$ is checked.

*Lemma (the continuous sweep never fails).* If $(N-1)\,s_{\min} \le x_N - x_1 \le (N-1)\,s_{\max}$ and $s_{\min} \le s_{\max}$, the window (3.83) is non-empty for every $j$, and the sweep returns a layout that satisfies all spacing limits.

*Proof.* Consider the invariant

$$ x_N - (N-j+1)\,s_{\max} \;\le\; x_{j-1} \;\le\; x_N - (N-j+1)\,s_{\min} \tag{3.84} $$

For $j = 2$ it is the precondition. If it holds, the four inequalities that make the window non-empty are satisfied: $x_{j-1} + s_{\min} \le x_{j-1} + s_{\max}$ because $s_{\min} \le s_{\max}$; $x_{j-1} + s_{\min} \le x_N - (N-j)\,s_{\min}$ by the right inequality of (3.84); $x_N - (N-j)\,s_{\max} \le x_{j-1} + s_{\max}$ by the left inequality; and $x_N - (N-j)\,s_{\max} \le x_N - (N-j)\,s_{\min}$. Any $x_j$ in the window satisfies (3.84) with $j$ replaced by $j+1$, by the second terms of the max and min in (3.83). By induction, the window is non-empty up to $j = N-1$, and for $j = N-1$ it gives $s_{\min} \le x_N - x_{N-1} \le s_{\max}$. Every spacing $x_j - x_{j-1}$ lies in $[s_{\min}, s_{\max}]$ by the first terms. $\square$

The sweep is also idempotent: a layout that already satisfies the limits lies in every window and is returned unchanged, so it acts as a projection onto the admissible layouts that changes the positions only where necessary. With a discrete set of admissible positions the sweep is greedy and can fail even when an admissible layout exists; this case is treated in Section 3.12.5.

*Fixed-point iteration.* Because the reactions depend on the positions, the map is applied repeatedly: analyse, re-space, analyse again. The full step of Equation (3.82) tends to overshoot, because blocks moved towards a load peak attract more load in the next analysis, which makes the iteration oscillate. The update is therefore under-relaxed,

$$ \mathbf{x} \leftarrow \mathcal{A}\big(\mathbf{x} + \omega\,(\mathbf{y} - \mathbf{x})\big), \qquad \omega = 0.5 \tag{3.85} $$

where $\mathcal{A}$ is the admissible sweep. Convergence of the fixed-point iteration is not guaranteed in general. The phase is therefore protected in three ways: it performs at most 15 iterations, it stops when the largest move is below 1 mm or the sweep fails, and it keeps the best layout found according to the score (3.75), not the last one. For the spacing example of Chapter 6 the phase lowered $U$ from 1.340 (given plan) to 0.787 within its 17 analyses, for the weight-curve example from 1.064 to 0.629 (Section 6.8).

### 3.12.4 Phase 2: sequential linear programming with a trust region

The second phase refines the layout by solving the minimax problem locally. The classical approach for non-smooth minimax problems is to keep the individual functions $u_i$, which are piecewise smooth, rather than the non-smooth maximum, to linearise each of them, and to solve the resulting linear minimax problem as a linear programme in epigraph form [@boyd2004; @fletcher1987]. This sequential linear programming goes back to Griffith and Stewart [@griffith1961] and is widely used in structural optimisation [@haftka1992]. The linear programme chooses the active (governing) blocks itself, so the non-differentiability of the maximum causes no difficulty. A trust region bounds the step, because the linearisation is valid only near the current layout [@nocedal2006; @conn2000].

*Linearisation.* At the current layout $\mathbf{x}$, with utilisations $u_i$, the utilisations of the moved layout $\mathbf{x} + \boldsymbol{\Delta}$ are approximated by

$$ u_i(\mathbf{x} + \boldsymbol{\Delta}) \approx u_i(\mathbf{x}) + \sum_{k} J_{ik}\,\Delta_k, \qquad J_{ik} \approx \frac{\partial u_i}{\partial x_k} \tag{3.86} $$

The sensitivities are computed by one-sided finite differences, one column per movable block,

$$ J_{ik} = \frac{u_i(\mathbf{x} + h\,\mathbf{e}_k) - u_i(\mathbf{x})}{h} \quad \text{or, if } \mathbf{x} + h\,\mathbf{e}_k \text{ is not admissible,} \quad J_{ik} = \frac{u_i(\mathbf{x} - h\,\mathbf{e}_k) - u_i(\mathbf{x})}{-h} \tag{3.87} $$

with $h = 0.05$ m; if neither perturbation is admissible, the column is zero. Each column costs one complete analysis, including the welding gaps at the perturbed position, so the sensitivities contain the effect of the gap change of Section 3.12.2. The step $h$ is large compared with the precision of the positions (1 µm) and of the analysis, and small compared with the floor spacing. Because the reactions are only piecewise smooth, a finite difference over a few centimetres is a secant slope that averages over small kinks, which is preferable here to an exact but local derivative.

*Linear programme.* With the movable blocks $\mathcal{M}$ (all blocks, or the interior blocks when the end blocks are fixed), the step $\boldsymbol{\Delta}$ ($\Delta_k = 0$ for $k \notin \mathcal{M}$), the new positions $x'_i = x_i + \Delta_i$, their spacings $\Delta x'_i = x'_{i+1} - x'_i$ and second differences $\Delta^2 x'_i = x'_{i+2} - 2x'_{i+1} + x'_i$, the linearised subproblem is

$$ \begin{aligned} \min_{\boldsymbol{\Delta},\,t} \quad & t + \mu \sum_{i=1}^{N-2} \left| \Delta^2 x'_i \right| + \varepsilon \sum_{k \in \mathcal{M}} |\Delta_k| \\ \text{s.t.} \quad & u_i + \sum_{k \in \mathcal{M}} J_{ik}\,\Delta_k \le t, && i = 1, \dots, N, \\ & s_{\min} \le \Delta x'_i \le s_{\max}, && i = 1, \dots, N-1, \\ & x'_1 \ge x_{\text{lo}}, \quad x'_N \le x_{\text{hi}}, \\ & |\Delta_k| \le \delta, && k \in \mathcal{M}. \end{aligned} \tag{3.88} $$

The variable $t$ is the predicted largest utilisation (epigraph variable of the maximum). The trust radius $\delta$ bounds the step. The second term, with $\mu = 0.002/\bar s$ per metre and $\bar s = (x_{\text{hi}} - x_{\text{lo}})/(N-1)$ the mean spacing, prefers a spacing that changes gradually: it trades $0.002$ of utilisation for one mean spacing of total variation in the spacing, which is negligible for the objective but selects, among layouts with practically the same largest utilisation, the one that is easiest to set out in the dock. The third term, with $\varepsilon = 10^{-6}$, avoids needless moves.

Problem (3.88) becomes a linear programme by the standard reformulation of absolute values [@boyd2004]. The step is split into its positive and negative parts, $\Delta_k = \Delta_k^{+} - \Delta_k^{-}$ with $\Delta_k^{\pm} \ge 0$, and every absolute second difference is replaced by an auxiliary variable $a_i$ with two inequalities:

$$ \begin{aligned} \min \quad & t + \mu \sum_{i} a_i + \varepsilon \sum_{k} \left(\Delta_k^{+} + \Delta_k^{-}\right) \\ \text{s.t.} \quad & \sum_k J_{ik}\,(\Delta_k^{+} - \Delta_k^{-}) - t \le -u_i, \\ & -(\Delta_{i+1} - \Delta_i) \le (x_{i+1} - x_i) - s_{\min}, \qquad (\Delta_{i+1} - \Delta_i) \le s_{\max} - (x_{i+1} - x_i), \\ & \pm\big(\Delta^2 x_i + \Delta_{i+2} - 2\Delta_{i+1} + \Delta_i\big) - a_i \le 0, \\ & 0 \le \Delta_k^{+} \le \delta_k^{+}, \quad 0 \le \Delta_k^{-} \le \delta_k^{-}, \quad a_i \ge 0, \quad t \ \text{free}, \end{aligned} \tag{3.89} $$

with $\Delta^2 x_i$ the second difference of the current positions, $\delta_k^{\pm} = \delta$ for interior blocks, and $\delta_1^{-} = \min(\delta, x_1 - x_{\text{lo}})$, $\delta_N^{+} = \min(\delta, x_{\text{hi}} - x_N)$ for the end blocks, which are the only blocks that can reach the limits of the extent. At an optimum, $a_i = |\Delta^2 x'_i|$ because $\mu > 0$, and at most one of $\Delta_k^{\pm}$ is non-zero because $\varepsilon > 0$, so Problem (3.89) is equivalent to Problem (3.88). With $m = |\mathcal{M}|$ movable blocks it has $2m + (N-2) + 1$ variables and at most $N + 2(N-1) + 2(N-2) = 5N - 6$ inequality rows. The programme is always feasible, because the current layout ($\boldsymbol{\Delta} = \mathbf{0}$, $t = \max_i u_i$) satisfies all constraints, and bounded, because the step is boxed and $t$ is bounded below by the first rows. It is solved with the HiGHS solver [@huangfu2018].

*Acceptance and trust-region update.* Let $t^{*}$ be the optimal $t$. If the predicted improvement $U(\mathbf{x}) - t^{*}$ is below $10^{-5}$, the linearised optimum is reached and the phase ends. Otherwise the trial layout $\mathbf{x} + \boldsymbol{\Delta}$ is formed, the end blocks are clipped to the extent, and, if round-off makes a spacing smaller than $s_{\min}$, the layout is repaired by the admissible sweep (3.83); if the repair fails, the radius is halved without an analysis. The trial is then analysed with the full model and accepted only if it improves the score (3.75):

$$ \text{accept if } \operatorname{score}(\mathbf{x} + \boldsymbol{\Delta}) < \operatorname{score}(\mathbf{x}): \quad \delta \leftarrow \min(2\delta,\ \bar s), \qquad \text{otherwise}: \quad \delta \leftarrow \delta / 2 \tag{3.90} $$

The initial radius is $\delta = \bar s/4$. After an accepted step the sensitivities are recomputed at the new layout; after a rejected step they are reused, and only the linear programme is solved again with the smaller radius, which costs no analyses. The phase stops when $\delta < 5$ mm, when the predicted improvement is below $10^{-5}$, when the linear programme fails, or after $I_S$ iterations (20 by default). Equation (3.90) is the simple-decrease variant of the trust-region method: instead of comparing the actual with the predicted reduction, every actual decrease is accepted. Since the score decreases strictly with every accepted step and the radius halves with every rejected one, the phase cannot cycle, and the number of analyses is bounded by $I_S\,(m + 1)$.

### 3.12.5 Phase 3: admissible positions at the floors and local search

Blocks must stand under transverse structure. The admissible positions $\mathcal{G}$ are therefore the structure (floor) positions of the input, sorted, when snapping to structure is enabled (the default); otherwise, if a grid step $\Delta_g$ is given, they are the multiples of $\Delta_g$ within the extent,

$$ \mathcal{G} = \left\{ m\,\Delta_g \ :\ m \in \mathbb{Z},\ x_{\text{lo}} \le m\,\Delta_g \le x_{\text{hi}} \right\} \tag{3.91} $$

and without either the continuous layout of phase 2 is the result.

*Rounding.* If the end blocks are free, each is moved to the nearest admissible position within the extent. The interior blocks are then rounded by the admissible sweep (3.83), each to the admissible position in its window nearest to its continuous position. A rounded layout satisfies all constraints of Problem (3.73) by construction. If the sweep finds no admissible position in some window, phase 3 is skipped and the continuous layout is kept; the position check of Section 3.8.6 then reports the blocks that are off the floors.

*Local search.* Rounding to floors spaced 0.8 m apart can raise $U$ noticeably, because a single block moved by up to 0.4 m changes its share of the load. The rounded layout is therefore improved by a best-improvement local search on $\mathcal{G}$. Its neighbourhood $\mathcal{N}(\mathbf{x})$ is built around the governing blocks:

$$ \mathcal{H} = \big\{ i \in \text{three largest } u_i \ :\ u_i \ge 0.98\,U(\mathbf{x}) \big\}, \qquad \mathcal{K} = \big\{ i-1,\ i,\ i+1 \ :\ i \in \mathcal{H} \big\} \cap \mathcal{M} \tag{3.92} $$

For every block $k \in \mathcal{K}$, every direction (aft or forward), every step size (one or two admissible positions) and every group, either block $k$ alone or blocks $k$ and $k+1$ together (if both are movable), a trial layout is formed by moving each block of the group by the given number of admissible positions from its own current position. A trial is discarded if a block runs off $\mathcal{G}$, if it leaves the extent or violates the spacing limits, or if it repeats an earlier trial. Hence $|\mathcal{N}(\mathbf{x})| \le 8\,|\mathcal{K}| \le 72$. The design of this neighbourhood follows from the minimax structure: only moves near the governing blocks can lower the maximum, so moving other blocks would spend analyses without effect on the objective. Moves of two positions allow a block to pass a position that is blocked by the spacing limits, and moving a block together with its forward neighbour shifts load along the line without changing the spacing between the two.

All trials of the neighbourhood are analysed, and the best is accepted if it improves the score (3.75); otherwise the layout is a local optimum with respect to $\mathcal{N}$ and the search ends. The search also ends when its budget of $25N$ analyses is used up. Because the score decreases strictly with every accepted move and $\mathcal{G}$ is finite, the search terminates even without the budget, which bounds the cost in the worst case. Section 6.8 reports the effect of this phase for the examples.

### 3.12.6 Infeasibility diagnosis and suggested number of blocks

When the final layout has $U > u_{\text{target}}$, the plan does not fit, and the result must say why, because the remedy depends on the cause. The diagnosis rests on a lower bound that holds for every support model and every layout.

*Lemma (strength bound).* For any layout in equilibrium with $\sum_i R_i = W$ and $R_i \ge 0$,

$$ U(\mathbf{x}) \;\ge\; \frac{\sum_i R_i}{\sum_i S_i} \;=\; \frac{W}{\sum_{i=1}^{N} S_i} \tag{3.93} $$

with equality if and only if all utilisations are equal.

*Proof.* By the definition of $U$, $R_i = u_i\,S_i \le U\,S_i$ for every $i$. Summing gives $W = \sum_i R_i \le U \sum_i S_i$. Equality requires $u_i = U$ for every block with $S_i > 0$, that is, for all blocks. $\square$

Equation (3.93) is the discrete form of the obvious requirement that the blocks together must be able to carry the ship. It also shows that equal utilisations are the ideal that the design aims at. For the spacing example of Chapter 6 the bound is $4500/8090 = 0.556$ and the design reached 0.761; for the weight-curve example the bound is $4643/8800 = 0.528$ and the design reached 0.568. The difference is caused by the end blocks, which also carry the overhangs, by the minimum spacing, which limits how densely blocks can be placed under heavy regions, and by the elastic hull.

The program distinguishes three causes and formulates the message accordingly:

1. **Total strength too small.** If $u_{\text{target}} \sum_i S_i < \sum_i \max(R_i, 0)$, Equation (3.93) proves that no spacing can satisfy the target. The message states the total strength and the load.
2. **End block at the limit of the extent.** If the governing block is the first or the last block and the end blocks are not fixed, the end block stands at the limit $x_{\text{lo}}$ or $x_{\text{hi}}$ (phase 1 placed it there), and its load is dominated by the overhang. The message recommends a longer block line, that is, a smaller $x_{\text{lo}}$ or a larger $x_{\text{hi}}$, to shorten the overhang.
3. **Blocks out of contact.** If blocks of the final layout carry no load ($R_i \le 0$), they are held off the keel by gaps, typically by the welding deformation. The gaps rather than the spacing then govern (Section 3.12.2), and more blocks may not help. The message recommends adjusting the block heights to the keel line, which is outside the scope of the present program, or using stronger blocks where the keel bears.

Causes 1 and 2 are combined in one reason text; cause 3 selects the wording of the recommendation. In every case of a failed design the number of blocks probably needed is suggested,

$$ N_{\text{sug}} = \max\!\left( N + 1,\ \left\lceil \frac{N\,U}{u_{\text{target}}} \right\rceil \right) \tag{3.94} $$

The estimate assumes that the load shared by the blocks scales inversely with their number, $U' \approx N U / N_{\text{sug}}$, so that $N_{\text{sug}}$ blocks of the same kind would bring $U'$ to the target. It is an estimate, not a guarantee: overhangs and gaps do not scale with $N$, which is why the message for cause 3 qualifies it ("if the load spread"). If the plan fits, the message states that all blocks are within strength, the value of $U$ and the number of analyses.

### 3.12.7 Closed-form optimum for a uniform load

A case with a known optimum shows what the method should achieve and serves as a verification benchmark (Section 5.6).

*Proposition.* Let a hull of length $L$ carry the uniform load $q$ on $[0, L]$, and let $N$ blocks of equal strength $S$ support it. Then for every layout and every support model in equilibrium

$$ U(\mathbf{x}) \;\ge\; U^{*} = \frac{q\,L}{N\,S} \tag{3.95} $$

The bound is attained by the layout with equal spacing and half-spacing overhangs,

$$ x_i^{*} = \left(i - \tfrac12\right) \frac{L}{N}, \qquad x_{i+1}^{*} - x_i^{*} = \frac{L}{N}, \qquad i = 1, \dots, N \tag{3.96} $$

in which every block carries exactly $qL/N$, both for the rigid hull and for the tributary-area model in which every block carries the load of its cell between the midpoints to its neighbours, the end blocks including the overhangs. All utilisations are then equal to $U^{*}$.

*Proof.* The lower bound is Equation (3.93) with $W = qL$ and $\sum_i S_i = N S$. For the tributary model, the cell of block $i$ in the layout (3.96) is $[(i-1)L/N,\ iL/N]$, with the overhangs $[0, L/(2N)]$ and $[L - L/(2N), L]$ included in the first and last cells, so every cell has the length $L/N$ and carries $qL/N$. For the rigid hull, the block centroid of the layout (3.96) is $L/2$, which coincides with the LCG of the uniform load, so $e = 0$ in Equation (3.30) and every block carries $W/N = qL/N$ (equal stiffness, no gaps). In both cases $U = qL/(NS) = U^{*}$. $\square$

The layout (3.96) is feasible for Problem (3.73) if $s_{\min} \le L/N \le s_{\max}$ and the extent contains $[L/(2N),\ L - L/(2N)]$. For the elastic hull, the bound (3.95) still holds, but equal reactions are attained only approximately, because the hull redistributes load near its free ends. The proposition shows in particular that the optimum is not the layout with the end blocks at the hull ends: each end block must stand half a spacing inboard so that its cell, including the overhang, has the same length as the others. The automated test of Section 5.6 starts from a random layout of 11 blocks and confirms that the method reaches $U^{*}$ within 1 % with all utilisations within 2 % of each other.

### 3.12.8 Algorithm

The complete method is summarised in Algorithm 3.3. The function `analyse` is the injected analysis: mesh and load model for the candidate positions, welding gaps by Equation (3.76), and the contact solution of the selected support model (Algorithms 3.1 and 3.2). Every call is counted. The default parameters are listed in Table 3.3.

**Algorithm 3.3 – Spacing design from block strengths.**

```text
input : blocks 1..N in order with strengths S_i and plan gaps g0_i, given positions x0,
        extent [x_lo, x_hi], s_min, s_max (default 2 × mean given spacing), u_target,
        fix_end_blocks, admissible positions G (floors or grid, optional), I_S
output: positions x, spacing x_{i+1} − x_i, status, message, history
 1  reject if (N−1) s_min > x_hi − x_lo, or fixed ends and (N−1) s_max < x_hi − x_lo
 2  if (N−1) s_max < x_hi − x_lo: x_hi ← x_lo + (N−1) s_max            (free ends)
 3  U(x0) ← analyse(x0)                                     (given plan, for reference)
 4  x ← x0; unless fix_end_blocks: x_1 ← x_lo, x_N ← x_hi, interior mapped affinely
 5  x ← A(x) (Eq. 3.83); analyse(x); best ← x
    — phase 1: strength-weighted equidistribution —
 6  repeat at most 15 times
 7      y ← equidistribution map of the reactions of x (Eqs. 3.80–3.82)
 8      y ← A(x + 0.5 (y − x)); stop if A fails or max|y − x| < 1 mm
 9      x ← y; analyse(x); if score(x) < score(best): best ← x        (Eq. 3.75)
10  x ← best
    — phase 2: sequential linear programming with trust region —
11  s̄ ← (x_hi − x_lo)/(N−1); δ ← s̄/4; μ ← 0.002/s̄; J ← none
12  for it = 1 … I_S
13      if J is none: J ← finite differences (Eq. 3.87)       (≤ m analyses)
14      (Δ, t*) ← solve LP (Eq. 3.89) with radius δ; stop if it fails
15      stop if U(x) − t* < 1e-5
16      y ← x + Δ; clip end blocks to the extent; if a spacing < s_min: y ← A(y)
17      if A failed: δ ← δ/2; continue
18      analyse(y)
19      if score(y) < score(x): x ← y; δ ← min(2δ, s̄); J ← none
20      else: δ ← δ/2; stop if δ < 5 mm
    — phase 3: admissible positions and local search —
21  if G is given:
22      y ← end blocks to nearest position of G in the extent (unless fixed); y ← A_G(y)
23      if y exists:
24          x ← y; analyse(x); budget ← 25 N
25          repeat
26              T ← neighbourhood of x (Eq. 3.92: governing blocks ± 1, moves of 1 or 2
                    positions, alone or with the next block, admissible only)
27              analyse every z in T while the budget lasts; z* ← best score
28              if z* is none or score(z*) ≥ score(x): stop
29              x ← z*
30  U ← max_i R_i(x)/S_i
31  if U ≤ u_target: status ← fits
32  else: status ← does not fit; N_sug ← max(N+1, ⌈N U/u_target⌉);
          reason ← total strength (Eq. 3.93) / end block at the extent limit / open blocks
33  return x, spacing, status, message, history of all phases
```

| Parameter | Symbol | Default | Purpose |
|---|---|---|---|
| Target utilisation | $u_{\text{target}}$ | 1.0 | acceptance of the layout |
| Minimum spacing | $s_{\min}$ | 1.0 m | practicable block distance |
| Maximum spacing | $s_{\max}$ | $2 \times$ mean given spacing | no long unsupported spans |
| Extent | $[x_{\text{lo}}, x_{\text{hi}}]$ | given block line | limits of the block line |
| Equidistribution iterations, relaxation | –, $\omega$ | 15, 0.5 | phase 1 |
| Finite-difference step | $h$ | 0.05 m | sensitivities, phase 2 |
| Trust radius: initial, maximum, minimum | $\delta$ | $\bar s/4$, $\bar s$, 5 mm | phase 2 |
| Smoothing weight, move penalty | $\mu$, $\varepsilon$ | $0.002/\bar s$, $10^{-6}$ | phase 2 |
| SLP iterations | $I_S$ | 20 | phase 2 |
| Governing tolerance, number of governing blocks | – | 2 %, 3 | phase 3 neighbourhood |
| Local-search budget | – | $25N$ analyses | phase 3 |

Table: Table 3.3 – Parameters of the spacing design and their default values.

## 3.13 Computational complexity

This section analyses the time and memory requirements of the algorithms of Sections 3.3–3.12 as they are implemented in the program (Chapter 4). The analysis uses the following symbols:

- $n_e$ number of beam elements and $n = 2(n_e + 1)$ number of degrees of freedom;
- $n_w$ number of weight items (including the items obtained from a weight curve) and $n_s$ number of hull sections;
- $N$ number of blocks and $m \le N$ number of movable blocks in the spacing design;
- $I$ number of iterations of the block-state iteration;
- $p$ number of survey stations and $q$ number of identified parameters;
- $T_A$ the cost of one complete support analysis of a layout.

For the example vessel, $n_e = 255$, $n = 512$, $N = m = 44$, $n_w = 9$ and $n_s = 3$. The measured run times that confirm the predictions of this section are reported in Section 5.9.

### 3.13.1 Load model, assembly and factorisation

**Mesh and weight distribution.** A weight curve of $n_c$ points is converted into at most $n_c - 1$ items by Equation (3.12) in $O(n_c)$. The mesh generator sorts the key points and subdivides the intervals between them, which costs $O(n_e + n_p \log n_p)$ for $n_p$ key points. Each weight item is then spread over the elements it covers with vectorised array operations over all elements, so the load model costs $O(n_w\,n_e)$ time and $O(n_e)$ memory. The exact integrals of Equations (3.16)–(3.18) used in the weight check are single passes over the elements, $O(n_e)$.

**Assembly.** The element matrices of Equation (3.38) are computed for all elements at once as an array of shape $n_e \times 4 \times 4$. The $16\,n_e$ entries are written in coordinate (COO) format and converted to compressed sparse row (CSR) format, which sums the overlapping entries of neighbouring elements. Both steps are linear, $O(n_e)$ in time and memory. The consistent load vector is accumulated in the same way.

**Factorisation.** The global matrix of the free–free beam couples only the two nodes of each element. With the degrees of freedom numbered node by node, every row has at most six non-zero entries, and all of them lie within a distance of three from the diagonal. The matrix is therefore banded with half-bandwidth $b = 3$, and the block springs only add to diagonal entries. LU factorisation of a banded matrix creates fill-in only inside the band [@davis2006]. The sparse direct solver SuperLU used by SciPy [@li2005superlu] does not exploit the band explicitly: it applies a fill-reducing column ordering and threshold partial pivoting. For a matrix of this one-dimensional structure, the factors remain of order $O(n)$ in size, so that both factorisation and solution cost

$$ T_{\text{LU}} = O(n\,b^2) = O(n), \qquad M_{\text{LU}} = O(n\,b) = O(n) \tag{3.97} $$

operations and memory. Section 5.9 confirms the near-linear growth by measurement.

### 3.13.2 Contact iteration and gap continuation

Every iteration of Algorithm 3.1 linearises the $N$ block laws ($O(N)$), adds the tangent stiffnesses to the diagonal and solves one sparse system. In the implementation, the matrix is factorised anew in each iteration, because the tangent stiffnesses change. One iteration of the elastic model therefore costs $O(n + N) = O(n)$, and $I$ iterations cost $O(I\,n)$. The rigid-hull model needs only the sums of Equation (3.29) and a $2 \times 2$ solve, $O(N)$ per iteration.

The number of iterations is bounded, but not by a small polynomial in general. Each block can be open, elastic or crushed, so there are at most $3^N$ state sets. Cycle detection prevents the simultaneous update from revisiting a state set, and the one-at-a-time fallback changes one block per iteration. In the implementation, a hard limit of $I_{\max} = 200$ iterations per solve is applied; if it is reached, or if the support becomes unstable, the gap continuation of Algorithm 3.2 is used. The continuation repeats the state iteration for a sequence of gap vectors; its cost is the sum of the iterations over all continuation steps. In practice the counts are small: the example vessel converged in one iteration, the welded hull in three (elastic) and five (rigid) iterations, and the hardest case found, a rigid hull with uniform 0.2 mrad joint kinks, in 11 iterations over one plus three continuation steps (Section 3.11).

**One analysis.** A complete support analysis of a layout, as performed for every candidate of the spacing design, rebuilds the mesh and the load model (the block positions are mesh nodes), evaluates the keel lift at the $N$ block positions with $J$ joints and $n_z$ shrinkage zones and fits the line of Equation (3.63), and solves the contact problem:

$$ T_A = O\big(n_w\,n_e + N\,(J + n_z) + I\,n\big) \tag{3.98} $$

which is linear in the mesh size for a bounded number of contact iterations.

### 3.13.3 Spacing design

The spacing design is dominated by the number of analyses. Counting the calls of `analyse` in Algorithm 3.3:

- the given plan and the start layout: 2 analyses;
- phase 1: at most 15 analyses;
- phase 2: per iteration at most $m$ analyses for the sensitivities (only after an accepted step) and one for the trial, so at most $I_S\,(m + 1)$ analyses;
- phase 3: one analysis of the rounded layout and at most $25N$ analyses in the local search, each step evaluating at most 72 trials.

The total number of analyses is therefore bounded by

$$ n_A \;\le\; 17 + I_S\,(m + 1) + 1 + 25\,N \tag{3.99} $$

which is linear in $N$ for a fixed iteration limit; for $N = m = 44$ and $I_S = 20$ the bound is 2018 analyses. In practice far fewer are needed, because the trust region converges before the iteration limit, columns of $\mathbf{J}$ are skipped where a block cannot move, and the local search stops at a local optimum: the spacing example of Chapter 6 used 193 analyses and the weight-curve example 501 (Section 6.8).

Every SLP iteration also solves one linear programme (3.89) with

$$ n_{\text{var}} = 2m + (N - 2) + 1 \le 3N - 1, \qquad n_{\text{row}} \le 5N - 6 \tag{3.100} $$

which is 131 variables and 214 rows for the example. The constraint matrix is assembled densely, $O(N^2)$ entries. The programme is solved by HiGHS, by default with the dual revised simplex method [@huangfu2018]. The simplex method has exponential worst-case complexity, but on practical problems the number of iterations typically grows roughly linearly with the number of rows and columns; at this size a solve takes milliseconds and is negligible compared with the analyses for the sensitivities. The finite-difference sensitivities cost $O(m\,T_A)$ per Jacobian, which is the dominant term of phase 2. Together, the spacing design costs

$$ T_{\text{spacing}} = O\big( n_A\,T_A + I_S\,T_{\text{LP}}(3N, 5N) \big), \qquad n_A = O(I_S\,N) \tag{3.101} $$

with memory $O(N^2)$ for the sensitivities and the linear programme, in addition to the $O(n)$ of one analysis. Because the analyses are independent, the finite-difference columns and the trials of the local search could be evaluated in parallel; the present implementation evaluates them sequentially.

### 3.13.4 Identification and post-processing

**Identification.** The least-squares problem of Equation (3.68) has $p$ survey stations and $q = J + 2$ unknowns ($J$ joint angles and the two straight-line coefficients), plus $J$ rows when regularisation is used. It is solved by an orthogonal factorisation in $O(p\,q^2)$ time and $O(p\,q)$ memory. For the sample survey, this is $25 \times 7$, which is negligible. The identification is performed once per case, not per candidate layout.

**Strength post-processing.** Shear force and bending moment are integrated element by element with a fixed number of sample points per element (eight in the implementation, with both one-sided values at every node), $O(n_e)$. The section properties and envelope limits are evaluated at all sample points, $O(n_s\,n_e)$. The distance of every block to the nearest transverse structure is computed with a dense $N \times n_{\text{st}}$ comparison for $n_{\text{st}}$ structure positions (151 floors in the example). The post-processing is performed for the given and for the designed layout, not for the candidates of the spacing design.

### 3.13.5 Summary of the complexity

Table 3.4 summarises the analysis. Every step of a single analysis is linear in the size of the mesh, apart from the number of contact iterations, which is small in practice. The spacing design multiplies this cost by the number of analyses, which grows linearly with the number of blocks.

| Algorithm | Time | Memory | Dominant term for the examples |
|---|---|---|---|
| Weight curve to items | $O(n_c)$ | $O(n_c)$ | 18 points, 12 items |
| Mesh and load model | $O(n_w\,n_e)$ | $O(n_e)$ | $9 \times 255$ element updates |
| Assembly (COO → CSR) | $O(n_e)$ | $O(n_e)$ | $16 \times 255$ matrix entries |
| Sparse LU and solve | $O(n)$ | $O(n)$ | $n = 512$, half-bandwidth 3 |
| Contact iteration | $O(I\,n)$ | $O(n)$ | $I = 1$ (example), $I = 3$ (welded hull) |
| Rigid-hull iteration | $O(I\,N)$ | $O(N)$ | $N = 44$ |
| One complete analysis | $T_A$, Eq. (3.98) | $O(n)$ | mesh, welding gaps, contact |
| Spacing design, phase 1 | $\le 17\,T_A$ | $O(N)$ | 17 analyses |
| Spacing design, phase 2 | $\le I_S\,(m+1)\,T_A + I_S\,T_{\text{LP}}$ | $O(N^2)$ | LP 131 variables × 214 rows |
| Spacing design, phase 3 | $\le (1 + 25N)\,T_A$ | $O(N)$ | ≤ 72 trials per step |
| Spacing design, total | $n_A \le 17 + I_S(m+1) + 1 + 25N$ analyses | $O(N^2 + n)$ | 193 and 501 analyses |
| Survey identification | $O(p\,q^2)$ | $O(p\,q)$ | $25 \times 7$ |
| Strength post-processing | $O(n_s n_e + N n_{\text{st}})$ | $O(n_e)$ | 2297 sample points |

Table: Table 3.4 – Time and memory complexity of the algorithms and their dominant terms for the examples.

## 3.14 Summary

This chapter has established the methodology used in this thesis:

- The ship's weight is represented by trapezoidally spread items, which preserve each item's weight and centre of gravity, subject to the middle-third condition. An item can be positioned by its LCG or by its longitudinal moment, with a consistency check when both are given. A weight-distribution curve with jumps is converted segment by segment into trapezoidal items; the conversion is exact and satisfies the middle-third condition automatically. On a mesh containing every item boundary, all load integrals are exact (Equations (3.4)–(3.18)).
- Each block is a compression-only spring whose stiffness follows from its layers in series, with an optional bilinear crushing branch; for a fixed contact state the block law is affine. The blocks are defined by a block plan that gives the number of blocks and the spacing, strength and gap of every block (Equations (3.19)–(3.26)).
- Method A treats the hull as rigid and reduces the problem to a $2 \times 2$ system; for equal stiffness and continuous block lines it reproduces the classical trapezoidal and triangular distributions (Equations (3.27)–(3.32)).
- Method B models the hull as a free–free Timoshenko or Euler–Bernoulli beam on block springs, discretised with Hermite elements and consistent loads and solved by a sparse direct method. Its nodal results are independent of mesh density (Equations (3.33)–(3.44)).
- The nonlinear support conditions are resolved by an active-set state iteration with cycle detection, a one-at-a-time fallback and a support-stability check (Algorithm 3.1, Equations (3.45)–(3.47)), and large, uneven gaps by a continuation in the gaps (Algorithm 3.2, Equations (3.70)–(3.71)).
- Shear force and bending moment are integrated exactly element by element, with a closure check. Stresses, envelope utilisations, the utilisation of every block with respect to its own strength, bearing pressures and positions relative to the floors are then evaluated (Equations (3.48)–(3.57)).
- The welding deformation of the hull is described by joint kinks, block shrinkage curvatures and an interpolated survey, and converted into block gaps; a straight-line part has no influence on the reactions. Joint angles can be identified from a keel survey by regularised least squares (Equations (3.58)–(3.69)).
- The spacing design computes the positions of the given blocks that minimise the largest ratio of reaction to strength, subject to spacing limits, the extent of the block line and the floor positions (Problem (3.73)). Every candidate is analysed with the full support model and the welding gaps at its own positions, so the welding deformation is a constraint of the design; erection joints act as hard points that limit what spacing can achieve (Equations (3.76)–(3.79)). The problem is solved by a strength-weighted equidistribution with a provably admissible sweep, a sequential linear programme with a trust region, and rounding to the floors with a best-improvement local search (Algorithm 3.3, Equations (3.80)–(3.92)). If the plan does not fit, a lower bound and the state of the final layout identify the cause, and a block count is suggested (Equations (3.93)–(3.94)). For a uniform load the optimum is known in closed form (Equations (3.95)–(3.96)).
- The time and memory cost of every algorithm was analysed; a single analysis is linear in the mesh size, and the number of analyses of the spacing design is linear in the number of blocks (Section 3.13, Table 3.4).

Chapter 4 describes how these models were implemented in software, and Chapter 5 how they were verified.
