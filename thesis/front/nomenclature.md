# Nomenclature {.unnumbered}

## Latin symbols {.unnumbered}

| Symbol | Description | Unit |
|---|---|---|
| $A$ | contact area of a block, $A = b\,\ell$ | m² |
| $A_s$ | effective shear area of the hull girder | m² |
| $a,\ b$ | rigid-body translation and rotation parameters of Method A | m, rad |
| $b,\ \ell$ | contact width and contact length of a block | m |
| $c$ | constant term of the linearised block law $R = k_t\delta + c$ | N |
| $E,\ G$ | Young's modulus, shear modulus | Pa |
| $e$ | eccentricity of the LCG from the centroid of the block line | m |
| $EI$ | bending stiffness of the hull girder | N·m² |
| $F_y$ | crush force of a block, $F_y = \sigma_c A$ | N |
| $\mathbf{F}$ | global load vector | N, N·m |
| $g_i$ | initial gap between keel line and top of block $i$ (plan gap plus welding gap at $x_i$) | m |
| $h_{\max}$ | maximum element length of the mesh | m |
| $I$ | second moment of area of the hull girder section | m⁴ |
| $I_S$ | maximum number of sequential-LP iterations of the spacing design | – |
| $J_{ik}$ | sensitivity of the utilisation of block $i$ to the position of block $k$, $\partial u_i / \partial x_k$ | 1/m |
| $\mathbf{K}$, $\mathbf{K}_e$ | global and element stiffness matrices | N/m, … |
| $k_i$ | initial stiffness of block $i$ | N/m |
| $k_2$ | post-crush stiffness of a block, $k_2 = \rho\,k$ | N/m |
| $k_f$ | equivalent foundation modulus of a block line, $k/s$ | N/m² |
| $k_t$ | tangent stiffness of a block in its current state | N/m |
| $L$ | element length; length of a weight item | m |
| $L_{pp}$ | length between perpendiculars | m |
| $LCG$, $\mathrm{LCG}_j$ | longitudinal centre of gravity of the ship, of weight item $j$ (from AP) | m |
| $M(x)$ | hull girder bending moment, sagging positive | N·m |
| $M_j$ | longitudinal moment of weight item $j$ about the AP, $M_j = W_j\,\mathrm{LCG}_j$ (input: t·m) | N·m |
| $N$ | number of blocks of the block plan | – |
| $N_{\text{sug}}$ | suggested number of blocks when a plan does not fit, $\lceil N\,U/u_{\text{target}} \rceil$ | – |
| $N_1 \dots N_4$ | Hermite shape functions | – |
| $n,\ n_e$ | number of degrees of freedom, number of beam elements | – |
| $P_j$ | point weight | N |
| $p$ | block bearing pressure, $R/A$; (Section 3.13) number of survey stations | Pa; – |
| $Q(x)$ | cumulative load along the block line used by the equidistribution step | N |
| $q(x),\ w(x)$ | load intensity (weight per unit length), downward positive; $q(x)$ also the input weight-distribution curve | N/m |
| $(x_k,\ q_k)$ | points of the piecewise-linear weight-distribution curve | m, N/m |
| $R_i$ | reaction of block $i$, upward on the hull | N |
| $S_i$ | strength of block $i$: the force it may carry (block-plan value, else capacity of its type) | N |
| $s$ | block spacing (generic); local coordinate in an element | m |
| $s_i$ | spacing between block $i$ and block $i+1$, $s_i = x_{i+1} - x_i$ (result of the spacing design) | m |
| $\bar s$ | mean spacing of the block line | m |
| $s_{\min},\ s_{\max}$ | minimum and maximum spacing between neighbouring blocks (default $s_{\max}$: twice the given mean spacing) | m |
| $t$ | epigraph variable of the minimax problem (bound on all utilisations) | – |
| $U$ | largest utilisation of the block line, $U = \max_i u_i$ (objective of the spacing design) | – |
| $U^\ast$ | optimal largest utilisation; closed form $qL/(NS)$ for a uniform load on $N$ equal blocks | – |
| $u_i$ | utilisation of block $i$, $u_i = R_i / S_i$ | – |
| $u_{\text{target}}$ | target utilisation; a plan fits when $U \le u_{\text{target}}$ (default 1) | – |
| $\mathbf{u}$ | global displacement vector | m, rad |
| $V(x)$ | hull girder shear force | N |
| $v(x)$ | vertical hull displacement, downward positive | m |
| $W$, $W_j$ | ship weight; weight of item $j$ | N |
| $w_a,\ w_b$ | end intensities of a trapezoidal weight item | N/m |
| $w_i$ | strength weight of the interval between blocks $i$ and $i+1$, $\tfrac12(S_i + S_{i+1})$ | N |
| $x$ | longitudinal coordinate from the aft perpendicular, forward positive | m |
| $x_i$ | position of block $i$ | m |
| $x_{a,j},\ x_{b,j}$ | aft and forward end of the extent of weight item $j$ | m |
| $x_{\text{lo}},\ x_{\text{hi}}$ | aft and forward limit of the allowed extent of the block line in the spacing design | m |
| $x_j$ | (Sections 3.9–3.10) position of erection joint $j$ | m |
| $z(x)$ | upward keel lift of the welded hull | m |
| $Z_{\text{deck}},\ Z_{\text{keel}}$ | section moduli at deck and keel | m³ |

## Greek symbols {.unnumbered}

| Symbol | Description | Unit |
|---|---|---|
| $\alpha$ | Tikhonov regularisation parameter (Section 2.11) | – |
| $\Delta_k$ | step of the position of block $k$ in a sequential-LP iteration | m |
| $\Delta^2 x_i$ | second difference of the block positions, $x_{i+2} - 2x_{i+1} + x_i$ | m |
| $\delta_i$ | compression of block $i$, $\delta_i = v_i - g_i$ | m |
| $\delta$ | trust radius of the sequential linear programming (spacing design) | m |
| $\delta_y$ | compression at which crushing starts, $F_y/k$ | m |
| $\varepsilon_{\text{deck}},\ \varepsilon_{\text{bottom}}$ | welding shrinkage strain of deck and bottom | – |
| $\eta_i$ | violation of the assumed state of block $i$ (contact iteration) | m |
| $\theta$ | rotation of the hull girder axis, $\mathrm{d}v/\mathrm{d}x$ | rad |
| $\kappa$ | welding-shrinkage curvature of a construction block; shear correction factor in $\kappa G A$ | 1/m; – |
| $\kappa GA$ | shear stiffness of the hull girder | N |
| $\lambda$ | characteristic parameter of a beam on elastic foundation, $(k_f/4EI)^{1/4}$; continuation parameter of the gap homotopy | 1/m; – |
| $\lambda_R$ | Tikhonov weight of the joint-angle identification from a keel survey (Section 3.10) | – |
| $\mu$ | smoothing weight of the spacing variation in the sequential LP, $0.002/\bar s$ (Section 3.12.4) | 1/m |
| $\xi$ | normalised element coordinate, $s/L$; integration variable | –; m |
| $\Pi$ | potential energy of the hull–block system | J |
| $\rho$ | post-crush stiffness ratio $k_2/k$; (Section 2.9.4) monitor function | –; – |
| $\sigma$ | normal stress, tension positive | Pa |
| $\sigma_c$ | crush stress of a block | Pa |
| $\tau$ | shear stress | Pa |
| $\Phi$ | shear-deformation parameter, $12EI/(\kappa G A L^2)$ | – |
| $\varphi_j$ | angular distortion (kink) at erection joint $j$ | rad |
| $\psi$ | section rotation of the Timoshenko beam | rad |
| $\psi(x)$ | unit-curvature shape function of a construction block | m² |
| $\omega$ | under-relaxation factor of the equidistribution fixed-point iteration (0.5) | – |

## Abbreviations {.unnumbered}

| Abbreviation | Meaning |
|---|---|
| AP / FP | aft / forward perpendicular |
| BM / SF | bending moment / shear force |
| CLI | command-line interface |
| COO / CSR | coordinate / compressed sparse row matrix formats |
| CSV | comma-separated values |
| EB | Euler–Bernoulli (beam theory) |
| FD | finite differences |
| FEM | finite element method |
| GUI | graphical user interface |
| HiGHS | high-performance linear optimisation solver |
| LCP | linear complementarity problem |
| LP | linear programme |
| PCHIP | piecewise cubic Hermite interpolating polynomial |
| SLP | sequential linear programming |
| SWBM / SWSF | still-water bending moment / shear force |
| YAML | YAML Ain't Markup Language (input file format) |
