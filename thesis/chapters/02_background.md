# 2 Theoretical Background and Literature Review

This chapter reviews the engineering practice and the theory on which the present work is built. It draws on two disciplines. The ship engineering background supplies the application domain: dry-docking practice (Section 2.1), block construction and welding distortion in shipbuilding (Section 2.2), hull girder longitudinal strength (Section 2.3), the classical docking calculations (Section 2.4) and the structural mechanics used to go beyond them: beam theory, beams on elastic foundations, the finite element method and unilateral contact (Sections 2.5 to 2.8). The computer science background — optimisation for design, including minimax problems, sequential linear programming and equidistribution (Section 2.9), algorithms and complexity in computational mechanics (Section 2.10), interpolation and inverse problems (Section 2.11), software engineering for engineering applications (Section 2.12) and verification and validation (Section 2.13) — supplies the algorithmic and software methods that form the core of the work. Section 2.14 identifies the gap that this thesis addresses.

## 2.1 Dry docking practice

### 2.1.1 Types of docking facility

Three principal types of facility are used to take ships out of the water: graving docks, floating docks and ship lifts or slipways [@mazurkiewicz1980]. A **graving dock** is a basin built into the shore and closed by a gate or caisson. The ship is floated in, the gate is closed and the water is pumped out, so that the ship settles onto supports laid on the dock floor. The floor of a graving dock is a massive structure and is normally treated as rigid in block-load calculations. Design criteria for such docks, including the floor and the loads it must sustain, are the subject of facility design standards such as @ufc421310. Graving docks are also widely used as **building docks**, in which a new ship is assembled on keel blocks and floated out on completion [@storch1995].

A **floating dock** is a pontoon structure with side walls that is ballasted down to receive the ship and then deballasted to lift it. Its pontoon deck is itself a flexible structure. The distribution of block loads therefore depends not only on the ship and the blocks but also on the stiffness of the dock and on how its ballast is distributed. **Ship lifts and slipways** transfer the ship onto cradles or platforms on land and raise similar questions of support distribution.

The present work is concerned with the support of the ship on its blocks. It treats the dock floor as rigid, which is appropriate for graving docks and building docks. For floating docks this is an idealisation, and the limitation was stated in Section 1.4.

### 2.1.2 The docking plan and the docking sequence

Each ship has a **docking plan**: a drawing that gives the ship's profile and bottom contours, the positions where blocks may and may not be placed, the locations of hull openings, appendages, sonar domes and docking keels, and the docking weight and centre of gravity. Docking practice uses this plan to lay out the blocks in advance of the docking [@nstm997]. Naval docking practice is codified in technical manuals such as the U.S. Navy's *Naval Ships' Technical Manual*, Chapter 997 [@nstm997]. The certification of the facilities and procedures used to dock naval ships is governed by @milstd1625.

A docking proceeds in stages. The ship is brought over the block line and aligned with it, the water level is lowered, and the ship settles onto the blocks, first at one end and then along its full length. Once the water has been removed, the full weight is carried by the blocks. This fully docked condition, which also describes a hull standing on its building blocks, is the condition considered in this thesis (Section 1.4).

### 2.1.3 Keel blocks, side blocks and timber caps

Keel blocks are placed along the centreline under the keel or flat keel plate. They typically consist of a base of concrete or steel, intermediate layers of hardwood, and a softwood **cap** that bears on the hull [@mazurkiewicz1980; @nstm997]. The cap has two functions. It conforms to small irregularities of the hull surface and so spreads the contact pressure. It is also the most flexible layer of the stack and therefore governs the vertical stiffness of the block. Timber loaded perpendicular to its grain is much less stiff and much weaker than timber loaded along the grain. Beyond a certain stress it deforms progressively and permanently, which is referred to as **crushing**.

Side blocks, or bilge blocks, are placed outboard of the centreline under the bottom or the turn of the bilge. They increase the support area, reduce the load on the keel line and provide transverse support. Their share of the load depends on the transverse stiffness of the bottom structure and on how closely they fit the hull. That share is a three-dimensional effect, which lies outside the longitudinal model of the present work.

Blocks are preferably positioned under transverse structure such as floors, transverse bulkheads and web frames. The concentrated reaction is then introduced into structure able to carry it, rather than into unsupported plating [@nstm997]. The admissible block positions therefore form a discrete set — in a typical hull, the floors at a regular frame spacing — and any design of the block spacing must end at positions of that set.

### 2.1.4 Why block loads must be calculated

Block loads must be known for three reasons. First, each block has a limited **capacity**, both as a whole and as an allowable **bearing pressure** on its contact area. Overloading a block crushes the cap and can transfer load progressively to neighbouring blocks. Second, the **hull girder** must carry the concentrated reactions together with its own weight distribution, and the resulting shear forces and bending moments must not exceed permissible values. Third, the **local structure** at each block, meaning the keel, floors and bottom plating, must be able to carry the reaction. Docking standards therefore require block-load calculations as part of the preparation of a docking [@milstd1625; @nstm997].

The blocks of a dock or yard are not all alike. Their capacity depends on their construction and their condition, so that the force a block may carry — its **strength** in the sense of this thesis — can differ from block to block along the line. Because the block loads depend on the block positions, the calculation is not only a check of a given arrangement: it is the basis for *choosing* the arrangement. The traditional practice of laying out the blocks by rules of thumb — a uniform spacing, closer spacing under heavy machinery and at the ends of the block line — and checking the result afterwards [@nstm997; @mazurkiewicz1980] leaves this choice to trial and error.

## 2.2 Block construction and welding distortion in shipbuilding

### 2.2.1 Block construction

Modern merchant and naval ships are built by **block construction**. The hull is divided into construction blocks (also called units, sections or grand blocks) that are fabricated, and to a large extent outfitted, in workshops and on assembly areas before being transported to the building dock or berth [@storch1995; @okumoto2009]. There the blocks are positioned on keel blocks, aligned and joined by welding at the **erection joints**, usually transverse butts through the bottom, side shell, decks and longitudinal members. Block construction shortens the building time and moves much of the welding to controlled workshop conditions, but it places high demands on dimensional accuracy, because each block must fit its neighbours at the erection joints [@storch1995].

### 2.2.2 Welding distortion

Fusion welding heats the material locally to above its melting point. As the weld and the surrounding metal cool, they contract, but the contraction is restrained by the colder material around them. The result is a combination of residual stresses and permanent deformation, known as **welding distortion** [@ueda2012]. The principal forms relevant to a ship hull are:

- **longitudinal shrinkage** along a weld line, which shortens the welded members;
- **transverse shrinkage** across a weld line, which shortens the structure in the direction normal to the weld;
- **angular distortion**, a rotation of the plates on either side of a weld caused by a variation of transverse shrinkage through the plate thickness; and
- **bending distortion** of members whose welds lie eccentrically to their neutral axis.

A widely used way of predicting such distortion in large structures is the concept of **inherent strain**: the permanent strains left in and near the weld are treated as an initial strain field, which is applied to an elastic model of the structure to compute the resulting deformation [@ueda2012]. The present work does not compute the inherent strains themselves. It adopts the idea that, once the welding is complete, its effect on the global shape of the hull can be represented by equivalent kinematic quantities.

### 2.2.3 Consequences for the keel line

At the scale of the hull girder, two effects of welding are important for the support on keel blocks. First, the welding of an erection joint produces a small **angular distortion** of the hull girder at that joint: the two adjoining blocks are no longer exactly in line, and the keel line acquires a kink. Second, if the upper part of the hull (decks and upper side shell) shrinks more than the bottom, for example because it contains more weld metal per unit length, each construction block **bends** so that its keel becomes concave upward. With a girder depth $D$ and longitudinal shrinkage strains $\varepsilon_\text{deck}$ and $\varepsilon_\text{bottom}$, simple beam kinematics gives the corresponding curvature

$$ \kappa = \frac{\varepsilon_\text{deck} - \varepsilon_\text{bottom}}{D} \tag{2.1} $$

which is positive (concave upward) when the deck shrinks more than the bottom.

Both effects accumulate along the hull. Even small kink angles of the order of a fraction of a milliradian at each of several joints, and small curvatures in each block, add up to deviations of the keel line of several centimetres at the ends of a long hull. If the hull is then supported on keel blocks whose tops lie on a straight line, the hull ends lift off the end blocks. Shipbuilding practice therefore pays considerable attention to accuracy control and to the measurement of the keel line [@storch1995; @okumoto2009]. For the calculation of block loads, the welding distortion appears as a set of **initial gaps** between the keel and the block tops, which can be large compared with the elastic deflections of the hull. Because the gap at a block is the keel lift *at the position of that block*, the gaps change whenever the blocks are moved. For a design of the block positions the welding deformation is therefore not a fixed input of the analysis but a function of the design variables.

## 2.3 Hull girder longitudinal strength

### 2.3.1 The ship as a beam

For the analysis of overall (primary) strength, a ship's hull is idealised as a thin-walled box girder, the **hull girder**, that bends in the vertical plane under the difference between the distributed weight and the distributed support [@hughes2010; @lewis1988]. In still water the support is buoyancy, and the net load per unit length is $q(x) = w(x) - b(x)$, where $w(x)$ is the weight distribution and $b(x)$ the buoyancy distribution. The shear force $V(x)$ and bending moment $M(x)$ follow by integration along the length. They are the **still-water shear force** and **still-water bending moment**. In waves, the wave-induced components are added to them. A bending moment that makes the hull concave upward is called **sagging**, and one that makes it convex upward **hogging**.

The weight distribution $w(x)$ is traditionally built from grouped weights — hull structure, machinery, outfit, tank contents — each spread over its extent so that its weight and its centre of gravity are preserved, or it is taken directly from the weight curve of the loading manual [@hughes2010; @lewis1988]. Weight lists commonly state the longitudinal moment of each item about a reference point instead of its centre of gravity, from which the centre follows by division by the weight.

On blocks the buoyancy is replaced by the block reactions, so the net load becomes the weight distribution minus a series of concentrated upward forces. The same integration then gives the hull girder loads on the blocks. The shear-force diagram is therefore saw-toothed, with a jump at every block, and the bending moment is governed by how closely the reaction distribution follows the weight distribution.

### 2.3.2 Bending stress and section modulus

Under the assumptions of simple beam theory (Section 2.5), the longitudinal bending stress at a distance $z$ from the neutral axis of the hull girder cross-section is

$$ \sigma = \frac{M\,z}{I} \tag{2.2} $$

where $I$ is the second moment of area of the effective longitudinal material about the neutral axis [@hughes2010]. The extreme stresses occur at the deck and at the keel. They are conveniently expressed through the **section moduli** $Z_\text{deck} = I/z_\text{deck}$ and $Z_\text{keel} = I/z_\text{keel}$. In sagging the deck is in compression and the keel in tension, and in hogging the reverse holds. The shear stress in the side shell and longitudinal bulkheads is governed by the shear force. To a first approximation it can be estimated by dividing the shear force by an effective shear area.

### 2.3.3 Permissible values

Classification societies set minimum requirements for hull girder strength. The common requirements for longitudinal strength of seagoing ships are laid down in the IACS Unified Requirement S11 [@iacss11]. For each ship, the approved **loading manual** gives permissible still-water bending moments and shear forces, usually as envelopes along the length with separate limits for hogging and sagging and, for some ships, separate harbour and seagoing values. In a docking calculation, the hull girder loads caused by the block reactions are compared with these permissible values. The comparison is meaningful only if the reactions themselves are realistic. A method that misrepresents how the reactions are distributed misrepresents the bending moment as well, as the results of this thesis show (Chapter 6). The same envelopes apply to any redesigned block arrangement, because moving the blocks redistributes the reactions and therefore changes the bending moment.

## 2.4 Classical docking calculations

### 2.4.1 The rigid-hull (trapezoidal) method

The traditional method for estimating keel-block loads treats the ship as a **rigid body** resting on a continuous, uniformly elastic block line [@nstm997]. Because the hull is rigid, the compression of the block line, and therefore the load per unit length on it, varies linearly along its length. The two unknowns of the linear distribution are fixed by the two equilibrium equations, Equations (1.1) and (1.2). For a continuous block line of length $L_b$ whose centroid lies a distance $e$ aft or forward of the ship's LCG, the load per unit length at the two ends of the block line is

$$ q_\text{end} = \frac{W}{L_b}\left(1 \pm \frac{6e}{L_b}\right) \tag{2.3} $$

provided $|e| \le L_b/6$. If the eccentricity exceeds one sixth of the block-line length, the formula would give a negative load at one end, which blocks cannot provide. The distribution then becomes triangular over only part of the block line, and the remainder of the blocks lifts off. This is the same "middle-third" result that governs the bearing pressure under a rigid footing.

The method is attractive because it requires only the weight, the LCG and the block-line geometry, and it can be evaluated by hand. Its essential limitation follows from its rigid-hull assumption: it cannot represent the redistribution of load that occurs when the hull bends. A real hull deflects more under heavy concentrated weights such as machinery and at the ends of the block line, where the overhanging bow and stern act as cantilevers. In practice such effects have to be covered by additional allowances or by engineering judgement. The method also has no means of representing a keel line that is not straight: a rigid hull with a curved keel touches only a few blocks, which is a poor description of a real, flexible hull. For the design of the block spacing the rigid method is of limited use for the same reason: the reactions it predicts depend on the block positions only through two parameters, so it cannot show where along the line a closer or wider spacing is needed. The present work addresses these effects directly by modelling the hull as an elastic beam with initial gaps.

### 2.4.2 Discrete blocks of unequal stiffness

The rigid-hull idea generalises readily to discrete blocks of unequal stiffness. If the hull displacement is $v(x) = a + bx$ and block $i$ has stiffness $k_i$, the two equilibrium equations become two linear equations in $a$ and $b$, and each reaction follows as $R_i = k_i(a + b x_i)$. With equal stiffnesses this reduces to the familiar result that the reaction varies linearly with the distance of the block from the centroid of the block line. The method is still a two-parameter model of the support distribution, however detailed the block data. Chapter 3 develops it fully, including gaps and crushing (Section 3.5).

## 2.5 Beam theory

### 2.5.1 Euler–Bernoulli beam

The classical theory of bending, attributed to Euler and Bernoulli, assumes that plane cross-sections remain plane and perpendicular to the deformed axis, and that deflections and rotations are small [@gere2013]. For a beam of bending stiffness $EI(x)$ under a distributed transverse load $q(x)$, the deflection $v(x)$ satisfies the fourth-order equation

$$ \frac{d^2}{dx^2}\!\left(EI\,\frac{d^2 v}{dx^2}\right) = q(x) \tag{2.4} $$

and the bending moment is related to the curvature by $M = -EI\,v''$ for the conventions of this thesis (downward displacement and sagging moment positive, Section 3.2). The shear force is the derivative of the bending moment. These relations underlie the hull girder idealisation of Section 2.3.

### 2.5.2 Timoshenko beam

For beams that are short or deep relative to their span, or in which the shear stiffness is low compared with the bending stiffness, the deformation caused by shear is no longer negligible. Timoshenko's theory [@timoshenko1921] relaxes the assumption that cross-sections remain perpendicular to the deformed axis. It introduces the cross-section rotation $\psi$ as an independent field, with the shear strain given by the difference between the slope of the axis and the rotation. The shear force is then proportional to this strain through the shear stiffness $\kappa G A$, where $G$ is the shear modulus, $A$ the cross-sectional area and $\kappa$ a shear correction factor accounting for the non-uniform distribution of shear stress. For a ship hull girder the effective shear area is essentially that of the side shell and longitudinal bulkheads [@hughes2010]. The relative importance of shear deformation over an element of length $L$ is measured by the dimensionless parameter

$$ \Phi = \frac{12\,EI}{\kappa G A\,L^2} \tag{2.5} $$

which vanishes in the Euler–Bernoulli limit. For a whole ship hull the parameter is small, but over the short distance between neighbouring blocks it is not. Shear deformation therefore influences how load is shared between adjacent blocks.

## 2.6 Beams on elastic foundations

A beam supported continuously by a medium whose reaction is proportional to the local deflection is called a **beam on an elastic (Winkler) foundation**. The classical treatment is that of @hetenyi1946. With a foundation modulus $k_f$ (force per unit length per unit deflection) and constant $EI$, the governing equation becomes

$$ EI\,\frac{d^4 v}{dx^4} + k_f\,v = q(x) \tag{2.6} $$

Its solutions are characterised by the parameter

$$ \lambda = \left(\frac{k_f}{4\,EI}\right)^{1/4} \tag{2.7} $$

which has the dimension of an inverse length. A concentrated load on such a beam produces a deflection that decays and oscillates with distance, and the load is spread over a length of the order of $\pi/\lambda$ on each side of its point of application [@hetenyi1946]. Two limits follow. If $1/\lambda$ is large compared with the length of the beam, the beam behaves as rigid on its foundation, which is the rigid-hull method of Section 2.4.1. If $1/\lambda$ is small, the foundation reaction follows the applied load closely and each load is carried locally.

A line of closely spaced keel blocks of stiffness $k$ at spacing $s$ approximates a Winkler foundation with $k_f \approx k/s$. Equation (2.7) therefore gives a direct physical measure of whether the hull of a docked ship behaves rigidly or flexibly relative to its blocks. For the illustrative example vessel of this thesis, the load-spreading length $\pi/\lambda$ is about 37 m, roughly one third of the block line. This explains why the rigid-hull method and the elastic model give markedly different results (Chapter 6). The relation $k_f \approx k/s$ also shows why the spacing is an effective design variable: a locally closer spacing makes the foundation locally stiffer and draws load to that region, but it spreads the load over more blocks, so that each block carries less. It also means that the hull can bend far enough under its own weight to close part of an initial gap at its ends, so that the effect of a welding-induced end lift cannot be judged from the gap alone. The continuous-foundation idealisation does not, however, represent individual blocks, gaps, block crushing or lift-off. These are represented exactly in a discrete model with block springs, which this thesis adopts.

## 2.7 Matrix structural analysis and the finite element method

The **displacement (stiffness) method** of matrix structural analysis expresses the equilibrium of a structure discretised into elements in the form

$$ \mathbf{K}\,\mathbf{u} = \mathbf{F} \tag{2.8} $$

where $\mathbf{u}$ is the vector of nodal displacements, $\mathbf{F}$ the vector of nodal loads and $\mathbf{K}$ the global stiffness matrix assembled from element stiffness matrices [@przemieniecki1968]. For a plane beam element with a transverse displacement and a rotation at each of its two nodes, the element stiffness matrix can be derived exactly from beam theory. This holds for both Euler–Bernoulli and Timoshenko assumptions [@przemieniecki1968; @cook2002].

The **finite element method** generalises this approach by approximating the displacement field within each element with shape functions and deriving the element matrices from a variational principle [@cook2002; @bathe2014]. For the Euler–Bernoulli beam, the cubic Hermite shape functions reproduce the exact homogeneous solution of the beam equation. When distributed loads are converted into **consistent nodal loads** using the same shape functions, the nodal displacements of a beam with prismatic elements are then exact, independent of the number of elements [@cook2002]. This property is relevant to the present work. It means that, for piecewise-linear loads on piecewise-prismatic hulls, mesh refinement serves only to resolve the deflected shape between nodes and does not change the nodal results.

Supports modelled as springs add their stiffness to the diagonal of $\mathbf{K}$. A beam that is held only by springs, as a docked hull is, has two rigid-body modes (vertical translation and rotation). It is therefore only stably supported if at least two springs at different positions are active. Because the stiffness matrix of a beam is banded and sparse, Equation (2.8) is solved efficiently by sparse direct methods [@davis2006]. The SciPy library [@virtanen2020scipy] provides such solvers, based on the SuperLU package [@li2005superlu]. The low cost of one solution is essential for the present work, because the design of the block spacing requires hundreds of complete analyses, each with a new mesh and new block positions.

## 2.8 Unilateral contact and complementarity

A block can push on the hull but cannot pull it. Its reaction is therefore subject to a **unilateral constraint**: either the block is in contact and its reaction is compressive, or it is out of contact and its reaction is zero. For rigid supports with an initial gap $g_i$, the conditions take the form of a **complementarity problem**:

$$ R_i \ge 0, \qquad g_i - v_i \ge 0, \qquad R_i\,(g_i - v_i) = 0 \tag{2.9} $$

that is, at each support either the reaction or the remaining gap is zero. When the structure is linear elastic, the problem of finding displacements and reactions that satisfy equilibrium together with Equation (2.9) is a **linear complementarity problem** (LCP). The mathematical theory of the LCP, including existence and uniqueness of solutions and pivoting and iterative solution algorithms, is treated comprehensively by @cottle1992. For elastic supports, the reaction is a piecewise-linear function of the compression $\delta_i = v_i - g_i$, as in Equation (1.3). The same structure appears, with the springs regularising the contact.

In computational contact mechanics more generally, unilateral constraints are enforced by **penalty**, **Lagrange multiplier** or **augmented Lagrangian** methods, often combined with an **active-set strategy** that iterates on which constraints are in contact [@wriggers2006]. For the docking problem the block springs are physical, not numerical penalties. Their stiffness is a property of the block, so a penalty-type formulation is in this case exact rather than approximate. Because the block law is piecewise linear, an active-set iteration on the state of every block (open, elastic or crushed) solves a linear problem for each trial set of states. It terminates with a solution that satisfies the block law exactly once the assumed states are consistent with the computed compressions. Simple active-set iterations are, however, known to be sensitive to the starting point and can cycle when many constraints change state at once [@wriggers2006]. This is precisely the situation created by large initial gaps, and Chapter 3 therefore combines the iteration with a continuation in the gaps (Section 3.11).

A consequence for design is that the block reactions are continuous but only piecewise smooth functions of the block positions: their derivatives jump whenever a block changes its state. Any optimisation of the positions must cope with this non-smoothness.

## 2.9 Optimisation for design

Structural optimisation seeks the values of design variables — member sizes, shapes, or the positions of supports — that minimise an objective subject to constraints on stresses, displacements and geometry, with every evaluation requiring a structural analysis [@haftka1992]. The optimal positioning of supports is a classical branch of this field: the support positions of beams and frames have been optimised for minimum bending moment, minimum compliance or maximum fundamental frequency [@mroz1975; @wang2004]. Two features of such problems recur. The objective is often the largest of many local responses, which leads to minimax problems (Section 2.9.2), and the responses are available only through an analysis, which leads to methods that solve a sequence of approximate subproblems (Section 2.9.3). A third, independent line of work — the equidistribution principle of adaptive mesh methods — addresses the closely related problem of placing a given number of points along a line so that each carries an equal share of a density (Section 2.9.4). All three are building blocks of the spacing design of Section 3.12; their common computational tool is linear programming.

### 2.9.1 Linear programming

A **linear programme** (LP) minimises a linear objective subject to linear equality and inequality constraints and bounds on the variables [@nocedal2006; @boyd2004]:

$$ \min_{\mathbf{x}}\ \mathbf{c}^\mathsf{T}\mathbf{x} \quad \text{subject to} \quad \mathbf{A}\,\mathbf{x} \le \mathbf{b}, \qquad \mathbf{l} \le \mathbf{x} \le \mathbf{u} \tag{2.10} $$

The feasible set of an LP is a convex polyhedron, every local minimum is a global minimum, and if an optimum exists it is attained at a vertex of the polyhedron. LPs with thousands of variables and constraints are solved reliably by **simplex** and **interior-point** methods [@nocedal2006]. LP problems are convex optimisation problems, and their optimality conditions and duality theory are well established [@boyd2004]. **HiGHS** is an open-source solver for large-scale LP problems; its dual revised simplex implementation is described by @huangfu2018. HiGHS is accessible from Python through the optimisation module of SciPy [@virtanen2020scipy], which is how it was used in the present work.

### 2.9.2 Minimax problems

Many engineering design criteria are of the **minimax** type: the largest of several quantities, for example the largest utilisation of a set of components, is to be made as small as possible,

$$ \min_{\mathbf{x}}\ F(\mathbf{x}), \qquad F(\mathbf{x}) = \max_{i = 1, \dots, N} f_i(\mathbf{x}) \tag{2.11} $$

Even when every $f_i$ is smooth, $F$ is not differentiable wherever two or more of the $f_i$ attain the maximum, and at a minimax optimum this is typically the case: several functions are *active* and equal, and no direction lowers all of them at once [@fletcher1987]. Gradient methods for smooth problems applied directly to $F$ therefore stall or zig-zag. The standard remedy is the **epigraph formulation**, which introduces an auxiliary variable $t$ as an upper bound of all functions [@boyd2004]:

$$ \min_{\mathbf{x},\,t}\ t \quad \text{subject to} \quad f_i(\mathbf{x}) \le t, \quad i = 1, \dots, N \tag{2.12} $$

This problem has a smooth objective and smooth constraints; when every $f_i$ is affine it is an LP of the form (2.10). An objective based on absolute values, $\sum_i |f_i(\mathbf{x})|$, is linearised in a similar way by bounding each absolute value by an auxiliary variable from above and below.

The characterisation of the optimum by several equal active functions has a familiar structural counterpart. A **fully stressed design**, in which every member of a structure is loaded to its allowable stress, is a classical design criterion and the basis of simple resizing algorithms [@haftka1992]. For the docking problem the analogue is a block arrangement in which every block uses the same share of its strength. Where such a layout is attainable, it minimises the largest utilisation; where constraints such as a minimum or maximum spacing prevent it, the optimum equalises the utilisation of the blocks that govern. Both views are used in Section 3.12.

### 2.9.3 Sequential linear programming and trust regions

When the functions $f_i$ are nonlinear, the epigraph problem (2.12) can be solved by **sequential linear programming** (SLP), introduced by @griffith1961 for process design and widely used in structural optimisation [@haftka1992]. At the current point $\mathbf{x}_k$ each function is replaced by its first-order Taylor expansion, and the step $\boldsymbol{\Delta}$ is obtained from the linear programme

$$ \min_{\boldsymbol{\Delta},\,t}\ t \quad \text{subject to} \quad f_i(\mathbf{x}_k) + \nabla f_i(\mathbf{x}_k)^\mathsf{T} \boldsymbol{\Delta} \le t, \qquad \mathbf{A}(\mathbf{x}_k + \boldsymbol{\Delta}) \le \mathbf{b}, \qquad \|\boldsymbol{\Delta}\|_\infty \le \delta_k \tag{2.13} $$

in which linear constraints of the original problem are kept exactly. The last condition is essential. Without it the linearised problem is often unbounded or produces steps far outside the region in which the linearisation is valid, and the iteration oscillates. In structural optimisation such bounds are known as **move limits** [@haftka1992]; in the theory of nonlinear optimisation they define a **trust region** of radius $\delta_k$ [@conn2000; @nocedal2006]. The trust radius is adapted from the agreement between the decrease predicted by the linear model, $F(\mathbf{x}_k) - t_k$, and the actual decrease $F(\mathbf{x}_k) - F(\mathbf{x}_k + \boldsymbol{\Delta})$: the step is accepted and the radius enlarged if the actual decrease is satisfactory, and the step is rejected and the radius reduced otherwise. For minimax problems this combination of an LP subproblem in the epigraph form with an $\ell_\infty$ trust region is a well-established method with global convergence properties [@fletcher1987; @conn2000].

The gradients $\nabla f_i$ are rarely available in closed form when each $f_i$ is the output of a nonlinear analysis. They can be approximated by **finite differences**, at the cost of one additional analysis per design variable [@nocedal2006; @haftka1992]. The step size is a compromise between truncation error and round-off error and, for piecewise-smooth responses such as contact reactions, between capturing the local slope and stepping across a change of contact state. The cost of an SLP iteration is thus dominated by the analyses for the sensitivities, which makes the number of analyses the natural measure of the efficiency of the method.

SLP converges to a local optimum. For non-convex problems a good starting point is therefore important, both for the quality of the result and for the number of iterations, and SLP is commonly preceded by a fast heuristic that brings the design into the neighbourhood of the optimum.

### 2.9.4 Equidistribution in adaptive mesh methods

The **equidistribution principle** was introduced by @deboor1974 for the choice of the knots of spline approximations and has become the basis of adaptive moving-mesh methods [@huang2011]. Given a positive **monitor function** $\rho(x)$ that measures where resolution is needed, the nodes $x_1 < x_2 < \dots < x_N$ of a one-dimensional mesh on $[x_1, x_N]$ are placed so that every cell contains the same share of its integral:

$$ \int_{x_i}^{x_{i+1}} \rho(\xi)\,\mathrm{d}\xi = \frac{1}{N-1} \int_{x_1}^{x_N} \rho(\xi)\,\mathrm{d}\xi, \qquad i = 1, \dots, N-1 \tag{2.14} $$

With the cumulative function $P(x) = \int_{x_1}^{x} \rho\,\mathrm{d}\xi$, which is increasing, the nodes follow by inverting $P$ at equally spaced values. The mesh is dense where $\rho$ is large and coarse where it is small. When the monitor function depends on the solution, and the solution on the mesh, equidistribution becomes a **fixed-point iteration**: compute the solution on the current mesh, evaluate the monitor function, redistribute the nodes, and repeat. Practical moving-mesh methods stabilise this iteration by smoothing the monitor function over neighbouring cells and by under-relaxing the node movement, and they impose bounds on the ratio of neighbouring cell sizes [@huang2011]. The analogy with the support problem is direct: the blocks play the role of the mesh nodes, the load to be carried plays the role of the monitor function, and equal shares of the load per cell correspond to equal utilisations of equal blocks. For blocks of different strength the shares are weighted by the strengths (Section 3.12.3). The fixed-point iteration is not an optimisation method and offers no guarantee of optimality, but it is global, cheap and robust, which makes it a natural starting point for the SLP refinement of Section 2.9.3.

## 2.10 Algorithms and complexity in computational mechanics

The models reviewed above become useful only through algorithms whose cost and reliability are understood. Computer science judges an algorithm by its correctness, its time and memory complexity as a function of the problem size, and its behaviour on the inputs that occur in practice. This section summarises the results that are relevant to the algorithms of Chapter 3; their complexity in the present application is analysed in Section 3.13.

### 2.10.1 Sparse direct solution of banded systems

The stiffness matrix of a beam discretised into elements couples each node only with its neighbours. It is therefore **banded**: all non-zero entries lie within a narrow band around the diagonal, whose width does not grow with the number of elements. Gaussian elimination (LU or Cholesky factorisation) of a banded matrix produces no fill-in outside the band. For a matrix of dimension $n$ and half-bandwidth $b$ the factorisation costs of the order of $n\,b^2$ operations and $n\,b$ words of storage, and every subsequent solution with a new right-hand side costs of the order of $n\,b$ operations [@davis2006]. For a one-dimensional beam model, where $b$ is a small constant, both factorisation and solution are thus linear in the number of nodes. General sparse direct solvers such as SuperLU [@li2005superlu] do not assume a band structure; they reduce fill-in by reordering the unknowns and gain efficiency by operating on dense blocks of columns (supernodes). For matrices that are already banded they retain the linear behaviour. Because every change of the block states changes the diagonal of the system, each step of the contact iteration requires a new factorisation; its linear cost is what makes the complete analysis, and hence the hundreds of analyses of a spacing design, affordable.

### 2.10.2 Complexity of linear programming

The **simplex method** moves from vertex to vertex of the feasible polyhedron. Its worst-case number of iterations grows exponentially with the problem size, but on practical problems it typically needs a number of iterations of the order of the number of constraints, and it is the method of choice for many medium-sized problems [@nocedal2006]. **Interior-point methods** have polynomial worst-case complexity and are competitive for large problems [@nocedal2006; @boyd2004]. Modern solvers such as HiGHS combine carefully engineered implementations of the dual revised simplex method with presolve techniques that reduce the problem before solving it [@huangfu2018]. The linear programmes of the spacing design have of the order of $3N$ variables and $5N$ constraints for $N$ blocks, that is, some hundreds for a typical block line, and such solvers return their solution in a small fraction of the time needed for the analyses that produce their data.

### 2.10.3 Active-set methods and continuation

Active-set methods for complementarity problems guess which constraints are active, solve the resulting linear system, and update the guess [@cottle1992; @wriggers2006]. Each iteration is cheap, but the number of possible active sets grows exponentially with the number of constraints, and a naive update rule can cycle between sets. Termination therefore has to be ensured by the update rule itself, for example by changing one constraint at a time once cycling is detected. **Continuation** (or homotopy) methods address the dependence on the starting point in a different way: they embed the difficult problem in a family of problems controlled by a parameter, start from a member whose solution is easy, and follow the solution as the parameter is changed in small steps, so that each step needs only a few iterations from a good starting point. Chapter 3 combines both techniques for the contact problem with large welding gaps (Sections 3.7 and 3.11). For the spacing design their robustness is a precondition: an analysis that fails for some candidate layout would make the objective undefined there and break the optimisation.

## 2.11 Interpolation and inverse problems

### 2.11.1 Interpolation of measured profiles

A surveyed keel line consists of heights measured at discrete stations. To evaluate it at the block positions it must be interpolated. Piecewise-linear interpolation is robust but has kinks at every station. Cubic splines are smooth but can overshoot between data points, producing artificial waves that the measured shape does not contain. **Monotone piecewise cubic Hermite interpolation** (PCHIP) chooses the derivatives at the data points so that the interpolant is monotone wherever the data are monotone [@fritsch1980]. It is therefore well suited to measured shapes in which spurious oscillations would be misinterpreted as physical deviations. For a design that moves the blocks, the interpolant must be evaluated at arbitrary positions, not only at the surveyed stations or at the blocks of the given plan.

### 2.11.2 Least squares and regularisation

Instead of interpolating a survey, one can fit the parameters of a physical model to it. If the model is linear in its parameters $\mathbf{p}$ — as a sum of kinks at known joint positions is — the fit is a **linear least-squares problem**. When the data are noisy and the columns of the model matrix are nearly dependent, the problem is **ill-conditioned**: small measurement errors produce large errors in the estimated parameters. **Tikhonov regularisation** stabilises the estimate by adding a penalty on the size of the parameters [@hansen1998]:

$$ \min_{\mathbf{p}}\ \left\| \mathbf{A}\,\mathbf{p} - \mathbf{d} \right\|_2^2 + \alpha\,\left\| \mathbf{p} \right\|_2^2 \tag{2.15} $$

where $\mathbf{d}$ is the vector of measured values and $\alpha \ge 0$ the regularisation parameter, which trades the fit to the data against the stability of the solution [@hansen1998]. The identification of welding parameters from a keel survey in Chapter 3 (Section 3.10) is such an inverse problem. Its reliability depends on the number and spacing of survey stations relative to the joint positions and on the measurement noise. An identified model has the further advantage over a direct interpolation that it describes the keel line by a few physical parameters, which can be evaluated at any candidate block position.

## 2.12 Software engineering for engineering applications

### 2.12.1 Requirements and architecture

Software engineering is concerned with the systematic specification, design, implementation, verification and evolution of software [@sommerville2016]. A first step is the formulation of **requirements**, usually divided into functional requirements (what the system must compute or do) and non-functional requirements such as accuracy, performance, usability, maintainability and portability [@sommerville2016]. For engineering calculation software, traceability of every result to its inputs and to the documented method is a further essential requirement.

The **architecture** of a system describes its decomposition into components and their dependencies. A **layered architecture** organises the components so that each layer uses only the services of the layers below it [@sommerville2016]. @martin2017 formulates the related **dependency rule**: source-code dependencies should point from the details — user interfaces, file formats, frameworks — towards the policies that embody the essential logic of the application, and never the reverse. Applied to an engineering application, the calculation core is then independent of the graphical interface and of the input format. It can be tested in isolation, reused from scripts, and given new interfaces without being changed. The same principle applies inside the core: an optimisation algorithm that receives the analysis as a function argument (**dependency injection**) does not depend on the structural model it optimises, and it can be tested with a simple substitute analysis whose optimum is known [@martin2017].

### 2.12.2 Testing and maintainability

Testing is the main technique for demonstrating that a program does what it is intended to do and for discovering defects [@sommerville2016]. **Unit tests** check individual components, **integration tests** check their interaction, and **regression tests** — usually the whole automated suite — are re-run after every change to ensure that it has not broken existing behaviour. For numerical software, meaningful tests compare computed results with analytical solutions, limiting cases or independent computations, with tolerances chosen to reflect round-off rather than to hide errors. For optimisation software, tests additionally check that the result satisfies all constraints and that it reaches a known optimum on problems constructed to have one. Frameworks such as pytest [@pytest] make such suites easy to write and run. Readable code, explicit units and conventions, input validation at the boundary and a documented data model contribute to maintainability, which matters particularly for software that will be extended by others, as is the case for academic engineering tools.

## 2.13 Verification and validation

The credibility of results from a computational model rests on two distinct activities [@oberkampf2010]. **Verification** establishes that the computer program solves the equations of the mathematical model correctly. It covers code verification, by comparison with analytical or manufactured solutions, and solution verification, by estimation of numerical errors such as discretisation and round-off. **Validation** establishes how accurately the mathematical model represents the physical reality it is meant to describe, by comparison with experimental or full-scale measurements.

The two activities answer different questions, and neither substitutes for the other. A program can be perfectly verified and yet represent reality poorly. A model can capture reality well and yet be implemented incorrectly. Systematic verification is a prerequisite for meaningful validation, because without it a disagreement with measurements cannot be attributed to the model rather than to the code [@oberkampf2010]. For a design method a third question arises: whether the design algorithm actually finds a good solution of the design problem. It is answered, as far as possible, by benchmark problems with a known optimum.

In software practice, verification is supported by **automated tests** that check computed results against known solutions every time the code changes (Section 2.12.2). In the present work such tests were written with the pytest framework [@pytest]. As stated in Section 1.4, this thesis addresses verification. Validation against measured block loads and keel shapes remains a subject for future work.

## 2.14 Research gap and positioning of this work

The literature and practice reviewed in this chapter provide all the theoretical ingredients for the analysis and the design of the support of a ship on keel blocks: the hull girder idealisation and its strength criteria [@hughes2010; @lewis1988; @iacss11], the classical docking calculations [@nstm997; @mazurkiewicz1980], beam theory including shear deformation [@gere2013; @timoshenko1921], beams on elastic foundations [@hetenyi1946], the finite element method [@przemieniecki1968; @cook2002; @bathe2014], the treatment of unilateral contact [@cottle1992; @wriggers2006], the understanding of welding distortion in block-built hulls [@ueda2012; @storch1995; @okumoto2009], optimisation for design, including minimax problems, sequential linear programming with trust regions, support-position optimisation and equidistribution [@haftka1992; @fletcher1987; @conn2000; @boyd2004; @nocedal2006; @huangfu2018; @deboor1974; @huang2011], the treatment of inverse problems [@fritsch1980; @hansen1998] and the principles of software engineering and verification [@sommerville2016; @martin2017; @oberkampf2010]. These ingredients are, however, usually encountered separately. Table 2.1 compares the typical approaches to determining block loads and block arrangements.

| Approach | Hull stiffness | Discrete blocks, gaps, lift-off | Block crushing | Welded (curved) keel line | Block spacing | Effort and transparency |
|---|---|---|---|---|---|---|
| Rigid-hull (trapezoidal) hand calculation | rigid | continuous line; lift-off by triangle rule | no | no | rules of thumb, checked afterwards | low effort; fully transparent |
| Beam on continuous elastic foundation | elastic | smeared; no individual blocks | no | no | no | moderate; analytical insight |
| Elastic beam on discrete contact springs | elastic, incl. shear | yes | yes (piecewise-linear) | if gaps are supplied | by trial and error | moderate; transparent if documented |
| 3-D finite element model of ship and dock | full 3-D | yes | possible | possible, with detailed data | by trial and error | high effort; requires detailed model |
| This work | elastic, incl. shear | yes, robust for large gaps | yes | kink and curvature model, or survey | designed automatically from the block strengths (minimax), welding gaps at every candidate position | moderate; open, verified software |

Table: Table 2.1 – Approaches to the determination of block loads and block arrangements.

The table shows the gap. Existing approaches *analyse* a given block arrangement; the spacing itself is chosen by rules of thumb and corrected by trial and error. **No existing tool designs the block spacing from the strengths of the individual blocks, with the reactions of an elastic hull on unilateral, crushable supports and with the welding deformation of the hull as a constraint.** The general literature on structural optimisation treats support positions, but for continuous structures with bilateral supports and smooth responses [@mroz1975; @wang2004], not for a free–free hull resting on compression-only blocks whose initial gaps follow a welded keel line and change with every move of a block.

From the point of view of computer science, the gap is not a missing physical theory but a missing **computational solution**: the algorithms, data structures and software that turn the theory into exact, robust, efficient and verifiable results for the problems that arise in practice. It lies at the intersection of four topics that are usually treated separately: **exact algorithms for unilateral contact** on a sparse elastic model of the hull girder, robust also for the large gaps of welded hulls; **data processing and inverse identification** of the as-built keel shape from survey data; **non-smooth design optimisation** of the block positions, in which every function evaluation is a complete contact analysis and the welding gaps move with the design; and **software engineering** that makes the combination verified, maintainable and usable. The ship engineering topics — the weight distribution, the elastic support model, the welding deformation of block-built hulls and the hull girder and block checks — define the problems these methods must solve. More specifically, what is missing is an implementation that is simultaneously

- **algorithmically exact and robust**, solving the piecewise-linear complementarity problem of the block law without artificial penalty or tolerance parameters that influence the result, with cycle detection and a continuation strategy for the large gaps caused by welding deformation, so that every candidate layout of a design can be analysed;
- **prescriptive**, computing the spacing of a given set of blocks of given strengths so that the largest ratio of reaction to strength is minimised, within spacing limits, the extent of the block line and the floor positions, by a combination of equidistribution, sequential linear programming with a trust region and a discrete local search;
- **constraint-aware**, treating the welding deformation as a constraint of the design by recomputing the gaps at every candidate block position, and diagnosing why a plan cannot be made to fit;
- **efficient**, exploiting the sparse, banded structure of the hull girder model, with a cost that is analysed and measured in terms of the number of analyses;
- **data-aware**, able to build the load from weight items, moments and a weight curve, and to use either a parametric welding model or a keel survey from which interpretable welding parameters are identified by regularised least squares;
- **verified**, with every major component, including the spacing design, tested against closed-form solutions, limiting cases and independent computations in an automated test suite, and with equilibrium checks evaluated at run time;
- **engineered**, with a layered architecture in which the calculation core is independent of its interfaces, and usable through a scriptable command-line interface and a graphical user interface; and
- **transparent and open**, with the derivations, the source code and the input format fully documented so that every number can be traced and reproduced.

The present work develops and documents such an implementation. Its models and algorithms, including their complexity, are derived in Chapter 3, its software engineering is described in Chapter 4, its verification and performance are presented in Chapter 5, and its application to the illustrative vessel — straight and welded, with a designed block spacing as the final result — is reported in Chapter 6.

## 2.15 Summary

Ships rest on keel blocks whose construction, especially the timber cap, governs their stiffness and their resistance to crushing, both in dry dock and in the building dock. Every block has a limited strength, and the blocks should stand under the floors of the hull. Block construction and welding leave the hull with angular distortion at the erection joints and shrinkage curvature within the blocks, which together lift the hull ends and create large initial gaps that depend on the position of each block. The hull girder must carry the concentrated reactions within permissible bending-moment and shear-force limits. The classical rigid-hull method gives a simple linear estimate of the block loads but cannot show where the spacing should change. Beam theory, beams on elastic foundations and the finite element method represent the flexibility of the hull, while complementarity and contact-mechanics theory describe the unilateral behaviour of the blocks, which makes the reactions piecewise-smooth functions of the block positions. Optimisation for design supplies the tools for choosing the positions: linear programming, the epigraph reformulation of minimax objectives and the analogy with fully stressed design, sequential linear programming with trust regions and finite-difference sensitivities, and the equidistribution principle of adaptive mesh methods as a fast global start. Banded sparse factorisation makes each analysis linear in cost with respect to the number of nodes, which makes the hundreds of analyses of a design affordable, and active-set iterations combined with continuation give a robust route to the exact solution of the contact problem. Monotone interpolation and regularised least squares provide the means to use keel surveys. Software engineering principles — requirements, layered architecture with an independent calculation core, dependency injection and automated testing — together with the distinction between verification and validation, govern how such a tool should be built and assessed. The chapter identified the research gap: no existing tool designs the block spacing from the block strengths with the welding deformation as a constraint. Closing it is the subject of the remainder of this thesis.
