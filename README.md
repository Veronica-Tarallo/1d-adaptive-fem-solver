# 1D Adaptive Finite Element Solver (AFEM) for Elliptic Boundary Value Problems

A from-scratch MATLAB implementation of a **1D Finite Element Method (FEM)** and **Adaptive Finite Element Method (AFEM)** solver for second-order elliptic boundary value problems with homogeneous Dirichlet boundary conditions:

$$-u''(x) = f(x) \quad \text{in } \Omega = (0, 1), \qquad u(0) = u(1) = 0$$

The solver supports both **linear ($P_1$)** and **quadratic ($P_2$)** Lagrangian finite elements and compares uniform mesh refinement against adaptive mesh refinement driven by an *a posteriori* error estimator.

## Mathematical Formulation & Adaptive Loop

To test convergence across different regularity regimes, the solver evaluates the benchmark analytical solution $u(x) = x^\alpha (1 - x)$ for:
* **Regular regime ($\alpha = 10$):** Smooth solution $u \in H^3(\Omega)$.
* **Singular / low-regularity regime ($\alpha = 4/3$):** Solution $u \in H^1(\Omega) \setminus H^2(\Omega)$ exhibiting a derivative singularity at $x = 0$, where uniform mesh refinement suffers from suboptimal convergence rates.

The adaptive algorithm executes the standard four-step loop:
1. **SOLVE:** Assembles the global stiffness matrix and right-hand side via element-wise Midpoint ($P_1$) or Simpson ($P_2$) quadrature rules and solves the resulting sparse linear system.
2. **ESTIMATE:** Computes a local residual-and-jump *a posteriori* error indicator $\eta_K^2$ on each element $K$, combining interior residual norms $\|f + u_h''\|_{0,K}^2$ (integrated via 3-point Simpson quadrature) with inter-element gradient jumps $[u_h']$ across interior nodes.
3. **MARK:** Marks the top **25%** of elements with the largest local error indicators $\eta_K^2$.
4. **REFINE:** Performs bisection refinement by inserting midpoint nodes into all marked elements.

## Key Numerical Results

* **Uniform Refinement:** Achieves optimal $H^1$-seminorm convergence rates $\mathcal{O}(h)$ for $P_1$ and $\mathcal{O}(h^2)$ for $P_2$ when $\alpha = 10$, but degrades to suboptimal rates for both $P_1$ and $P_2$ when $\alpha = 4/3$ due to the singularity at $x = 0$.
* **Adaptive Refinement (AFEM):** Dynamically concentrates degrees of freedom near $x = 0$, recovering the optimal asymptotic convergence rates $\mathcal{O}(N^{-1})$ for $P_1$ and $\mathcal{O}(N^{-2})$ for $P_2$ with respect to the number of Degrees of Freedom ($N$) even in the low-regularity case $\alpha = 4/3$.

## Repository Structure

### Main Scripts & Visualization
* `Main_FEM1D.m` - Main driver script running both uniform and adaptive FEM loops for configurable regularity $\alpha \in \{4/3, 10\}$ and polynomial degree $k \in \{1, 2\}$, saving simulation outputs to `.mat` files.
* `main_plot_confronto.m` - Generates the $2 \times 2$ log-log comparison grid evaluating $H^1$-seminorm error vs. Degrees of Freedom (Uniform vs. Adaptive) against theoretical $\mathcal{O}(N^{-1})$ and $\mathcal{O}(N^{-2})$ slopes.
* `grafici_rispetto_h.m` - Generates log-log convergence plots comparing $P_1$ vs. $P_2$ elements with respect to mesh size $h$ on uniform grids.
* `dati_confronto_*.mat` - Pre-computed benchmark datasets for all 4 parameter combinations ($\alpha \in \{4/3, 10\}$, $k \in \{1, 2\}$).

### Core FEM & Error Estimation Routines
* `FEM1D_byElements.m` - Element-by-element assembly of the sparse stiffness matrix and load vector for $P_1$ and $P_2$ Lagrangian elements.
* `stimatore_errore_1D.m` - Computes the element-wise *a posteriori* residual-and-jump error estimator $\eta_K^2$ (`ESTIMATE` phase).
* `H1s_1D_error.m` - Evaluates the global approximation error in the $H^1$-seminorm $|u - u_h|_{1,\Omega}$.
* `quadrature_formula_1D.m` - Provides quadrature nodes and weights on the reference interval $[0,1]$ (`gauss2`, `gauss3`, and `simpson`).

### Basis Functions & Problem Definition
* `phiP1.m`, `phiP2.m` - Linear ($P_1$) and quadratic ($P_2$) Lagrangian shape functions on $[0,1]$.
* `gradphiP1.m`, `gradphiP2.m` - Derivatives of the $P_1$ and $P_2$ reference basis functions.
* `exactu.m`, `gradexactu.m`, `effe.m` - Exact analytical solution $u(x)$, its first derivative $u'(x)$, and the source term $f(x)$.

## How to Run

1. Clone or download this repository and open MATLAB in the project directory.
2. Run `main_plot_confronto.m` and `grafici_rispetto_h.m` directly using the included `.mat` benchmark datasets, or re-run `Main_FEM1D.m` across all four parameter combinations (`alpha = 4/3` or `10`, and `kk = 1` or `2`) to regenerate the simulation outputs from scratch.
