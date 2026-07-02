# AFEM — Adaptive P1 Finite Element Code

Minimal Matlab/Octave code for solving a linear second-order elliptic problem with continuous
piecewise-linear (P1) finite elements on an **adaptively refined** triangular mesh, driven by a
residual a posteriori error estimator. A Python script reproduces the method and runs the
convergence study. It accompanies Section 6 ("Implementation of an Adaptive Method") of the short
course *Adaptive Finite Element Methods for Elliptic PDEs*, and is the adaptive counterpart of the
fixed-mesh code at <https://github.com/morinpedro/fixed-mesh-fem>.

## The model problem

The code solves

$$-\nabla\cdot(a\,\nabla u) + b\cdot\nabla u + c\,u = f \quad \text{in } \Omega, \qquad
u = g_D \text{ on } \Gamma_D, \qquad a\,\partial_n u = g_N \text{ on } \Gamma_N,$$

with `a > 0`, `b`, `c ≥ 0` piecewise constant. The default example is the classical reentrant-corner
problem on the L-shaped domain, with the singular harmonic solution
`u(r,θ) = r^(2/3) sin(2θ/3)` (so `f = 0` and `g_D = u`), for which uniform refinement is suboptimal
and adaptivity pays off.

## The adaptive loop

The driver `afem.m` runs the standard loop, one function per module:

```
SOLVE  ->  ESTIMATE  ->  MARK  ->  REFINE
```

- **SOLVE** (`assemble_and_solve.m`): P1 assembly on the current mesh — element stiffness and mass
  matrices via the reference map, midpoint quadrature for the load vector, Simpson's rule on Neumann
  segments, Dirichlet rows replaced by the identity.
- **ESTIMATE** (`estimate.m`): the residual indicator
  `η_T² = C₁ h_T² ‖r_T‖²_{L²(T)} + C₂ h_T ‖J_T‖²_{L²(∂T)}`, with interior residual
  `r_T = b·∇U + cU − f` and normal-flux jumps `J_T` across interior/Neumann sides.
- **MARK** (`mark_elements.m`): three strategies — `GR` (global/uniform), `MS` (maximum strategy),
  and `GERS` (Dörfler / guaranteed error reduction).
- **REFINE** (`refine_mesh.m`, `refine_element.m`): recursive newest-vertex bisection with
  completion, keeping the mesh conforming.

## Mesh representation

Each domain lives in its own folder holding four plain-text files, all indexed **by element** — the
richer connectivity (neighbours and per-side boundary flags) is what makes local refinement and the
jump term `O(1)` per element:

| File | Contents |
|------|----------|
| `vertex_coordinates.txt` | one `x y` per vertex |
| `elem_vertices.txt` | three global vertex indices `v1 v2 v3` per triangle (side `i` is opposite vertex `i`; side 3, the edge `v1 v2`, is the refinement edge) |
| `elem_neighbours.txt` | the triangle across each side, or `0` on the boundary |
| `elem_boundaries.txt` | per-side flag: `0` interior, `>0` Dirichlet, `<0` Neumann |

Three ready-made meshes are included: `L_shape_dirichlet/`, `square_all_dirichlet/`, and
`square_mixed/`. To start the adaptive loop from a finer (uniform) L-shape mesh, generate one with
`gen_mesh_L_shape` (see below) instead of relying on the default 6-triangle mesh.

## Requirements

- **Matlab or GNU Octave** for the `.m` files (Octave ≥ 5 recommended; it is free).
- **Python 3** with **NumPy** and **SciPy** for `adaptive_Lshape.py`:
  ```bash
  pip install numpy scipy
  ```

## Files

| File | Description |
|------|-------------|
| `afem.m` | Driver: runs the SOLVE–ESTIMATE–MARK–REFINE loop and plots mesh and solution. |
| `init_data.m` | Problem data (`a, b, c, f, gD, gN`) and adaptivity parameters (strategy, tolerance, marking constants, `n_refine`). |
| `assemble_and_solve.m` | SOLVE: fixed-mesh P1 assembly on the current mesh. |
| `estimate.m` | ESTIMATE: residual-type a posteriori error indicators. |
| `mark_elements.m` | MARK: the `GR`, `MS`, and `GERS` strategies. |
| `refine_mesh.m`, `refine_element.m` | REFINE: recursive newest-vertex bisection with completion. |
| `get_dirichlet_neumann.m` | Reconstructs the Dirichlet vertex list and Neumann segments from `elem_boundaries`. |
| `H1_err.m`, `L2_err.m` | Errors against a known solution (midpoint quadrature). |
| `u_ex3.m`, `grdu_ex3.m` | The exact reentrant-corner solution and its gradient. |
| `gen_mesh_L_shape.m` | Generate a uniform L-shape starting mesh (any refinement level) in the element-based format. |
| `gen_mesh_L_shape.py` | Same generator in Python (identical output). |
| `adaptive_Lshape.py` | Python reproduction: runs the full loop and prints the convergence table. |
| `L_shape_dirichlet/`, `square_all_dirichlet/`, `square_mixed/` | Initial-mesh data folders. |

## Quick start (Octave/Matlab)

From this folder, start Octave (or Matlab) and run:

```matlab
afem
```

The driver reads `init_data.m`, loads the initial mesh named by the `domain` variable there
(default `L_shape_dirichlet`), applies `global_refinements` uniform refinements, and then iterates
the adaptive loop, pausing at each step so you can watch the mesh concentrate near the reentrant
corner. Change the marking strategy (`adapt.strategy = 'MS' | 'GERS' | 'GR'`) and its parameters in
`init_data.m`.

### Finer starting meshes

Instead of the default 6-triangle mesh you can generate a uniform L-shape mesh at any level and let
`afem` load it (the generator writes into `L_shape_dirichlet/` by default):

```matlab
gen_mesh_L_shape(4)     % 96 triangles, 65 vertices  ->  L_shape_dirichlet/
afem                    % runs the adaptive loop from that finer starting mesh
```

`N = 1` reproduces the shipped mesh exactly; larger `N` refines uniformly (`h = 1/N`). Pass a second
argument to write to a different folder, e.g. `gen_mesh_L_shape(8, 'L_shape_fine')`, and set
`domain` in `init_data.m` accordingly. A Python version with identical output is provided for those
who prefer it: `python3 gen_mesh_L_shape.py 4`. The generated triangulation uses the same
newest-vertex labeling as the shipped mesh, so recursive bisection stays conforming.

## Convergence study (Python)

`adaptive_Lshape.py` reproduces the loop on the L-shape and prints the error table and observed rate
automatically:

```bash
python3 adaptive_Lshape.py
```

### Expected output

The energy error decays like `N_p^{-1/2}` — the **optimal** rate — even though `u ∈ H^{1+2/3-ε}`
only, in contrast to the suboptimal `N_p^{-1/3}` of uniform refinement:

```
     Np       H1err     rate
     84   1.0557e-01      ---
    230   6.2106e-02    0.53
    793   3.2974e-02    0.51
   2256   1.9189e-02    0.52
   6584   1.1137e-02    0.51
  18762   6.5434e-03    0.51
  53615   3.8587e-03    0.50
```

## Notes

- The problem-data functions in `init_data.m` (and `afem.m`) are anonymous functions
  (e.g. `prob_data.f = @(x) 0;`), so the code runs unchanged in both GNU Octave and current MATLAB.
- The mesh folder names use underscores (`L_shape_dirichlet`); make sure the `domain` variable in
  `init_data.m` matches the folder you want to load.
- `n_refine = 2` bisects each marked element twice, adding an interior node — the interior-node
  property the convergence theory relies on.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
