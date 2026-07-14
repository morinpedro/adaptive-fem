# AFEM — Adaptive P1 Finite Element Code

Minimal Matlab/Octave code for solving a linear second-order elliptic problem with continuous
piecewise-linear (P1) finite elements on an **adaptively refined** triangular mesh, driven by a
residual a posteriori error estimator. It accompanies Section 6 ("Implementation of AFEM") of the
short course *Adaptive Finite Element Methods for Elliptic PDEs*, and is the adaptive counterpart of
the fixed-mesh code at <https://github.com/morinpedro/fixed-mesh-fem>.

There is no Python code in this repository — the whole adaptive loop (mesh generation, assembly,
estimator, marking, bisection) is Matlab/Octave only.

## The model problem

The code solves

$$-\nabla\cdot(a\,\nabla u) + b\cdot\nabla u + c\,u = f \quad \text{in } \Omega, \qquad
u = g_D \text{ on } \Gamma_D, \qquad a\,\partial_n u = g_N \text{ on } \Gamma_N,$$

with `a > 0`, `b`, `c ≥ 0` piecewise constant. Four examples are built into `init_data.m` (pick one
by editing the `example = '...'` line near the top):

| `example` | Domain | Solution |
|---|---|---|
| `'L-shape'` (default) | `L_shape_dirichlet` | reentrant-corner singular solution `u(r,θ) = r^(2/3) sin(2θ/3)` (`f = 0`, `g_D = u`) — uniform refinement is suboptimal here and adaptivity pays off |
| `'square-Dirichlet'` | `square_all_dirichlet` | smooth Gaussian bump `u(x) = exp(-10‖x‖²)`, fully Dirichlet |
| `'square-mixed'` | `square_mixed` | `f = 1`, homogeneous data, mixed Dirichlet/Neumann boundary |
| `'big-square-Dirichlet'` | `big_square_all_dirichlet` | the same Gaussian bump on a larger square domain |

## The adaptive loop

The driver `afem.m` runs the standard loop, one function per module, passing `prob_data`, `adapt`,
`mesh`, `uh`, `fh` explicitly between them as function arguments:

```
SOLVE  ->  ESTIMATE  ->  MARK  ->  REFINE
```

- **SOLVE** (`[uh,fh] = assemble_and_solve(prob_data,mesh,uh,fh)`): P1 assembly on the current mesh —
  element stiffness and mass matrices via the reference map, midpoint quadrature for the load
  vector, Simpson's rule on Neumann segments, Dirichlet rows replaced by the identity.
- **ESTIMATE** (`[global_est, mesh] = estimate(prob_data, adapt, mesh, uh)`): the residual indicator
  `η_T² = C₁ h_T² ‖r_T‖²_{L²(T)} + C₂ h_T ‖J_T‖²_{L²(∂T)}`, with interior residual
  `r_T = b·∇U + cU − f` and normal-flux jumps `J_T` across interior/Neumann sides.
- **MARK** (`mesh = mark_elements(adapt, mesh)`): three strategies — `GR` (global/uniform), `MS`
  (maximum strategy), and `Doerfler` (bulk/Dörfler strategy, the shipped default).
- **REFINE** (`[mesh,uh,fh,n_refined] = refine_mesh(mesh,uh,fh)`, using `refine_element.m`):
  recursive newest-vertex bisection with completion, keeping the mesh conforming. These two
  functions are the one deliberate exception to the argument-passing style: to avoid Matlab
  copying the whole mesh struct on every recursive call, they communicate through the globals
  `adapt_global_mesh`, `adapt_global_uh`, `adapt_global_fh` internally; `refine_mesh` still takes
  and returns `mesh`, `uh`, `fh` normally, so this is invisible to the caller.

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

Four ready-made meshes are included: `L_shape_dirichlet/`, `square_all_dirichlet/`,
`square_mixed/`, and `big_square_all_dirichlet/`. To start the adaptive loop from a finer (uniform)
L-shape mesh, generate one with `gen_mesh_L_shape` (see below) instead of relying on the default
6-triangle mesh.

## Requirements

**Matlab or GNU Octave** for the `.m` files (Octave ≥ 5 recommended; it is free). No other
dependencies.

## Files

| File | Description |
|------|-------------|
| `afem.m` | Driver: loads the mesh and problem data, runs the SOLVE–ESTIMATE–MARK–REFINE loop, and plots mesh and solution. |
| `init_data.m` | `[prob_data, adapt] = init_data()` — problem data (`a, b, c, f, gD, gN`, exact solution) for whichever of the four built-in `example`s is selected, plus the adaptivity parameters (strategy, tolerance, marking constants, `n_refine`), which apply across all examples. |
| `assemble_and_solve.m` | SOLVE: fixed-mesh P1 assembly on the current mesh. |
| `estimate.m` | ESTIMATE: residual-type a posteriori error indicators. |
| `mark_elements.m` | MARK: the `GR`, `MS`, and `Doerfler` strategies. |
| `refine_mesh.m`, `refine_element.m` | REFINE: recursive newest-vertex bisection with completion (uses Matlab globals internally for performance — see above). |
| `get_dirichlet_neumann.m` | Reconstructs the Dirichlet vertex list and Neumann segments from `elem_boundaries`. |
| `H1_err.m`, `L2_err.m` | Errors against a known solution (midpoint quadrature); called from `afem.m` using the exact solution of the active example. |
| `gen_mesh_L_shape.m` | Generate a uniform L-shape starting mesh (any refinement level) in the element-based format. |
| `L_shape_dirichlet/`, `square_all_dirichlet/`, `square_mixed/`, `big_square_all_dirichlet/` | Initial-mesh data folders. |

There are no `u_ex3.m`/`grdu_ex3.m` files: each example's exact solution and gradient are now
defined inline, as anonymous functions, inside its `case` block in `init_data.m`.

## Quick start (Octave/Matlab)

From this folder, start Octave (or Matlab) and run:

```matlab
afem
```

The driver calls `init_data.m` (which selects the `example` set near the top of that file — default
`'L-shape'`), loads the corresponding initial mesh, applies `prob_data.initial_global_refinements`
uniform refinements if requested, and then iterates the adaptive loop, pausing at each step so you
can watch the mesh concentrate near the reentrant corner (or wherever the estimator is largest for
the chosen example). Change the marking strategy (`adapt.strategy = 'MS' | 'Doerfler' | 'GR'`) and
its parameters in `init_data.m`.

### Finer starting meshes

Instead of the default 6-triangle L-shape mesh you can generate a uniform L-shape mesh at any level
and let `afem` load it (the generator writes into `L_shape_dirichlet/` by default):

```matlab
gen_mesh_L_shape(4)     % 96 triangles, 65 vertices  ->  L_shape_dirichlet/
afem                    % runs the adaptive loop from that finer starting mesh
```

`N = 1` reproduces the shipped mesh exactly; larger `N` refines uniformly (`h = 1/N`). Pass a second
argument to write to a different folder, e.g. `gen_mesh_L_shape(8, 'L_shape_fine')`, and set
`prob_data.domain` in `init_data.m` accordingly. The generated triangulation uses the same
newest-vertex labeling as the shipped mesh, so recursive bisection stays conforming.

## Reproducing a convergence study

The shipped defaults (`adapt.Doerfler_theta = 0.3`, `adapt.max_iterations = 10`) are tuned for a
quick, interactive run: `theta = 0.3` is a conservative bulk-marking target, so relatively few
elements get marked each iteration and, combined with the 10-iteration cap, the mesh stays small
(on the order of a few hundred vertices). That is fine for watching the loop converge in the plot
window, but it will *not* reach the mesh sizes needed to see the asymptotic convergence rate
clearly.

To reproduce a multi-order-of-magnitude convergence study like the one in Section 6 of the notes,
edit `init_data.m` and raise the marking target and iteration budget, e.g.

```matlab
adapt.Doerfler_theta = 0.5;
adapt.max_iterations = 25;
adapt.tolerance = 1e-7;
```

then run `afem` and read off `n_elements`, `n_vertices`, and the printed `H1-error` at each
iteration. On the `'L-shape'` example this recovers the optimal energy rate `Np^{-1/2}`, in contrast
to the `Np^{-1/3}` rate of uniform refinement (`adapt.strategy = 'GR'`).

## Notes

- The problem-data functions in `init_data.m` are anonymous functions (e.g.
  `prob_data.f = @(x) 0;`), so the code runs unchanged in both GNU Octave and current MATLAB.
- The mesh folder names use underscores (`L_shape_dirichlet`); make sure `prob_data.domain` in
  `init_data.m` matches the folder you want to load.
- `n_refine = 2` bisects each marked element twice, adding an interior node — the interior-node
  property the convergence theory relies on.
- The Neumann jump term in `estimate.m` is written for homogeneous Neumann data (`gN = 0`); to run a
  problem with nonzero `gN`, extend the `- 0` term in the Neumann branch to `- gN(midpoint)`.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
