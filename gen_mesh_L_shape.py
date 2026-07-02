#!/usr/bin/env python3
"""
gen_mesh_L_shape.py

Generate a uniform starting mesh of the L-shaped domain

    Omega = (-1,1)^2  \\  ( (0,1] x [-1,0) )

in the *element-based* format read by the adaptive code (afem.m):

    vertex_coordinates.txt   one  "x y"  per vertex
    elem_vertices.txt        three global vertex indices  v1 v2 v3  per triangle
    elem_neighbours.txt      triangle across each side  (0 on the boundary)
    elem_boundaries.txt      per-side flag: 0 interior, >0 Dirichlet, <0 Neumann

Every unit sub-square is split by its bottom-left--to--top-right diagonal into two
triangles laid out for newest-vertex bisection: side 3 (the edge v1 v2) is the
refinement edge and v3 is the newest vertex.  This is the same consistent labeling
as the shipped 6-triangle mesh (which this script reproduces exactly for N = 1), so
the recursive bisection in refine_element.m stays conforming and terminates.

The boundary is treated as fully Dirichlet (flag +1), matching L_shape_dirichlet.

Usage
-----
    python3 gen_mesh_L_shape.py N [folder]

N       subdivisions per unit length (h = 1/N).  N = 1 gives the shipped mesh.
folder  output directory (created if needed).  Default: L_shape_dirichlet, i.e.
        the directory afem.m loads by default, so afterwards you can simply run
        afem with a finer starting mesh and no code change.

Examples
--------
    python3 gen_mesh_L_shape.py 4                 # writes L_shape_dirichlet/  (65 vertices)
    python3 gen_mesh_L_shape.py 8 L_shape_fine    # writes L_shape_fine/       (225 vertices)
"""
import os
import sys


def gen_mesh_L_shape(N):
    """Return (coords, elem, neigh, bdry) with 1-based vertex indices."""
    h = 1.0 / N
    npt = 2 * N + 1
    idx = lambda i, j: i + j * npt  # flat key for a grid point (i = x-index, j = y-index)

    # --- pass 1: keep the cells whose centre lies in Omega; mark the grid points used ---
    used = {}
    kept = []
    for j in range(2 * N - 1, -1, -1):          # top row of cells first
        for i in range(0, 2 * N):
            cx, cy = -1.0 + (i + 0.5) * h, -1.0 + (j + 0.5) * h
            if cx > 0 and cy < 0:               # the missing lower-right quadrant
                continue
            kept.append((i, j))
            for (a, b) in ((i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)):
                used[idx(a, b)] = (a, b)

    # --- pass 2: number the used vertices, scanning top -> bottom, left -> right ---
    gid = {}
    coords = []
    for j in range(2 * N, -1, -1):
        for i in range(0, 2 * N + 1):
            if idx(i, j) in used:
                gid[idx(i, j)] = len(coords) + 1
                coords.append((-1.0 + i * h, -1.0 + j * h))

    # --- pass 3: two triangles per kept cell (newest-vertex layout) ---
    elem = []
    for (i, j) in kept:
        bl = gid[idx(i, j)]; br = gid[idx(i + 1, j)]
        tl = gid[idx(i, j + 1)]; tr = gid[idx(i + 1, j + 1)]
        elem.append((bl, tr, tl))               # A: refinement edge bl-tr, newest vertex tl
        elem.append((tr, bl, br))               # B: refinement edge tr-bl, newest vertex br

    # --- neighbours and boundary flags via an edge -> (triangle, side) map ---
    # side i is opposite vertex i:  side1 = (v2,v3), side2 = (v3,v1), side3 = (v1,v2)
    side_local = ((1, 2), (2, 0), (0, 1))
    edge_map = {}
    for t, tri in enumerate(elem):
        for s, (pa, pb) in enumerate(side_local):
            va, vb = tri[pa], tri[pb]
            edge_map.setdefault((min(va, vb), max(va, vb)), []).append((t, s))

    nt = len(elem)
    neigh = [[0, 0, 0] for _ in range(nt)]
    bdry = [[0, 0, 0] for _ in range(nt)]
    for key, lst in edge_map.items():
        if len(lst) == 2:
            (t1, s1), (t2, s2) = lst
            neigh[t1][s1] = t2 + 1
            neigh[t2][s2] = t1 + 1
        elif len(lst) == 1:
            (t1, s1), = lst
            bdry[t1][s1] = 1                     # fully Dirichlet boundary
        else:
            raise RuntimeError("non-manifold edge %s" % (key,))
    return coords, elem, neigh, bdry


def write_mesh(folder, coords, elem, neigh, bdry):
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "vertex_coordinates.txt"), "w") as f:
        for (x, y) in coords:
            f.write("%g %g\n" % (x, y))
    with open(os.path.join(folder, "elem_vertices.txt"), "w") as f:
        for t in elem:
            f.write("%d %d %d\n" % t)
    with open(os.path.join(folder, "elem_neighbours.txt"), "w") as f:
        for r in neigh:
            f.write("%d %d %d\n" % tuple(r))
    with open(os.path.join(folder, "elem_boundaries.txt"), "w") as f:
        for r in bdry:
            f.write("%d %d %d\n" % tuple(r))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python3 gen_mesh_L_shape.py N [folder]")
    N = int(sys.argv[1])
    folder = sys.argv[2] if len(sys.argv) > 2 else "L_shape_dirichlet"
    coords, elem, neigh, bdry = gen_mesh_L_shape(N)
    write_mesh(folder, coords, elem, neigh, bdry)
    print("wrote %d vertices and %d triangles to '%s/' (N = %d, h = 1/%d)"
          % (len(coords), len(elem), folder, N, N))
