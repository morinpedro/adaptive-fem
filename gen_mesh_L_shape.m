function gen_mesh_L_shape(N, folder)
% gen_mesh_L_shape(N, folder)
%   Generate a uniform starting mesh of the L-shaped domain
%       Omega = (-1,1)^2  \  ( (0,1] x [-1,0) )
%   in the element-based format read by the adaptive code (afem.m):
%
%       vertex_coordinates.txt   one  "x y"  per vertex
%       elem_vertices.txt        three vertex indices  v1 v2 v3  per triangle
%       elem_neighbours.txt      triangle across each side  (0 on the boundary)
%       elem_boundaries.txt      per-side flag: 0 interior, >0 Dirichlet, <0 Neumann
%
%   Every unit sub-square is split by its bottom-left--to--top-right diagonal into
%   two triangles laid out for newest-vertex bisection: side 3 (the edge v1 v2) is
%   the refinement edge and v3 is the newest vertex.  This is the same consistent
%   labeling as the shipped 6-triangle mesh (reproduced exactly for N = 1), so the
%   recursive bisection in refine_element.m stays conforming and terminates.
%   The boundary is treated as fully Dirichlet (flag +1), matching L_shape_dirichlet.
%
%   Inputs
%     N       subdivisions per unit length (h = 1/N).  N = 1 gives the shipped mesh.
%     folder  output directory (created if needed).  Default 'L_shape_dirichlet',
%             i.e. the directory afem.m loads by default, so afterwards you can run
%             afem with a finer starting mesh and no code change.
%
%   Examples
%     gen_mesh_L_shape(4)                 % writes L_shape_dirichlet/  (65 vertices)
%     gen_mesh_L_shape(8, 'L_shape_fine') % writes L_shape_fine/       (225 vertices)

if (nargin < 2)
  folder = 'L_shape_dirichlet';
end

h   = 1 / N;
npt = 2*N + 1;                       % grid points per direction

% --- pass 1: keep cells whose centre is in Omega; mark the grid points they use ---
used = false(npt, npt);             % used(i+1, j+1) for grid indices i,j = 0..2N
kept = zeros(0, 2);
for j = 2*N-1 : -1 : 0              % top row of cells first
  for i = 0 : 2*N-1
    cx = -1 + (i+0.5)*h;  cy = -1 + (j+0.5)*h;
    if (cx > 0 && cy < 0)          % the missing lower-right quadrant
      continue;
    end
    kept(end+1, :) = [i j];
    used(i+1, j+1) = true;  used(i+2, j+1) = true;
    used(i+1, j+2) = true;  used(i+2, j+2) = true;
  end
end

% --- pass 2: number used vertices, scanning top -> bottom, left -> right ---
gid    = zeros(npt, npt);
nv     = 0;
coords = zeros(0, 2);
for j = 2*N : -1 : 0
  for i = 0 : 2*N
    if (used(i+1, j+1))
      nv = nv + 1;
      gid(i+1, j+1) = nv;
      coords(nv, :) = [-1 + i*h, -1 + j*h];
    end
  end
end

% --- pass 3: two triangles per kept cell (newest-vertex layout) ---
nc   = size(kept, 1);
nt   = 2 * nc;
elem = zeros(nt, 3);
for k = 1 : nc
  i = kept(k, 1);  j = kept(k, 2);
  bl = gid(i+1, j+1);  br = gid(i+2, j+1);
  tl = gid(i+1, j+2);  tr = gid(i+2, j+2);
  elem(2*k-1, :) = [bl tr tl];      % A: refinement edge bl-tr, newest vertex tl
  elem(2*k,   :) = [tr bl br];      % B: refinement edge tr-bl, newest vertex br
end

% --- neighbours and boundary flags via a sorted edge list ---
% side i is opposite vertex i:  side1 = (v2,v3), side2 = (v3,v1), side3 = (v1,v2)
sl = [2 3; 3 1; 1 2];
E  = zeros(3*nt, 4);                 % [minv maxv tri side]
r  = 0;
for t = 1 : nt
  for s = 1 : 3
    va = elem(t, sl(s,1));  vb = elem(t, sl(s,2));
    r = r + 1;
    E(r, :) = [min(va,vb) max(va,vb) t s];
  end
end
[~, ord] = sortrows(E(:, 1:2));      % group identical edges
E = E(ord, :);

ne = zeros(nt, 3);                   % 0 on the boundary
bd = zeros(nt, 3);                   % 0 interior, +1 Dirichlet
r  = 1;  R = 3*nt;
while (r <= R)
  if (r < R && E(r,1) == E(r+1,1) && E(r,2) == E(r+1,2))
    t1 = E(r,3);  s1 = E(r,4);  t2 = E(r+1,3);  s2 = E(r+1,4);
    ne(t1, s1) = t2;  ne(t2, s2) = t1;
    r = r + 2;
  else
    t1 = E(r,3);  s1 = E(r,4);
    bd(t1, s1) = 1;                  % fully Dirichlet boundary
    r = r + 1;
  end
end

% --- write the four text files ---
if (exist(folder, 'dir') == 0)
  mkdir(folder);
end
dlmwrite(fullfile(folder, 'vertex_coordinates.txt'), coords, 'delimiter', ' ', 'precision', '%g');
dlmwrite(fullfile(folder, 'elem_vertices.txt'),      elem,   'delimiter', ' ', 'precision', '%d');
dlmwrite(fullfile(folder, 'elem_neighbours.txt'),    ne,     'delimiter', ' ', 'precision', '%d');
dlmwrite(fullfile(folder, 'elem_boundaries.txt'),    bd,     'delimiter', ' ', 'precision', '%d');

fprintf(1, 'wrote %d vertices and %d triangles to %s (N = %d, h = 1/%d)\n', nv, nt, folder, N, N);
end
