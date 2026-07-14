function [mesh, uh, fh, n_refined] = refine_mesh(mesh, uh, fh)
% function [mesh, uh, fh, n_refined] = refine_mesh(mesh, uh, fh)
%   This function refines the mesh with the marked elements
%   from  mesh.mark  and interpolates  uh  and  fh

if (max(mesh.mark)==0)
  % no elements marked, doing nothing
  return
end

%% To make the recursion efficient we copy the mesh and 
% the vectors uh, fh into global variables
% (to avoid several copies)

global adapt_global_mesh adapt_global_uh adapt_global_fh

adapt_global_mesh = mesh;
adapt_global_uh = uh;
adapt_global_fh = fh;

%% We first extend the matrices and vectors to save 
% time in memory allocation.
% Estimation of final number of elements.
nelem_new = adapt_global_mesh.n_elem*2*max(adapt_global_mesh.mark);
adapt_global_mesh.elem_vertices(nelem_new,1) = 0;
adapt_global_mesh.elem_neighbours(nelem_new,1) = 0;
adapt_global_mesh.elem_boundaries(nelem_new,1) = 0;
adapt_global_mesh.mark(nelem_new) = 0;

% Estimation of final number of vertices
nvertices_new = 4*adapt_global_mesh.n_vertices;
adapt_global_mesh.vertex_coordinates(nvertices_new,1) = 0;
adapt_global_uh(nvertices_new,1) = 0;
adapt_global_fh(nvertices_new,1) = 0;

n_refined = 0;

first_marked = find(adapt_global_mesh.mark > 0, 1); 
while (first_marked)
  n_refined = n_refined + refine_element(first_marked);
  first_marked = find(adapt_global_mesh.mark > 0, 1); 
end


% we now clean the unused space
adapt_global_mesh.elem_vertices(adapt_global_mesh.n_elem+1:nelem_new,:) = [];
adapt_global_mesh.elem_neighbours(adapt_global_mesh.n_elem+1:nelem_new,:) = [];
adapt_global_mesh.elem_boundaries(adapt_global_mesh.n_elem+1:nelem_new,:) = [];
adapt_global_mesh.mark(adapt_global_mesh.n_elem+1:nelem_new) = [];

adapt_global_mesh.vertex_coordinates(adapt_global_mesh.n_vertices+1:nvertices_new,:) = [];
adapt_global_uh(adapt_global_mesh.n_vertices+1:nvertices_new,:) = [];
adapt_global_fh(adapt_global_mesh.n_vertices+1:nvertices_new,:) = [];

mesh = adapt_global_mesh;
uh = adapt_global_uh;
fh = adapt_global_fh;