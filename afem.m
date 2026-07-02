%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFEM in 2d for the elliptic problem
%  -div(a grad u) + b . grad u + c u = f   in Omega 
%       u = gD  in Gamma_D
%   du/dn = gN  in Gamma_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we first read all the initial data from the file init_data.m
init_data
% at this point we have filled structures
%       prob_data  and  adapt
% that are also global variables, and thus
% visible in all functions

% we now define three more global variables for the mesh to save memory
global  mesh  uh  fh

%% initialize mesh and finite element functions
% read the mesh from 'domain'
mesh.elem_vertices      = load([domain '/elem_vertices.txt']); 
mesh.elem_neighbours    = load([domain '/elem_neighbours.txt']);
mesh.elem_boundaries    = load([domain '/elem_boundaries.txt']);
mesh.vertex_coordinates = load([domain '/vertex_coordinates.txt']);
mesh.n_elem     = size(mesh.elem_vertices, 1);
mesh.n_vertices = size(mesh.vertex_coordinates, 1);
uh = zeros(mesh.n_vertices, 1);
fh = zeros(mesh.n_vertices, 1);

grd_u_exact = @(x) grdu_ex3(x);
    u_exact = @(x) u_ex3(x);
% plot the mesh
figure(2)
clf()
figure(1)
triplot(mesh.elem_vertices, ...
	mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
axis equal
title('initial mesh');
pause

if (global_refinements)
  mesh.mark = global_refinements*ones(mesh.n_elem,1);
  refine_mesh;
  figure(1)
  triplot(mesh.elem_vertices, ...
	  mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
  axis equal
  title('initial mesh after global refinements');
  pause
end

%% adaptive strategy
iter_counter = 1;
all_est = [];
all_err = [];
while (1)
  assemble_and_solve;
  figure(2)
  trimesh(mesh.elem_vertices, ...
	  mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
	  uh);
  title('approximate solution uh')
  pause

  est = estimate(prob_data, adapt);
  fprintf('n_elements: %5d  H1-error: %.6f  L2-error: %.6f  error estimation: %.6f\n',...
	 mesh.n_elem, ...
	 H1_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, grd_u_exact), ...
	 L2_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, u_exact), ...
	 est);

  if ((iter_counter >= adapt.max_iterations) | (est < adapt.tolerance))
    break;
  else
    mark_elements(adapt);
    refine_mesh;
    figure(1)
    triplot(mesh.elem_vertices, ...
	    mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
    axis equal
    title('new mesh')
    pause
  end
  iter_counter = iter_counter + 1;
  
end


