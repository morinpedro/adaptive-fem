%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFEM in 2d for the elliptic problem
%  -div(a grad u) + b . grad u + c u = f   in Omega 
%       u = gD  in Gamma_D
%   du/dn = gN  in Gamma_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we first read all the initial data from the file init_data.m
[prob_data, adapt] = init_data();
% at this point we have filled structures
%       prob_data  and  adapt

%% initialize mesh and finite element functions
% read the mesh from 'domain'
mesh.elem_vertices      = load([prob_data.domain '/elem_vertices.txt']); 
mesh.elem_neighbours    = load([prob_data.domain '/elem_neighbours.txt']);
mesh.elem_boundaries    = load([prob_data.domain '/elem_boundaries.txt']);
mesh.vertex_coordinates = load([prob_data.domain '/vertex_coordinates.txt']);
mesh.n_elem     = size(mesh.elem_vertices, 1);
mesh.n_vertices = size(mesh.vertex_coordinates, 1);
uh = zeros(mesh.n_vertices, 1);
fh = zeros(mesh.n_vertices, 1);

% plot the mesh
figure(2)
clf()
figure(1)
triplot(mesh.elem_vertices, ...
	mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
axis equal
title('initial mesh');
pause

if (prob_data.initial_global_refinements)
  mesh.mark = prob_data.initial_global_refinements*ones(mesh.n_elem,1);
  [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
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

fprintf('n_elements  n_vertices   H1-error    L2-error    error_est \n');
%        -00000000-  -00000000-  00.0000000  00.0000000  00.0000000
while (1)
  [uh,fh] = assemble_and_solve(prob_data,mesh,uh,fh);
  figure(2)
  trimesh(mesh.elem_vertices, ...
  mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
  uh);
  title('approximate solution uh')
  pause(0.1)
  
  [est, mesh] = estimate(prob_data, adapt, mesh, uh);
  fprintf(' %8d    %8d  %10.7f  %10.7f  %10.7f  \n',...
    mesh.n_elem, mesh.n_vertices, ...
    H1_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, prob_data.grad_u_exact), ...
    L2_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, prob_data.u_exact), ...
    est);

  if ((iter_counter >= adapt.max_iterations) || (est < adapt.tolerance))
    break;
  else
    mesh = mark_elements(adapt, mesh);
    [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
    figure(1)
    triplot(mesh.elem_vertices, ...
	    mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
    axis equal
    title('new mesh')
    pause(0.1)
  end
  iter_counter = iter_counter + 1;
  
end


