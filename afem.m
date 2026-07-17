%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFEM in 2d for the elliptic problem
%  -div(a grad u) + b . grad u + c u = f   in Omega 
%       u = gD  in Gamma_D
%   du/dn = gN  in Gamma_N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_figures = 1;
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
if (show_figures)
  figure(2)
  clf()
  figure(1)
  triplot(mesh.elem_vertices, ...
    mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
  axis equal
  title('initial mesh');
  pause(1)
end

if (prob_data.initial_global_refinements)
  mesh.mark = prob_data.initial_global_refinements*ones(mesh.n_elem,1);
  [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
  if (show_figures)
    figure(1)
    triplot(mesh.elem_vertices, ...
      mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
    axis equal
    title('initial mesh after global refinements');
    pause
  end
end

%% adaptive strategy
iter_counter = 1;
all_N_vert = []
all_est = [];
all_H1_err = [];
all_L2_err = [];

fprintf('n_elements  n_vertices   H1-error    L2-error    error_est \n');
%        -00000000-  -00000000-  00.0000000  00.0000000  00.0000000
while (1)
  [uh,fh] = assemble_and_solve(prob_data,mesh,uh,fh);
  if (show_figures)
    figure(2)
    trimesh(mesh.elem_vertices, ...
    mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), ...
    uh);
    view(70,12)
    title('approximate solution uh')
    pause(0.1)
  end

  [est, mesh] = estimate(prob_data, adapt, mesh, uh);
  N_vert = mesh.n_vertices;
  h1 = H1_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, prob_data.grad_u_exact);
  l2 = L2_err(mesh.elem_vertices, mesh.vertex_coordinates, uh, prob_data.u_exact);
  fprintf(' %8d    %8d  %10.7f  %10.7f  %10.7f  \n', ...
  mesh.n_elem, N_vert, h1, l2, est);
  all_est(iter_counter) = est;
  all_H1_err(iter_counter) = h1;
  all_L2_err(iter_counter) = l2;
  all_N_vert(iter_counter) = N_vert;
  if ((iter_counter >= adapt.max_iterations) ...
    || (est < adapt.tolerance) ...
    || (N_vert > adapt.max_vertices))
    break;
  else
    mesh = mark_elements(adapt, mesh);
    [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
    if (show_figures)
      figure(1)
      triplot(mesh.elem_vertices, ...
        mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) )
      axis equal
      graph_title=sprintf('mesh at iteration %d, N_p = %d',iter_counter,mesh.n_vertices);
      title(graph_title)
      % print('-dpng',sprintf('pictures/%s-mesh-%03d',prob_data.example,iter_counter))
      % print('-deps',sprintf('pictures/%s-mesh-%03d',prob_data.example,iter_counter))
      pause(0.1)
    end
  end
  iter_counter = iter_counter + 1;
  if(iter_counter > 5) 
    figure(3)
    loglog( ...
      all_N_vert,all_H1_err, ...
      all_N_vert,all_est, ...
      all_N_vert,all_L2_err)
    legend('H1-err', 'estimator', 'L2-err')
    pause(0.1)
  end
end


