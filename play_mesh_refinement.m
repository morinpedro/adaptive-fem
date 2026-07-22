% play with a mesh

[prob_data, adapt] = init_data();
mesh.elem_vertices = load([prob_data.domain '/elem_vertices.txt']);
mesh.elem_neighbours = load([prob_data.domain '/elem_neighbours.txt']);
mesh.elem_boundaries = load([prob_data.domain '/elem_boundaries.txt']);
mesh.vertex_coordinates = load([prob_data.domain '/vertex_coordinates.txt']);
mesh.n_elem = size(mesh.elem_vertices, 1);
mesh.n_vertices = size(mesh.vertex_coordinates, 1);
uh = zeros(mesh.n_vertices, 1);
fh = zeros(mesh.n_vertices, 1);

% plot the mesh
figure(4)
triplot(mesh.elem_vertices, ...
mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) ,'LineWidth',2 )
axis equal
title('initial mesh','FontSize',20);
% pause

if (prob_data.initial_global_refinements)
  mesh.mark = prob_data.initial_global_refinements*ones(mesh.n_elem,1);
  [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
  triplot(mesh.elem_vertices, ...
  mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2), 'LineWidth', 2)
  % plot_element_numbers(mesh)
  axis equal
  title('initial mesh after global refinements (hit ENTER)','FontSize',20);
  % pause
end
oldmesh = mesh;
while 1
  title('Click on elements (Enter or right-click to stop)', 'FontSize',20);

  el_to_refine = click_elements(mesh); 
  mesh.mark = zeros(mesh.n_vertices,1);
  mesh.mark(el_to_refine) = 1;
  [mesh,uh,fh] = refine_mesh(mesh,uh,fh);
  triplot(mesh.elem_vertices, ...
  mesh.vertex_coordinates(:,1), mesh.vertex_coordinates(:,2) ,'r', 'LineWidth', 2)
  axis equal
  % plot_element_numbers(mesh)
  hold on
  triplot(oldmesh.elem_vertices, ...
  oldmesh.vertex_coordinates(:,1), oldmesh.vertex_coordinates(:,2), 'LineWidth',2)
  hold off
  oldmesh = mesh;
end