function n_refined = refine_element(elem_to_refine)
% function refine_element(elem_to_refine)
% this function refines the indicated element
% and its corresponding neighbor, if the marked 
% sides coincide. 
% Otherwise, asks for a refinement of the neighbor

global adapt_global_mesh adapt_global_uh adapt_global_fh

n_refined = 0;

neighbour = adapt_global_mesh.elem_neighbours(elem_to_refine, 3);

if (neighbour > 0) % not a boundary side
  if (adapt_global_mesh.elem_neighbours(neighbour, 3) ~= elem_to_refine)  % not compatible
    n_refined = n_refined + refine_element(neighbour);
  end
end

n_elem = adapt_global_mesh.n_elem;
n_vertices = adapt_global_mesh.n_vertices;

% now the neighbour is compatible and comes the refinement
% the neighbour information may have changed
neighbour = adapt_global_mesh.elem_neighbours(elem_to_refine, 3);

v_elem = adapt_global_mesh.elem_vertices(elem_to_refine,:);
neighs_elem = adapt_global_mesh.elem_neighbours(elem_to_refine,:); 
bdries_elem = adapt_global_mesh.elem_boundaries(elem_to_refine,:); 

n_vertices = n_vertices + 1;
adapt_global_mesh.vertex_coordinates(n_vertices, :) = 0.5*(adapt_global_mesh.vertex_coordinates(v_elem(1),:) + adapt_global_mesh.vertex_coordinates(v_elem(2),:));
adapt_global_uh(n_vertices, :) = 0.5*(adapt_global_uh(v_elem(1), :) + adapt_global_uh(v_elem(2), :));
adapt_global_fh(n_vertices, :) = 0.5*(adapt_global_fh(v_elem(1), :) + adapt_global_fh(v_elem(2), :));
 
if (neighbour == 0) % marked side is at the boundary
  n_elem = n_elem + 1;
  adapt_global_mesh.elem_vertices(elem_to_refine, :) = [v_elem(3) v_elem(1) n_vertices];
  adapt_global_mesh.elem_vertices(n_elem, :)         = [v_elem(2) v_elem(3) n_vertices];

  % update neighbour information
  adapt_global_mesh.elem_neighbours(elem_to_refine, :) = [neighs_elem(3) n_elem neighs_elem(2)];
  adapt_global_mesh.elem_neighbours(n_elem, :) = [elem_to_refine neighs_elem(3) neighs_elem(1)];
  % update boundary information
  adapt_global_mesh.elem_boundaries(elem_to_refine, :) = [bdries_elem(3) 0 bdries_elem(2)];
  adapt_global_mesh.elem_boundaries(n_elem, :) = [0 bdries_elem(3) bdries_elem(1)];

  % update neighbour information of the neighbours
  if (neighs_elem(1) > 0)
    f = find(adapt_global_mesh.elem_neighbours(neighs_elem(1),:) == elem_to_refine);
    adapt_global_mesh.elem_neighbours(neighs_elem(1),f) = n_elem;
  end

  adapt_global_mesh.mark(n_elem) = max(adapt_global_mesh.mark(elem_to_refine) - 1, 0);
  adapt_global_mesh.mark(elem_to_refine) = max(adapt_global_mesh.mark(elem_to_refine) - 1, 0);
  n_refined = n_refined + 1;
else
  v_neigh = adapt_global_mesh.elem_vertices(neighbour,:);
  neighs_neigh = adapt_global_mesh.elem_neighbours(neighbour,:); 
  bdries_neigh = adapt_global_mesh.elem_boundaries(neighbour,:); 

  n_elem = n_elem + 2;
  adapt_global_mesh.elem_vertices(elem_to_refine, :) = [v_elem(3)  v_elem(1) n_vertices];
  adapt_global_mesh.elem_vertices(n_elem-1, :)       = [v_elem(2)  v_elem(3) n_vertices];
  adapt_global_mesh.elem_vertices(neighbour, :)      = [v_neigh(3) v_neigh(1) n_vertices];
  adapt_global_mesh.elem_vertices(n_elem, :)         = [v_neigh(2) v_neigh(3) n_vertices];

  % update neighbour information
  adapt_global_mesh.elem_neighbours(elem_to_refine, :) = [n_elem n_elem-1 neighs_elem(2)];
  adapt_global_mesh.elem_neighbours(n_elem-1, :) = [elem_to_refine neighbour neighs_elem(1)];
  adapt_global_mesh.elem_neighbours(neighbour, :) = [n_elem-1 n_elem neighs_neigh(2)];
  adapt_global_mesh.elem_neighbours(n_elem, :) = [neighbour elem_to_refine neighs_neigh(1)];
  % update boundary information
  adapt_global_mesh.elem_boundaries(elem_to_refine, :) = [0 0 bdries_elem(2)];
  adapt_global_mesh.elem_boundaries(n_elem-1, :)       = [0 0 bdries_elem(1)];
  adapt_global_mesh.elem_boundaries(neighbour, :)      = [0 0 bdries_neigh(2)];
  adapt_global_mesh.elem_boundaries(n_elem, :)         = [0 0 bdries_neigh(1)];

  % update neighbour information of the neighbours
  if (neighs_elem(1) > 0)
    f = find(adapt_global_mesh.elem_neighbours(neighs_elem(1),:) == elem_to_refine);
    adapt_global_mesh.elem_neighbours(neighs_elem(1),f) = n_elem-1;
  end
  if (neighs_neigh(1) > 0)
    f = find(adapt_global_mesh.elem_neighbours(neighs_neigh(1),:) == neighbour);
    adapt_global_mesh.elem_neighbours(neighs_neigh(1),f) = n_elem;
  end

  adapt_global_mesh.mark(n_elem-1) = max(adapt_global_mesh.mark(elem_to_refine) - 1, 0);
  adapt_global_mesh.mark(elem_to_refine) = max(adapt_global_mesh.mark(elem_to_refine) - 1, 0);

  adapt_global_mesh.mark(n_elem) = max(adapt_global_mesh.mark(neighbour) - 1, 0);
  adapt_global_mesh.mark(neighbour) = max(adapt_global_mesh.mark(neighbour) - 1, 0);
  n_refined = n_refined + 2;
end

adapt_global_mesh.n_elem = n_elem;
adapt_global_mesh.n_vertices = n_vertices;

