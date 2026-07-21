function elems = click_elements(mesh)
%CLICK_ELEMENTS  Interactively click on a triangulated mesh and return
%   the element numbers of the triangles clicked.
%
%   elems = click_elements(mesh) assumes the mesh is already plotted in
%   the current figure, e.g. via
%     triplot(mesh.elem_vertices, mesh.vertex_coordinates(:,1), ...
%                                   mesh.vertex_coordinates(:,2))
%   Left-click inside a triangle to select it (the element is highlighted
%   and its number is printed); press Enter or right-click to stop.
%   Returns a column vector of the element indices selected, in the
%   order they were clicked (clicks outside the mesh are ignored).
%
%   mesh.elem_vertices        n_elem x 3 vertex indices
%   mesh.vertex_coordinates   n_vertices x 2 coordinates

  elems = [];
  hold on
  fprintf('Click on elements (Enter or right-click to stop)...\n');
  while true
    [x, y, button] = ginput(1);
    if isempty(x) || button ~= 1
      break
    end
    T = find_element(mesh, x, y);
    if isempty(T)
      fprintf('  (%.4g, %.4g): outside the mesh\n', x, y);
      continue
    end
    fprintf('  (%.4g, %.4g): element %d\n', x, y, T);
    elems(end+1,1) = T; %#ok<AGROW>

    v = mesh.elem_vertices(T,:);
    patch(mesh.vertex_coordinates(v,1), mesh.vertex_coordinates(v,2), ...
          'blue', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    ctr = mean(mesh.vertex_coordinates(v,:), 1);
    text(ctr(1), ctr(2), num2str(T), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', 'FontSize', 14, 'FontWeight', 'bold');
  end
  hold off
end

function T = find_element(mesh, x, y)
%FIND_ELEMENT  Index of the triangle containing (x,y), or [] if none.
%   Vectorized barycentric-coordinate point-in-triangle test over all
%   elements at once.

  v1 = mesh.vertex_coordinates(mesh.elem_vertices(:,1), :);
  v2 = mesh.vertex_coordinates(mesh.elem_vertices(:,2), :);
  v3 = mesh.vertex_coordinates(mesh.elem_vertices(:,3), :);

  d = (v2(:,2)-v3(:,2)).*(v1(:,1)-v3(:,1)) + (v3(:,1)-v2(:,1)).*(v1(:,2)-v3(:,2));
  a = ((v2(:,2)-v3(:,2)).*(x-v3(:,1)) + (v3(:,1)-v2(:,1)).*(y-v3(:,2))) ./ d;
  b = ((v3(:,2)-v1(:,2)).*(x-v3(:,1)) + (v1(:,1)-v3(:,1)).*(y-v3(:,2))) ./ d;
  c = 1 - a - b;

  tol = -1e-10;
  idx = find(a >= tol & b >= tol & c >= tol, 1);
  if isempty(idx)
    T = [];
  else
    T = idx;
  end
end
