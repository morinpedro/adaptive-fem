function plot_element_numbers(mesh, varargin)
%PLOT_ELEMENT_NUMBERS  Label each triangle with its element number.
%
%   plot_element_numbers(mesh) overlays the element index at the
%   centroid of every triangle in MESH, where
%     mesh.elem_vertices        n_elem x 3 vertex indices
%     mesh.vertex_coordinates   n_vertices x 2 coordinates
%   i.e. the same convention used by
%     triplot(mesh.elem_vertices, mesh.vertex_coordinates(:,1), ...
%                                   mesh.vertex_coordinates(:,2))
%
%   Extra name/value pairs are passed through to TEXT and override the
%   defaults below, e.g.
%     plot_element_numbers(mesh, 'FontSize', 20, 'Color', 'r')

  v1 = mesh.vertex_coordinates(mesh.elem_vertices(:,1), :);
  v2 = mesh.vertex_coordinates(mesh.elem_vertices(:,2), :);
  v3 = mesh.vertex_coordinates(mesh.elem_vertices(:,3), :);
  centroids = (v1 + v2 + v3) / 3;

  n_elem = size(mesh.elem_vertices, 1);
  labels = arrayfun(@num2str, (1:n_elem)', 'UniformOutput', false);

  hold on
  text(centroids(:,1), centroids(:,2), labels, ...
       'HorizontalAlignment', 'center', ...
       'VerticalAlignment', 'middle', ...
       'FontSize', 14, varargin{:});
  hold off
end
