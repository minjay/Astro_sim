function [cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X)
% conduct some initial computations
% Input variables:
% 
% X: an n-by-2 matrix
%
% Output variables:
%
% cx: the first column of X
% cy: the second column of X
% n: the number of points
% DT: a Delaunay triangulation object
% E: a matrix of edges
% cell_log_intensity: the log intensity of voronoi cells
% cell_area: the area of voronoi cells

cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation(cx, cy);

[V, R] = voronoiDiagram(DT);
E = edges(DT);

n = length(cx);
intensity = zeros(n, 1);
cell_area = zeros(n, 1);
for i = 1:n
    cell_area(i) = polyarea(V(R{i}, 1), V(R{i}, 2));
    intensity(i) = 1/cell_area(i);
end

cell_log_intensity = log(intensity);

end