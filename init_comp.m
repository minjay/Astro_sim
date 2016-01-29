function [cx, cy, n, DT, E, flux, area_cell] = init_comp(X)

cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation(cx, cy);

[V, R] = voronoiDiagram(DT);
E = edges(DT);

n = length(cx);
flux = zeros(n, 1);
area_cell = zeros(n, 1);
for i = 1:n
    area_cell(i) = polyarea(V(R{i}, 1), V(R{i}, 2));
    flux(i) = 1/area_cell(i);
end

flux = log(flux);

end