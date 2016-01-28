X = sim_inhomo_Pois_const([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.1 0.1], [20 20]);
X = sim_inhomo_Pois_Gauss([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.05 0.05], [20 20]);
cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation([cx cy]);
triplot(DT)

[V, R] = voronoiDiagram(DT);
E = edges(DT);

axis image

n = length(cx);
flux = zeros(n, 1);
area_cell = zeros(n, 1);
for i = 1:n
    area_cell(i) = polyarea(V(R{i}, 1), V(R{i}, 2));
    flux(i) = 1/area_cell(i);
end

flux = log(flux);

hold on
scatter(cx, cy, [], flux, 'o', 'filled')
colorbar

s1 = [0.3 0.3];
s2 = [0.7 0.7];
init_sets = cell(3, 1);
dist1 = (cx-s1(1)).^2+(cy-s1(2)).^2;
dist2 = (cx-s2(1)).^2+(cy-s2(2)).^2;
[~, index1] = min(dist1);
[~, index2] = min(dist2);
% not on the border!
index3 = 1;
init_sets{1} = index1;
init_sets{2} = index2;
init_sets{3} = index3;

figure
triplot(DT)
hold on
scatter([cx(index1) cx(index2) cx(index3)], [cy(index1) cy(index2) cy(index3)], 'r*')
axis image

% graph-based SRG
invalid = find(isnan(flux))';
adj_mat = get_adj_mat( E, n );
[init_sets, labeled_set] = SRG_graph( init_sets, flux, area_cell, n, adj_mat, invalid);

% plot
figure
triplot(DT)
hold on
scatter(cx(init_sets{1}), cy(init_sets{1}), 'ro', 'filled')
scatter(cx(init_sets{2}), cy(init_sets{2}), 'go', 'filled')
scatter(cx(init_sets{3}), cy(init_sets{3}), 'co', 'filled')
axis image
viscircles([0.3 0.3], 0.1, 'EdgeColor', 'y', 'LineWidth', 1.5)
viscircles([0.7 0.7], 0.1, 'EdgeColor', 'y', 'LineWidth', 1.5)
