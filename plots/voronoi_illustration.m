rng(0)

n = 50;
X = rand([n, 2]);

fig = figure;
cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation(cx, cy);
% default blue color: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#2
triplot(DT, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2)
hold on
scatter(cx, cy, 48, [0.8500, 0.3250, 0.0980], 'filled')
title('Delaunay triangulation')
axis image
min_white_margin(gca);

fig = figure;
voronoi(cx, cy);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
hold on
scatter(cx, cy, 48, [0.8500, 0.3250, 0.0980], 'filled')
title('Voronoi tessellation')
axis image
min_white_margin(gca);
