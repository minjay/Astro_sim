% numerical experiment 5
% several point sources on an extended source
% L-shaped extended source

clear
close all

GRAY = [0.6 0.6 0.6];

range_x = [0.2 0.4; 0.4 0.6; 0.6 0.8];
range_y = [0.2 0.4; 0.2 0.8; 0.6 0.8];
base_num_in_L_shape = [2 6 2];
base_num_in_circle = ones(1, 3);
factor = 30;
sample_factor = 1;
lambda = 1000;
seed = 1;
% generate simulated data (inhomogeneous Poisson point process)
X = sim_inhomo_Pois_const_L_shape(range_x, range_y, factor * base_num_in_L_shape, seed);

loc = [0.35 0.3; 0.5 0.5; 0.65 0.7];
radius = [0.025 0.025 0.025];
X = [X; sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle)];

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

subplot(1, 3, 1)
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', GRAY, 'LineWidth', 1);
end
hold on
% horizontal lines
for i = 1:3
    plot(range_x(i, :), [range_y(i, 1) range_y(i, 1)], 'Color', GRAY, 'LineWidth', 1)
    plot(range_x(i, :), [range_y(i, 2) range_y(i, 2)], 'Color', GRAY, 'LineWidth', 1)
end
% vertical lines
plot([range_x(1, 1) range_x(1, 1)], range_y(1, :), 'Color', GRAY, 'LineWidth', 1)
plot([range_x(3, 2) range_x(3, 2)], range_y(3, :), 'Color', GRAY, 'LineWidth', 1)
plot([range_x(1, 2) range_x(1, 2)], [range_y(1, 2) range_y(3, 2)], 'Color', GRAY, 'LineWidth', 1)
plot([range_x(2, 2) range_x(2, 2)], [range_y(1, 1) range_y(3, 1)], 'Color', GRAY, 'LineWidth', 1)
scatter(X(:, 1), X(:, 2), 'k.')
axis([0 1 0 1])
axis square
box on

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
subplot(1, 3, 2)
triplot(DT, 'Color', GRAY)
hold on
% specify the colormap
colors = lines(num);
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
for i = 1:num_s_pt
    scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, colors(i + num_s, :), 'd', 'filled')
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 12, 'k', '^', 'filled')
end
axis image

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(1, 3, 3)
triplot(DT, 'Color', GRAY)
hold on
% the final result
selected = sets_all{index_BIC};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}), 12,  colors(index, :), 'filled')
    end
end
axis image
hold on
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', GRAY, 'LineWidth', 1, 'EnhanceVisibility', false);
end
% horizontal lines
for i = 1:3
    plot(range_x(i, :), [range_y(i, 1) range_y(i, 1)], 'Color', GRAY, 'LineWidth', 1)
    plot(range_x(i, :), [range_y(i, 2) range_y(i, 2)], 'Color', GRAY, 'LineWidth', 1)
end
% vertical lines
plot([range_x(1, 1) range_x(1, 1)], range_y(1, :), 'Color', GRAY, 'LineWidth', 1)
plot([range_x(3, 2) range_x(3, 2)], range_y(3, :), 'Color', GRAY, 'LineWidth', 1)
plot([range_x(1, 2) range_x(1, 2)], [range_y(1, 2) range_y(3, 2)], 'Color', GRAY, 'LineWidth', 1)
plot([range_x(2, 2) range_x(2, 2)], [range_y(1, 1) range_y(3, 1)], 'Color', GRAY, 'LineWidth', 1)

set(h, 'Position', [0, 0, 800, 250]);
