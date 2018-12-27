% Circular extended source

clear
close all

GRAY = [0.6 0.6 0.6];

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];

base_num_in_circle = [10 ones(1, 4)];

factor = 30;
sample_factor = 1;
lambda = 1000;
seed = 0;

X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle, seed);

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

subplot(1, 3, 1)
plot_circles(loc, radius)
hold on
scatter(X(:, 1), X(:, 2), 'k.')
axis([0 1 0 1])
axis square
box on
set(gca, 'fontsize', 12)

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
% specify the colormap
colors = lines(num);
plot_seeds(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
set(gca, 'fontsize', 12)

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(1, 3, 3)
plot_segmentation(DT, index_BIC, sets_all, cx, cy, colors)
hold on
plot_circles(loc, radius)
axis image
set(gca, 'fontsize', 12)

set(h, 'Position', [0, 0, 800, 260]);
