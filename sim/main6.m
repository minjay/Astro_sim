% numerical experiment 6
% several point sources on an extended source
% ring-shaped extended source

clear
close all

GRAY = [0.6 0.6 0.6];

loc_ring = [0.3 0.5];
radius_out = 0.4;
radius_in = 0.2;
% generate simulated data (inhomogeneous Poisson point process)
X = sim_inhomo_Pois_const_ring([0 1], [0 1], 500, loc_ring, radius_out, radius_in, 250, 0);

loc = [0.5 0.7; 0.6 0.5; 0.5 0.3];
radius = 0.025*ones(1, 3);
num = 25*ones(1, 3);
X = [X; sim_inhomo_Pois_const([0 1], [0 1], 500, loc, radius, num)];

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.075], [0.05 0.05], [0.05 0.075]);

subplot(2, 3, 1)
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', 'k', 'LineWidth', 1);
end
axis([0 1 0 1])
axis square
box on
title('(a) True')
hold on
ang(loc_ring, radius_out, [-pi/2 pi/2], 'k', 1);
ang(loc_ring, radius_in, [-pi/2 pi/2], 'k', 1);
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)-radius_out loc_ring(2)-radius_in], 'k', 'LineWidth', 1)
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)+radius_in loc_ring(2)+radius_out], 'k', 'LineWidth', 1)

subplot(2, 3, 2)
scatter(X(:, 1), X(:, 2), 'k.')
axis square
box on
title('(b) Simulated Data')

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% plot log intensity
subplot(2, 3, 3)
triplot(DT, 'Color', GRAY)
hold on
scatter(cx, cy, 12, cell_log_intensity, 'filled')
hb = colorbar;
set(hb, 'position', [0.94 0.54 0.02 0.41])
colormap(jet)
axis image
title('(c) Constructed Graph')

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
subplot(2, 3, 4)
triplot(DT, 'Color', GRAY)
hold on
% specify the colormap
colors = lines(num_s);
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
for i = 1:num_s_pt
    scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, 'r', 'd', 'filled')
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 12, 'k', '^', 'filled')
end
axis image
title('(d) Seeds')

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

% plot the over-segmented image
subplot(2, 3, 5)
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :), 'filled')
end
for i = 1:num_s_pt
    scatter(cx(region_sets{i+num_s}), cy(region_sets{i+num_s}), 12, 'r', 'filled')
end
axis image
title('(e) Oversegmented Graph')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(2, 3, 6)
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
    viscircles(loc(i, :), radius(i), 'EdgeColor', 'k', 'LineWidth', 0.75, 'EnhanceVisibility', false);
end
ang(loc_ring, radius_out, [-pi/2 pi/2], 'k', 1);
ang(loc_ring, radius_in, [-pi/2 pi/2], 'k', 1);
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)-radius_out loc_ring(2)-radius_in], 'k', 'LineWidth', 1)
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)+radius_in loc_ring(2)+radius_out], 'k', 'LineWidth', 1)
title('(f) Segmentation after Region Merging')

set(h, 'Position', [0, 0, 875, 500]);
