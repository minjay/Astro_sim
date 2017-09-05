addpath(genpath('/home/minjay/G-SRG'))
addpath(genpath('/home/minjay/Astro_sim'))

GRAY = [0.6 0.6 0.6];

% read data
filename='photon_loc.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X = A.data;
clear A

% normalize locations to [0, 1] x [0, 1]
min_x = min(X(:, 1));
max_x = max(X(:, 1));
min_y = min(X(:, 2));
max_y = max(X(:, 2));
X(:, 1) = (X(:, 1) - min_x) ./ (max_x - min_x);
X(:, 2) = (X(:, 2) - min_y) ./ (max_y - min_y);

fig = figure;
scatter(X(:, 1), X(:, 2), '.')
axis image
saveas(fig, 'data', 'epsc')

disp('Conducting some initial computations...')
% init comp
bound_x = [0 1];
bound_y = [0 1];
n = size(X, 1);
count = ones(n, 1);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, bound_x, bound_y, count);

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
invalid = find(isnan(cell_log_intensity));
valid = setdiff(1:n, invalid);
scatter(cx(valid), cy(valid), 12, cell_log_intensity(valid), 'filled')
colorbar
colormap(jet)
axis image
saveas(fig, 'log_intensity', 'epsc')

disp('Getting initial seeds...')
% get seeds
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt, invalid] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.1, 0.1, 10, cell_log_intensity, cell_area, cx, cy, 2, 100, 10);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
fig = figure;
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
saveas(fig, 'seeds', 'epsc')

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

disp('Seeded region growing...')
% graph-based SRG
adj_mat = get_adj_mat( E, n );
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

% plot the over-segmented image
fig = figure;
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :), 'filled')
end
for i = 1:num_s_pt
    scatter(cx(region_sets{i+num_s}), cy(region_sets{i+num_s}), 12, 'r', 'filled')
end
axis image
saveas(fig, 'over_seg', 'epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Region merging...')
[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

fig = figure;
plot(num-1:-1:0, BIC_all, '-o', 'MarkerSize', 3)
saveas(fig, 'BIC', 'epsc')

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% the final result
selected = sets_all{index_BIC};
for i = 1:num
    if ~isempty(selected{i})
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i})));
        scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
saveas(fig, 'over_seg', 'epsc')
