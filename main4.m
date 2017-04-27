% numerical experiment 4
% several point sources on an extended source

clear
close all

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 5)];
num = [500 50*ones(1, 5)];
% generate simulated data (inhomogeneous Poisson point process)
X = sim_inhomo_Pois_const([0 1], [0 1], 1000, loc, radius, num, 0);

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.075], [0.075 0.05], [0.05 0.02]);
subplot(2, 3, 1)
scatter(X(:, 1), X(:, 2), '.')
axis square

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1]);

% plot log intensity
subplot(2, 3, 2)
triplot(DT)
hold on
scatter(cx, cy, [], cell_log_intensity, 'o', 'filled')
%colorbar
colormap(jet)
axis image

% get seeds
[seeds, num, invalid] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
    0.1, 0.1, 3, cell_log_intensity, cell_area, cx, cy, 2);
disp(['Number of regions is ', num2str(num)])

% plot the seeds
subplot(2, 3, 3)
triplot(DT)
hold on
% specify the colormap
colors = lines(num);
for i = 1:num
    scatter(cx(seeds{i}), cy(seeds{i}), [], colors(i, :), '*')
end
axis image

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
adj_mat = get_adj_mat( E, n );
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

% plot the over-segmented image
subplot(2, 3, 4)
triplot(DT)
hold on
for i = 1:num
    scatter(cx(region_sets{i}), cy(region_sets{i}), [],  colors(i, :), 'filled')
end
axis image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
subplot(2, 3, 5)
plot(num-1:-1:0, BIC_all, '-o', 'MarkerSize', 3)
xlabel('Number of clusters/sources')
ylabel('BIC')
[min_BIC, index_BIC] = min(BIC_all);
text(num-index_BIC-1, min_BIC-150, num2str(num-index_BIC), 'FontSize', 12)
axis square

subplot(2, 3, 6)
triplot(DT)
hold on
% the final result
selected = sets_all{index_BIC};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}), [],  colors(index, :), 'filled')
    end
end
axis image

set(h, 'Position', [0, 0, 800, 500]);
