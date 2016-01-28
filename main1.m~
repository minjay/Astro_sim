clear
close all

% X = sim_inhomo_Pois_const([0 1], [0 1], 200, [0.3 0.3; 0.7 0.7], [0.1 0.1], [50 50]);

% generate simulated data (inhomogeneous Poisson point process)
X = sim_inhomo_Pois_Gauss([0 1], [0 1], 200, [0.3 0.3; 0.7 0.7], [0.05 0.05], [50 50]);

cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation([cx cy]);
figure
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

% apply logrithm to fluxes so that the bump is somewhat flattened
flux = log(flux);

hold on
scatter(cx, cy, [], flux, 'o', 'filled')
colorbar

num = 0;
index_set = [];
% union_index_set records the points that have been chosen
union_index_set = [];
% the points (on the boundary) need to be excluded during the process of
% G-SRG
invalid = find(isnan(flux))';
for i = 0.1:0.2:0.9
    for j = 0.1:0.2:0.9
        dist = (cx-i).^2+(cy-j).^2;
        % sort them ascendingly
        [~, index] = sort(dist);
        num = num+1;
        index_set{num} = [];
        % recording the number of points needed for each init region set
        times = 0;
        for k = 1:n
            % choose the points that are not on the boundary and have not
            % been chosen yet
            if ~ismember(index(k), union(invalid, union_index_set))
                times = times+1;
                index_set{num} = [index_set{num} index(k)];
                union_index_set = [union_index_set index(k)];
            end
            % the size of each init region set
            % it could happen that there are not enough points to be
            % assigned when the background is very sparse
            if times==5
                break
            end
        end
    end
end    

% plot of the seeds
figure
triplot(DT)
hold on
% specify the colormap
colors = lines(num);
for i = 1:num
    scatter(cx(index_set{i}), cy(index_set{i}), [], colors(i, :), '*')
end
axis image

% keep the variable index_set and copy it to the variable init_sets
init_sets = cell(num, 1);
for i = 1:num
    init_sets{i} = index_set{i};
end

% graph-based SRG
adj_mat = get_adj_mat( E, n );
[init_sets, labeled_set] = SRG_graph( init_sets, flux, area_cell, n, adj_mat, invalid);

% plot
figure
triplot(DT)
hold on
for i = 1:num
    scatter(cx(init_sets{i}), cy(init_sets{i}), [],  colors(i, :), 'filled')
end
axis image
viscircles([0.3 0.3], 0.1, 'EdgeColor', 'r', 'LineWidth', 1.5);
viscircles([0.7 0.7], 0.1, 'EdgeColor', 'r', 'LineWidth', 1.5);

% note that the flux of region here has not been taken logrithm
flux_region = zeros(num, 1);
area_region = zeros(num, 1);
for i = 1:num
    area_region(i) = sum(area_cell(init_sets{i}));
    flux_region(i) = sum(area_cell(init_sets{i}).*exp(flux(init_sets{i})))/area_region(i);
end

% source detection
% to remove false positive
figure
h_bp = boxplot(flux_region);
h_ol = findobj(h_bp, 'Tag', 'Outliers');
ols = get(h_ol, 'YData');
n_ol = length(ols);

% plot
figure
triplot(DT)
hold on
for i = 1:n_ol
    index_region = find(flux_region==(ols(i)));
    scatter(cx(init_sets{index_region}), cy(init_sets{index_region}), [],  colors(index_region, :), 'filled')
end
axis image
viscircles([0.3 0.3], 0.1, 'EdgeColor', 'r', 'LineWidth', 1.5);
viscircles([0.7 0.7], 0.1, 'EdgeColor', 'r', 'LineWidth', 1.5);


