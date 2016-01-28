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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note that the flux of region here has not been taken logrithm
flux_region = zeros(num, 1);
area_region = zeros(num, 1);
for i = 1:num
    area_region(i) = sum(area_cell(init_sets{i}));
    flux_region(i) = sum(area_cell(init_sets{i}).*exp(flux(init_sets{i})))/area_region(i);
end

% get the connectivity among regions
adj_mat_region = zeros(num);
for i = 1:num-1
    for j = i+1:num
        adj_mat_region(i, j) = check_connect(i, j, init_sets, adj_mat);
        adj_mat_region(j, i) = adj_mat_region(i, j);
    end
end

% construct a region merging tree structure
% t represents the index of region merging
region_set = 1:num;
rec_min_diff = zeros(num-1, 1);
rec_sets = cell(num-1, 1);
rec_sets{1} = init_sets;
for t = 1:num-1
    % find the pair of regions with minimal flux difference
    min_diff = realmax;
    for i = region_set
        for j = region_set
            if i<j && adj_mat_region(i, j) && abs(flux_region(i)-flux_region(j))<min_diff
                min_diff = abs(flux_region(i)-flux_region(j));
                merge1 = i;
                merge2 = j;
            end
        end
    end
    % del merge2 from region_set since it will be merged with set merge1
    region_set(region_set==merge2) = [];
    % the i-th element of rec_min_diff represents the min diff when there
    % are i+1 regions
    rec_min_diff(num-t) = min_diff;
    % merge two regions
    flux_region(merge1) = (area_region(merge1)*flux_region(merge1)+area_region(merge2)*flux_region(merge2))/...
        (area_region(merge1)+area_region(merge2));
    area_region(merge1) = area_region(merge1)+area_region(merge2);
    % update the connectivity
    for i = 1:num
        if adj_mat_region(merge2, i)
            adj_mat_region(merge1, i) = 1;
            adj_mat_region(i, merge1) = 1;
        end
    end
    % record the merged sets
    tmp = rec_sets{t};
    tmp{merge1} = [tmp{merge1} tmp{merge2}];
    tmp{merge2} = [];
    % the t-th element of rec_sets represents the clustering when there are
    % num-t+1 regions
    rec_sets{t+1} = tmp;
end

% decide threshold by L-method
figure
plot(2:num, rec_min_diff, '.')
min_rmse = realmax;
for c = 3:num-2
    [b1, ~, r1] = regress(rec_min_diff(1:c-1), [ones(c-1, 1) (2:c)']);
    [b2, ~, r2] = regress(rec_min_diff(c:num-1), [ones(num-c, 1), (c+1:num)']);
    rmse1 = sqrt(mean(r1.^2));
    rmse2 = sqrt(mean(r2.^2));
    rmse = (c-1)/(num-1)*rmse1+(num-c)/(num-1)*rmse2;
    if rmse<min_rmse
        min_rmse = rmse;
        % stop merging when there are knee regions
        knee = c;
        b1_min = b1;
        b2_min = b2;
    end
end
hold on
plot(2:knee, b1_min(1)+b1_min(2)*(2:knee), 'r')
plot(knee+1:num, b2_min(1)+b2_min(2)*(knee+1:num), 'g')

figure
triplot(DT)
hold on
% the final result
selected = rec_sets{num-knee+1};
for i = 1:num
    if ~isempty(selected{i})
        scatter(cx(selected{i}), cy(selected{i}), [],  colors(i, :), 'filled')
    end
end
axis image