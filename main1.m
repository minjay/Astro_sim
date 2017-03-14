clear
close all

% generate simulated data (inhomogeneous Poisson point process)
% there is only one extended source
X = sim_inhomo_Pois_Gauss([0 1], [0 1], 200, [0.5 0.5], 0.1, 200, 1);

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1]);

% plot log intensity
figure
triplot(DT)
hold on
scatter(cx, cy, [], cell_log_intensity, 'o', 'filled')
colorbar
colormap(jet)
axis image

% get seeds
[seeds, num, invalid] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 3, cell_log_intensity, cx, cy);

% plot the seeds
figure
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
figure
triplot(DT)
hold on
for i = 1:num
    scatter(cx(region_sets{i}), cy(region_sets{i}), [],  colors(i, :), 'filled')
end
viscircles([0.5 0.5], 0.2, 'EdgeColor', 'r', 'LineWidth', 1.5);
axis image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the intensity of regions and the connectivity among regions
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

% compute the log-likelihood of the over-segmented graph before region merging
log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);

% construct a region merging tree
% t represents the iter # of region merging
region_indices = 1:num;
log_like_all = zeros(num, 1);
log_like_all(1) = log_like;
sets_all = cell(num, 1);
sets_all{1} = region_sets;
for t = 1:num-1
    % find the pair of regions with the maximal log-likelihood increasement
    max_inc = -realmax;
    for i = region_indices
        for j = region_indices
            if i<j && adj_mat_region(i, j)
                numerator = region_intensity(i)*region_area(i)+region_intensity(j)*region_area(j);
                sum_area = region_area(i)+region_area(j);
                delta_log_like = region_num_cells(i)*log(numerator/region_intensity(i)/sum_area)+...
                    region_num_cells(j)*log(numerator/region_intensity(j)/sum_area);
                % delta_log_like is defined as new log-like minus old
                % log-like
                % the larger, the better
                if delta_log_like>max_inc
                    max_inc = delta_log_like;
                    merge_region1 = i;
                    merge_region2 = j;
                end
            end
        end
    end
    % del merge_region2 from region_indices since it will be merged with merge_region1
    region_indices(region_indices==merge_region2) = [];
    % update log-likelihood
    log_like_all(t+1) = log_like_all(t)+max_inc;
    % merge merge_region1 and merge_region2
    % only keep merge_region1
    region_intensity(merge_region1) = (region_intensity(merge_region1)*region_area(merge_region1)+...
        region_intensity(merge_region2)*region_area(merge_region2))/...
        (region_area(merge_region1)+region_area(merge_region2));
    region_area(merge_region1) = region_area(merge_region1)+region_area(merge_region2);
    region_num_cells(merge_region1) = region_num_cells(merge_region1)+region_num_cells(merge_region2);
    % update the connectivity
    index = setdiff(find(adj_mat_region(merge_region2, :)==1), merge_region1);
    adj_mat_region(merge_region1, index) = 1;
    adj_mat_region(index, merge_region1) = 1;
    % record the sets after region merging
    tmp = sets_all{t};
    tmp{merge_region1} = [tmp{merge_region1} tmp{merge_region2}];
    tmp{merge_region2} = [];
    sets_all{t+1} = tmp;
end

figure
plot(log_like_all, '-o')

BIC_all = -2*log_like_all+(num:-1:1)'*log(n);
figure
plot(BIC_all, '-o')

figure
triplot(DT)
hold on
% the final result
selected = sets_all{24};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}), [],  colors(index, :), 'filled')
    end
end
axis image