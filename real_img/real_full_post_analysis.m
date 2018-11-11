load('real_full_result_2018_11_1_23_4_37.mat')

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
min_white_margin(gca);
saveas(fig, 'seeds', 'png')

% drop empty entries in the cell array
num = length(selected);

area_all = [];
log_int_all = [];
selected_nonempty = {};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        selected_nonempty{index} = selected{i};
        area = sum(cell_area(selected{i}));
        area_all = [area_all, area];
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
        log_int_all = [log_int_all, log_int];
    end
end

num_nonempty = length(selected_nonempty);

% get source boundaries formed by Voronoi cell edges
[V, R] = voronoiDiagram(DT);
vx_edges_all = {};
vy_edges_all = {};
for j = 1:num_nonempty
    cells_in_region = R(selected_nonempty{j});
    edges = [];
    for i = 1:length(cells_in_region)
        % each row records two vertex indices of an edge
        new_edges = [cells_in_region{i}' [cells_in_region{i}(2:end)'; cells_in_region{i}(1)]];
        % sort vertex indices to avoid ambiguity
        new_edges = sort(new_edges, 2);
        edges = [edges; new_edges];
    end
    [unique_edges, ~ , ind] = unique(edges, 'rows');
    % get the count of each unique edge
    counts = histc(ind, unique(ind));
    % remove edges that appear at least twice to get the boundary edges
    edges = unique_edges(counts==1, :);
    % x-axis
    vx_edges = [V(edges(:, 1), 1)'; V(edges(:, 2), 1)'];
    % y-axis
    vy_edges = [V(edges(:, 1), 2)'; V(edges(:, 2), 2)'];
    nume = size(vx_edges, 2);
    vx_edges = [vx_edges; NaN(1, nume)];
    vx_edges = vx_edges(:);
    vx_edges_all{j} = vx_edges;
    vy_edges = [vy_edges; NaN(1, nume)];
    vy_edges = vy_edges(:);
    vy_edges_all{j} = vy_edges;
end

fig = figure;
plot(cx, cy, '.')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'color', 'r')
end
axis image

left_corner1 = [0.2 0.45];
right_corner1 = [0.45 0.6];
rectangle('Position', [left_corner1 right_corner1 - left_corner1])
text(left_corner1(1), right_corner1(2)+0.01, 'Region 1')

left_corner2 = [0.5 0.35];
right_corner2 = [0.75 0.6];
rectangle('Position', [left_corner2 right_corner2 - left_corner2])
text(left_corner2(1), right_corner2(2)+0.01, 'Region 2')

left_corner3 = [0.3 0.05];
right_corner3 = [0.55 0.2];
rectangle('Position', [left_corner3 right_corner3 - left_corner3])
text(left_corner3(1), right_corner3(2)+0.01, 'Region 3')

min_white_margin(gca);
saveas(fig, 'segment_result', 'png')

fig = figure;
plot(cx, cy, '.')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'color', 'r')
end
axis image
axis([left_corner1(1) right_corner1(1) left_corner1(2) right_corner1(2)])
title('Region 1')
min_white_margin(gca);
saveas(fig, 'region1', 'png')

fig = figure;
plot(cx, cy, '.')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'color', 'r')
end
axis image
axis([left_corner2(1) right_corner2(1) left_corner2(2) right_corner2(2)])
title('Region 2')
min_white_margin(gca);
saveas(fig, 'region2', 'png')

fig = figure;
plot(cx, cy, '.')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'color', 'r')
end
axis image
axis([left_corner3(1) right_corner3(1) left_corner3(2) right_corner3(2)])
title('Region 3')
min_white_margin(gca);
saveas(fig, 'region3', 'png')

% point sources (to be compared with the result from the wavdetect
% algorithm)
% not to be fancy, choose a cutoff point
cutoff_point_source = 0.002;

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_nonempty
    % if area is less than the threshold
    if area_all(i) < cutoff_point_source
        scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int_all(i)*ones(length(selected_nonempty{i}), 1), 'filled')
    end
end
colorbar('EastOutside')
colormap(hsv)
axis image
min_white_margin(gca);
saveas(fig, 'point_sources', 'png')

fig = figure;
plot(num-1:-1:0, BIC_all, '-o', 'MarkerSize', 3)
y_range = get(gca, 'ylim');
hold on
plot([num-index_BIC num-index_BIC],y_range, 'LineWidth', 1.5)
axis tight
position = get(gcf, 'Position');
% shrink width and height
set(gcf, 'Position', [position(1), position(2), position(3)/1.5, position(4)/1.5])
min_white_margin(gca);

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
colorbar('EastOutside')
colormap(hsv)
axis image
min_white_margin(gca);
saveas(fig, 'segment_result_points', 'png')
