fig = figure;
% 10 replicates
cmap = hsv(10);
num_nonempty_all = [];
for ii = 1:10
    load(['real_full_result_', num2str(ii), '.mat'])
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
    % record number of nonempty
    num_nonempty_all = [num_nonempty_all, num_nonempty];

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

    for j = 1:num_nonempty
        line(vx_edges_all{j}, vy_edges_all{j}, 'Color', cmap(ii, :))
    end
end
axis image
