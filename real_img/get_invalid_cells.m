function [invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n)
% Gets invalid and valid cells.
% Invalid cells are those with NaN log intensity and isolated cells.

invalid = find(isnan(cell_log_intensity));
valid = setdiff(1:n, invalid)';

constr_graph = graph(adj_mat);
components = conncomp(constr_graph, 'OutputForm', 'cell');
num_comps = length(components);
comp_size = zeros(num_comps, 1);
for i = 1:num_comps
    comp_size(i) = length(components{i});
end
max_comp_size = max(comp_size);
for i = 1:num_comps
    % the connected components that have size smaller than the largest size
    % are supposed to be isolated
    if comp_size(i) < max_comp_size
        isolated_cells = components{i};
        invalid = [invalid; isolated_cells];
        valid = setdiff(valid, isolated_cells);
    end
end

end

