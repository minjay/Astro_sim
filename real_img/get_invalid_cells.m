function [invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n)
% Gets invalid and valid cells.
% Invalid cells are those with NaN log intensity and isolated cells.

invalid = find(isnan(cell_log_intensity));
valid = setdiff(1:n, invalid)';

adj_mat_valid = adj_mat(valid, valid);
num_nbrs = sum(adj_mat_valid, 2);
isolated_cells = valid(num_nbrs == 0);

invalid = [invalid; isolated_cells];
valid = setdiff(valid, isolated_cells);

end

