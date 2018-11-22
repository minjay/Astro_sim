function [invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n)
% Gets invalid and valid cells.
% Invalid cells are those with NaN log intensity and isolated cells.

invalid = find(isnan(cell_log_intensity));
valid = setdiff(1:n, invalid)';

%constr_graph = graph(adj_mat);
%components = conncomp(constr_graph, 'OutputForm', 'cell');
%num_comps = length(components);
%comp_size = zeros(num_comps, 1);
%for i = 1:num_comps
%    comp_size(i) = length(components{i});
%end
% graph function is not available before R2015b.
% instead, following 
% https://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/
adj_mat = adj_mat + eye(size(adj_mat, 1));
[p, ~, r, ~] = dmperm(adj_mat);
num_comps = length(r) - 1;
comp_size = zeros(num_comps, 1);
components = cell(num_comps, 1);
for i = 1:num_comps
    comp_size(i) = r(i + 1) - r(i);
    components{i} = p(r(i):r(i + 1) - 1);
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

