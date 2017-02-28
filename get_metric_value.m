function value = get_metric_value(metric, N, region_sets, region_area, region_intensity, region_num_cells)
% compute the value of metric
%
% Input variables:
%
% metric: metric of interest
% N: total number of points
% region_sets: region sets
% region_area: area of regions
% region_intensity: intensity of regions
% region_num_cells: number of cells in each region
%
% Output variables:
%
% value: value of the metric

% number of regions
K = length(region_sets);

log_like = -sum(log(1:N));
for i = 1:K
    log_like = log_like+region_num_cells(i)*log(region_intensity(i))-...
        region_intensity(i)*region_area(i);
end

if strcmp(metric, 'BIC')
    value = -2*log_like+K*log(N);
elseif strcmp(metric, 'log-like')
    value = log_like;
end
    
end