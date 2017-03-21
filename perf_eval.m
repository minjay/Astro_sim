function [dr, far, err, res_source_area_sort, source_x, source_y, n_S_correct] =...
    perf_eval(n_source, loc, radius, num_of_photon, cx, cy, selected, cell_log_intensity, cell_area)
% performance evaluation

% compute/estimate the intensity of the fitted regions
region_area = zeros(n_source+1, 1);
region_intensity = zeros(n_source+1, 1);
for i = 1:n_source+1
    region_area(i) = sum(cell_area(selected{i}));
    region_intensity(i) = sum(cell_area(selected{i}).*exp(cell_log_intensity(selected{i})))/region_area(i);
end
% sort regions by their intensities
[~, res_index_sort] = sort(region_intensity, 'descend');
% obtain the fitted source sets, sorted in descending order of
% intensities
res_source_sort = selected(res_index_sort(1:end-1));
res_source_area_sort = region_area(res_index_sort(1:end-1));
% compute/estimate the intensity of the true sources
% only works for round-shaped sources
est_source_int = num_of_photon./(pi*radius.^2);
% sort sources by their intensities
[~, data_index_sort] = sort(est_source_int, 'descend');
% obtain the true location, radius of the sources, sorted in descending 
% order of intensities
loc_sort = loc(data_index_sort, :);
radius_sort = radius(data_index_sort);
% n1 records the number of successfully detected source photons for
% each source
% n2 records the number of false alarms for each source
n1 = zeros(n_source, 1);
n2 = zeros(n_source, 1);
% n_S records the number of photons in S_i
% n_Sc records the number of photons in S_i^c
% also include photons from the background
n_S = zeros(n_source, 1);
n_Sc = zeros(n_source, 1); 
% estimated source position
source_x = zeros(n_source, 1);
source_y = zeros(n_source, 1);
n_S_correct = zeros(n_source, 1);
for i = 1:n_source
    dist_sq = (cx(res_source_sort{i})-loc_sort(i, 1)).^2+(cy(res_source_sort{i})-loc_sort(i, 2)).^2;
    n1(i) = length(find(dist_sq<=radius_sort(i)^2));
    n2(i) = length(find(dist_sq>radius_sort(i)^2));
    dist_sq = (cx-loc_sort(i, 1)).^2+(cy-loc_sort(i, 2)).^2;
    n_S(i) = length(find(dist_sq<=radius_sort(i)^2));
    n_Sc(i) = length(find(dist_sq>radius_sort(i)^2));
    denom = sum(exp(cell_log_intensity(res_source_sort{i})));
    source_x(i) = sum(exp(cell_log_intensity(res_source_sort{i})).*cx(res_source_sort{i}))/denom;
    source_y(i) = sum(exp(cell_log_intensity(res_source_sort{i})).*cy(res_source_sort{i}))/denom;
    % with background correction
    % min(region_intensity) is supposed to be the intensity of the
    % background
    n_S_correct(i) = length(res_source_sort{i})-min(region_intensity)*res_source_area_sort(i);
end
% detection rate
dr = n1./n_S;
% false alarm rate
far = n2./n_Sc;
% total misclassification rate
err = 1-dr+far;
   
end