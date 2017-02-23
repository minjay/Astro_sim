function [seeds, num_s, invalid] = get_seeds_sim(st_x, en_x, st_y, en_y,...
    step_x, step_y, set_size, cell_log_intensity, cx, cy)
% get initial seeds, which are uniformly spreaded in the region
%
% Input variables:
%
% st_x: the starting value of x
% en_x: the end value of x
% st_y: the starting value of y
% en_y: the end value of y
% step_x: the step size of x
% step_y: the step size of y
% set_size: the size of each seed set
% cell_log_intensity: the log intensity of voronoi cells
% cx: the x coordinate of points
% cy: the y coordinate of points
%
% Output variables:
%
% seeds: seed sets
% num_s: the number of seed sets
% invalid: the set of points that are on the boundary

n = length(cx);

% counter of number of seeds
num_s = 0;
% empty cell
seeds = {};

% set of points that have been selected
unselected_points = 1:n;
% set of points that are on the boundary
invalid = find(isnan(cell_log_intensity))';
unselected_points = setdiff(unselected_points, invalid);

% place seeds uniformly
for i = st_x:step_x:en_x
    for j = st_y:step_y:en_y
        % compute the distance between the photons/points and the grid point
        dist = (cx(unselected_points)-i).^2+(cy(unselected_points)-j).^2;
        % sort them ascendingly
        [~, index] = sort(dist);
        num_s = num_s+1;
        seeds{num_s} = unselected_points(index(1:set_size));
        % remove selected points
        unselected_points = setdiff(unselected_points, seeds{num_s});
    end
end    

end