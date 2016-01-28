function X = sim_homo_Pois(range_x, range_y, lambda)

len_x = range_x(2)-range_x(1);
len_y = range_y(2)-range_y(1);
n = poissrnd(len_x*len_y*lambda);
X = zeros(n, 2);
X(:, 1) = rand(n, 1)*len_x+range_x(1);
X(:, 2) = rand(n, 1)*len_y+range_y(1);

% scatter(X(:, 1), X(:,2), '.')
% axis square

end