function X = sim_inhomo_Pois_Gauss(range_x, range_y, lambda, loc, sigma, num)
% ex:
% X = sim_inhomo_Pois_Gauss([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.05 0.05], [100 100]);

X = sim_homo_Pois(range_x, range_y, lambda);

n_s = length(sigma);
max_iter = 1e3;
for i = 1:n_s
    rand_num = mvnrnd(loc(i, :)', sigma(i)^2*eye(2), max_iter);
    index = find(((rand_num(:, 1)-loc(i, 1)).^2+(rand_num(:, 2)-loc(i, 2)).^2)<=2*sigma(i));
    index = index(1:num(i));
    rand_num = rand_num(index, :);
    X = [X; rand_num];
end

scatter(X(:, 1), X(:,2), '.')
axis square

end