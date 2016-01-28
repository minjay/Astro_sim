function X = sim_inhomo_Pois_const(range_x, range_y, lambda, loc, radius, num)
% ex:
% X = sim_inhomo_Pois_const([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.1 0.1], [100 100]);

X = sim_homo_Pois(range_x, range_y, lambda);

n_s = length(radius);
max_iter = 1e3;
for i = 1:n_s
    rand_num = rand(max_iter, 2);
    rand_num(:, 1) = rand_num(:, 1)*radius(i)*2+loc(i, 1)-radius(i);
    rand_num(:, 2) = rand_num(:, 2)*radius(i)*2+loc(i, 2)-radius(i);
    index = find(((rand_num(:, 1)-loc(i, 1)).^2+(rand_num(:, 2)-loc(i, 2)).^2)<=radius(i)^2);
    index = index(1:num(i));
    rand_num = rand_num(index, :);
    X = [X; rand_num];
end

scatter(X(:, 1), X(:,2), '.')
axis square

end