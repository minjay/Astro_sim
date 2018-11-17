function [sx, sy] = uniform_sample_triangle(vx, vy, n)
% Reference: https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle

r1 = rand(n, 1);
r2 = rand(n, 1);

coef1 = 1 - sqrt(r1);
coef2 = sqrt(r1) .* (1 - r2);
coef3 = r2 .* sqrt(r1);

sx = coef1 * vx(1) + coef2 * vx(2) + coef3 * vx(3);
sy = coef1 * vy(1) + coef2 * vy(2) + coef3 * vy(3);

end