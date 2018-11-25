g1 = [ones(1500, 1); 2*ones(1500, 1); 3*ones(1500, 1)];
g2 = [ones(500, 3); 2*ones(500, 3); 3*ones(500, 3)];
g2 = g2(:);
shift = 0.01;

load('sim4_result.mat')
figure
boxplot_group(metrics', g1, g2, '(a) Circular extended source')
min_white_margin(gca, shift)

load('sim5_result.mat')
figure
boxplot_group(metrics', g1, g2, '(b) Z-shaped extended source')
min_white_margin(gca, shift)

load('sim6_result.mat')
figure
boxplot_group(metrics', g1, g2, '(c) Arc-shaped extended source')
min_white_margin(gca, shift)


