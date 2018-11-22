load('real_full_result_2018_11_1_23_4_37.mat')
load('real_full_result_boot.mat')

n_boot = size(impute_log_int, 2);

boot_std = std(impute_log_int, 0, 2);

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
scatter(cx, cy, 12, log(boot_std), 'filled')
colorbar
colormap(jet)
axis image
min_white_margin(gca);
saveas(fig, 'log_boot_std', 'png')

histogram(n_region, 10)
xlabel('Number of segments')
ylabel('Frequency')
min_white_margin(gca);
