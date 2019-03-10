% load fitted result
load('real_full_result_2018_11_1_23_4_37.mat')

% read data (in its original scale)
filename='photon_loc.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X = A.data;
clear A

% init the last column (recording the segment id)
X(:, 3) = 0;

num = length(selected);
area_all = [];
log_int_all = [];
selected_nonempty = {};
index = 0;

for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        selected_nonempty{index} = selected{i};
        area = sum(cell_area(selected{i}));
        area_all = [area_all, area];
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
        log_int_all = [log_int_all, log_int];
    end
end

% sorted by log intensity descendingly
[~, sorted_index] = sort(log_int_all, 'descend');
selected_nonempty = selected_nonempty(sorted_index);

num_segments = length(selected_nonempty);
for segment_id = 1:num_segments
    X(selected_nonempty{segment_id}(1, :), 3) = segment_id;
end

save('segment_id.txt', 'X', '-ascii')

for segment_id = 1:num_segments
    h = figure;
    scatter(X(:, 1), X(:, 2), '.')
    hold on
    x_segment = X(X(:, 3) == segment_id, 1);
    y_segment = X(X(:, 3) == segment_id, 2);
    convhull_points = convhull(x_segment, y_segment);
    scatter(x_segment, y_segment, '.')
    plot(x_segment(convhull_points), y_segment(convhull_points), '-')
    axis tight
    title(['Segment ', num2str(segment_id)])
    saveas(h, ['convhull_segment_', num2str(segment_id), '.png'])
    close all
end

% write to movie
writerObj = VideoWriter('segments.avi');
writerObj.FrameRate = 1;
open(writerObj);
for segment_id = 1:num_segments
  thisimage = imread(['convhull_segment_', num2str(segment_id), '.png']);
  writeVideo(writerObj, thisimage);
end
close(writerObj);
