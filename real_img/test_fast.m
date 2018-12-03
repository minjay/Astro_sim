% test whether fast region merging is the same as region merging
load('real_full_fast_result_2018_12_2_21_20_30.mat')
selected_fast = selected;

load('real_full_result_2018_11_1_23_4_37.mat')

index = 0;
selected_nonempty = {};
selected_size = [];
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        selected_nonempty{index} = selected{i};
        selected_size = [selected_size length(selected{i})];
    end
end

index = 0;
selected_fast_nonempty = {};
selected_fast_size = [];
for i = 1:length(selected_fast)
    if ~isempty(selected_fast{i})
        index = index + 1;
        selected_fast_nonempty{index} = selected_fast{i};
        selected_fast_size = [selected_fast_size length(selected_fast{i})];
    end
end

assert(length(selected_nonempty) == length(selected_fast_nonempty))

assert(all(sort(selected_size) == sort(selected_fast_size)))