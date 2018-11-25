function [] = add_annotation(loc, subplot_pos, message)
% add custom annotation

abs_loc = [...
    subplot_pos(1) + loc(1) * subplot_pos(3),...
    subplot_pos(2) + loc(2) * subplot_pos(4),...
    loc(3),...
    loc(4)...
    ];
annotation('textbox', abs_loc, 'String', message,...
    'FitBoxToText', 'on', 'FontSize', 12, 'HorizontalAlignment', 'left',...
    'VerticalAlignment', 'middle', 'EdgeColor', 'none');

end