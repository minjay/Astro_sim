function plot_arc(loc_ring, radius_in, radius_out)

GRAY = [0.6 0.6 0.6];

ang(loc_ring, radius_out, [-pi/2 pi/2], GRAY, 1);
ang(loc_ring, radius_in, [-pi/2 pi/2], GRAY, 1);
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)-radius_out loc_ring(2)-radius_in], 'Color', GRAY, 'LineWidth', 1)
plot([loc_ring(1) loc_ring(1)], [loc_ring(2)+radius_in loc_ring(2)+radius_out], 'Color', GRAY, 'LineWidth', 1)

end