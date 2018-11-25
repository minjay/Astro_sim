function [] = min_white_margin(gca, shift, reduction)

if nargin == 1
    shift = 0;
    reduction = 0;
elseif nargin == 2
    reduction = 0;
end

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left-shift bottom ax_width-reduction ax_height];

end