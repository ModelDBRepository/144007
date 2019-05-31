function setPresentationDefaults(figHandle, blackbg, varargin)
% This is a utility function changes the defaults for the current axis


allchildren = allchild(figHandle);
ax = findobj(allchildren, 'Type', 'axes');
set(ax, 'Fontsize', 16, 'LineWidth', 1, varargin{:}, 'TickDir', 'out');
tx = findobj(allchildren, 'Type', 'text');
set(tx, 'Fontsize', 16);
if (length(ax) > 1)
    for i = 1:length(ax)
        th = findobj(ax(i), 'Type', 'text');
        set(th, 'Fontsize', 16);
    end
end

if blackbg
    turnPlotBlack(figHandle);
end


function turnPlotBlack(figHandle)
% Turns all of the axes on a plot black with white labeling

allchildren = allchild(figHandle);
txt = findobj(allchildren, 'Type', 'text');
set(txt, 'Color', 'w');
ax = findobj(allchildren, 'Type', 'axes');
set(ax, 'Color', 'k', 'Xcolor', 'w', 'Ycolor', 'w', 'zcolor', 'w');
lh = findobj(ax, 'Type', 'line', '-and', 'Color', 'k');
th = findobj(ax, 'Type', 'text', '-and', 'Color', 'k');
set(lh, 'Color', [1 1 1]);
set(th, 'Color', 'w');
if (length(ax) > 1)
    for i = 1:length(ax)
        sa = allchild(ax(i));
        tc = findobj(sa, 'Type', 'text', '-and', 'Color', 'k');
        set(tc, 'color', 'w');
    end
end
