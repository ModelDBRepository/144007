function [peakwidth, peak_edge_times] = computePeakWidth(time, trace, peaktime, peakval, direction, baseline)
% function [peakwidth, peak_edge_times] = computePeakWidth(time, trace, peaktime, peakval, direction, baseline)
%
% Function finds the peak width.  The default width measurement is full width at half height (FWHH).
% All of the arguments are self-explanatory, except for maybe direction. 1 is for detection of maxima,
% and -1 is for minima.
pb=0;

peaki = find(time >= peaktime, 1, 'first');
before = trace(peaki:-1:1);
after = trace(peaki:end);
[begin_height, begin_i] = findPercentThresholdCrossing(.5, before, direction, baseline);
begin_i = peaki-(begin_i-1);
[end_height, end_i] = findPercentThresholdCrossing(.5, after, direction, baseline);
end_i = peaki + end_i - 1;
peak_edge_i = [begin_i(2) end_i(2)];
if (sum(isnan(peak_edge_i)) > 0) 
    peakwidth = NaN;
    peak_edge_times = [NaN, NaN];
else
    peak_edge_times = time(peak_edge_i);
    peakwidth = peak_edge_times(2) - peak_edge_times(1);
    heights = [begin_height(2) end_height(2)];
    if pb
        figure; hold on;
        plot(time, trace, 'k');
        plot(peaktime, peakval, 'rx');
        plot(peak_edge_times, heights, 'bx-');  
    end
    if sum(peak_edge_i == 1 | peak_edge_i == length(trace)) %if either of the indices are the edges of the vector
        peakwidth = NaN;
    end
end


% --------------------------------------------------------------------------
function [val, ind] = findPercentThresholdCrossing(percent, data, varargin)
%
% function [val, ind] = findThresholdCrossing(percent)
%
% Finds the value and index of the point which a vector goes from an
% initial mean to a maximum value.  Vector must be greater than 2000
% elements.  Takes the baseline from the first elements, so that should be
% the baseline period.  So, it finds the threshold crossing, and then finds
% the return below the threshold. 
% Varargin{1} is the direction of the thresholding.  1 is positive relative
% to baseline, -1 is negative, and 0 is autodetect based on relative
% amplitudes of signal from baseline.
% Varargin{2} is the baseline level.  Normally this is automatically calculated 
% from the beginning of the trace. This gives a way to set it manually.
pb = 0;

dir = 0;
if(~isempty(varargin))
    n_vararg = nargin - 2;
    dir = varargin{1};   
else
    n_vararg = 0;
end
ymin = nanmin2(data);
ymax = nanmax2(data);
if n_vararg >=2 %Arg 2 is the optional baseline argument
    y0 = varargin{2};
else %things to do without a baseline value
    try %compute it from the beginning of the trace
        y0 = nanmean2(data(1:500)); 
    catch % or if that doesn't work, set it at zero
        disp('Defaulting baseline value to zero');
        y0 = 0;
    end
end
if ((dir == 0) && (abs(ymax-y0) > abs(y0 - ymin)))
    dir = 1;
    threshold = (ymax-y0)*percent + y0;
    maxi = find(data == ymax);
elseif (dir == 0) 
    dir = -1;
    threshold = (ymin-y0)*percent + y0;
    maxi = find(data == ymin);
elseif (dir == 1)
    threshold = (ymax-y0)*percent + y0;
    maxi = find(data == ymax);
else
    threshold = -1*(y0-ymin)*percent + y0;
    maxi = find(data == ymin);
end

if dir==1
    tempval = find(data >= threshold);
else
    tempval = find(data <= threshold);
end

if(isempty(tempval)) %in case there is nothing that meets threshold
    ind = [NaN NaN];
    val = [NaN NaN];
    return;
end
try
    ind(1) = tempval(1);
    %now find the breaks in the area above threshold
    spacing = diff(tempval);
    breaks = find(spacing > 1 & tempval(2:end) > maxi(1)); %find the dips under thresh after the peak
    if ~isempty(breaks)
        ind(2) = tempval(breaks(1));
    else
        ind(2) = tempval(end);
    end
    val = data(ind);
    ind = ind';
catch
    ind(2) = NaN;
    val(2) = NaN;
end

if (pb)
    figure; hold on;
    plot(data, 'k');
    plot(ind, val, 'rx');
end
