function trace_filt = median_filt(trace, filterTime, dt, varargin)
% function median_filt(trace, filterTime, dt) 
% written by Peter W Jones
%
% Gives a slightly more
% appropriate version of the function which lets the caller specify a 
% filter time width, and a sampling rate for the trace.  The function
% then filters with the correct length filter
% Varargin is a dimension number, to be passed to medfilt1 in order to filter 
% a whole matrix of traces at once, but only if you pass a matrix.

width = round(filterTime / dt);
if (rem(width, 2) == 0) width = width+1; end

if (~isempty(varargin) && min(size(trace)) > 1)
    dim = varargin{1};
    blksiz = size(trace, dim);
    trace_filt = medfilt1(trace, width, [], dim);
else
    trace_filt = medfilt1(trace, width);
end
