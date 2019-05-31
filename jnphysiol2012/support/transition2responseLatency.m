function latencies = transition2responseLatency(transitionDurations, cell_type, disp_type)
% function lgmd_resps = transitionDuration2lgmdResponse(transitionDuration, type)
%
% 2010, Peter W Jones
% This function is used to translate a transition duration of luminance change 
% on the monitor (from an actual looming stimulus), to a response latency (time to peak from start of the 
% stimulus for most measures). 
% TYPE is a string specifying either, 'lmc', 'cc' or 'vc' depending on the relationship the caller wants.

monitor_phot_linp = [0.4540 18.0145];
phot.linp =  [0.4740 23.4226]; %The latencies of photoreceptor onsets, for the projector.  
adjustTransitionDurations = @(x)(polyval(monitor_phot_linp, x) - phot.linp(2))./ phot.linp(1); % Transform between latencies of projector and monitor responses due to different luminances.  Idea being to calculate the equivilant projector responses for monitor stimuli in order to use the fits from projector-based data.

% Parameters for the linear fits between projector luminance changes and their latencies
lmc.linp= [0.40  37]; %latencies of LMC response peaks, in ms 
vc.linp = [0.52  66]; %latencies of the medullary inputs from the LGMD VC data, in ms
cc.linp = [0.49  84]; %latencies for the LGMD current clamp responses, in ms

if(strcmp(disp_type, 'monitor'))
    transitionDurations = adjustTransitionDurations(transitionDurations);
end
ts = eval(cell_type);
latencies = polyval(ts.linp, transitionDurations);