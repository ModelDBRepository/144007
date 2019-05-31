% --------------------------------------------------------------
function hist_vect = makeWeightedHist(bin_edges, vals, weights)
% ----------------------------------------------------------------
%
% function hist_vect = makeWeightedHistogram(bins, vals, weights)
%
% This is a function that works like a normal histogram function, except
% each of the events in the histogram have a different weight.  I'm using it
% for making a trace of the summed synaptic input over time, given a set of 
% synaptic timings and weights.  
% 
% Inputs and outputs are self-explanatory.  The vals and weights vectors should be the same size.
% The length of hist_vect will be length(bin_edges)-1

hist_vect = zeros(length(bin_edges)-1, 1);
for j=1:length(vals) %for each event
        ii = find(vals(j) >= bin_edges, 1, 'last'); %find the time bin it's in
        if (~isempty(ii) && ii <= length(hist_vect)) 
            hist_vect(ii) = hist_vect(ii) + weights(j); %add the weight
        end 
end