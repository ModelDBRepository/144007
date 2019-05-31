function synTrace = makeAlphaSynapseTrace(t, events, alpha_onset, alpha_tau)
% synTrace = makeAlphaSynapseTrace(t, events, alpha_onset, alpha_tau)
% ------------------------------------------------------------------------
%
% This function convolves an alpha synapse trace with a binary
% event trace. synTrace is always returned as a column vector.

t = t(:);
events = events(:);
if size(t) ~= size(events)
    error('vectors t and events are not the same size.  They must be the same size');
end

gmax = 1;
if (nargin < 3)
    alpha_onset = 0;
    alpha_tau = 3;
end
dt = mean(diff(t));
alpha_t = 0:dt:(10*(alpha_onset + alpha_tau)); alpha_t = alpha_t';
alpha_func = @(t,gmax, onset, tau) max(0, gmax .* ((t - onset)./tau) .* exp(-(t - onset - tau)./tau));
ay = alpha_func(alpha_t, gmax, alpha_onset, alpha_tau);
synTrace = conv(events, ay);
synTrace = synTrace(1:length(events));
