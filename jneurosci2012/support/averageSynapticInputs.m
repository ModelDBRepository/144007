function [mean_weighted_exc_inputs, mean_weighted_inh_inputs, time_vect] = averageSynapticInputs(synInputFile, time_vect, loverv, pb)
% function [mean_weighted_exc_inputs, mean_weighted_inh_inputs] = averageSynapticInputs(time_vect, pb)
% plot the averaged synaptic inputs, but need to construct vectors for each of them first

load(synInputFile); %this is the location of a file with a saved set of synaptic inputs
plotblue = [.3, .6, 1]; % a nice blue color
colors = {[0 1 0], [1 0 0], plotblue, [0 0 0]};
if isempty(time_vect)
    time_vect = -1850:.05:150; %default time period and sampling
end
if isempty(loverv)
    loverv = [10 40 80];
end
time_vect = time_vect(:);
dt = mean(diff(time_vect));
edges = [time_vect - dt./2; time_vect(end) + dt./2]; %edges of the time bins
lh=[];
ntrials = round(length(synapticInputs)/length(loverv));
weighted_exc_inputs = NaN.*ones(length(time_vect), ntrials, length(loverv));
weighted_inh_inputs = NaN.*ones(length(time_vect), ntrials, length(loverv));
for ii = 1:length(loverv)
    si = find([synapticInputs.loverv] == loverv(ii));
    for jj=1:length(si)
        exc_hist = makeWeightedHist(edges, synapticInputs(si(jj)).exc_times, synapticInputs(si(jj)).exc_gsyn_vec);
        weighted_exc_inputs(:,jj,ii) = makeAlphaSynapseTrace(time_vect, exc_hist, 0, .3);
    end
    for jj=1:length(si)
        inh_hist = makeWeightedHist(edges, synapticInputs(si(jj)).inh_times, synapticInputs(si(jj)).inh_gsyn_vec);
        weighted_inh_inputs(:,jj,ii) = makeAlphaSynapseTrace(time_vect, inh_hist, 0, 3);
    end
end
mean_weighted_exc_inputs = squeeze(nanmean2(weighted_exc_inputs, 2));
mean_weighted_inh_inputs = squeeze(nanmean2(weighted_inh_inputs, 2));

if pb
    figure;
    subplot(2,1,1);
    title('Excitatory Inputs');
    plot(time_vect, mean_weighted_exc_inputs);
    subplot(2,1,2);
    title('Inhibitory Inputs');
    plot(time_vect, mean_weighted_inh_inputs);
end


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
    alpha_onset = .1;
    alpha_tau = 3;
end
dt = mean(diff(t));
alpha_t = 0:dt:(10*(alpha_onset + alpha_tau)); alpha_t = alpha_t';
alpha_func = @(t,gmax, onset, tau) max(0, gmax .* ((t - onset)./tau) .* exp(-(t - onset - tau)./tau));
ay = alpha_func(alpha_t, gmax, alpha_onset, alpha_tau);
synTrace = conv(events, ay);
synTrace = synTrace(1:length(events));
