function [opt_delay, peak_corr] = findOptimalDelay(sig1, sig2, range, dt, varargin)
% function [opt_delay, peak_corr] = findOptimalDelay(t, sig1, sig2)
% 
% This function finds the optimal delay between the two input signals to maximize
% the correlation between the two signals. This means that it finds the delay of peak 
% of the cross-correlogram. The returned delay is the amount of time you need to shift signal
% 2 to achieve a maximal correlation with signal 1.  This, if signal 1 lags signal 2, the function returns
% a positive delay.
% Sig1 and 2 are the two signals.  Range is the time range of measured correlations, and
% dt is the sampling interval of the signals.

if (length(varargin) >= 1)
    pb=varargin{1};
else
    pb=0;
end
if (length(varargin) >= 2)
    tvect=varargin{2};
else
    tvect=1:length(sig1);
end


%ww = randn(1000,1);
%baseline subtract and normalize
sig1 = sig1 - sig1(1); sig2 = sig2 - sig2(1);
sig1 = sig1/sum(sig1); sig2 = sig2/sum(sig2);

[c12,lags12] = xcorr(sig1,sig2,'coeff');
t_lags12 = lags12*dt;
% autocorrelation
[c11,lags11] = xcorr(sig1,sig1,'coeff');
t_lags11 = lags11*dt;
[c22,lags22] = xcorr(sig2,sig2,'coeff');
t_lags22 = lags22*dt;
% find the xcorr peak
[max_c12, maxi] = max(c12);
max_lag = t_lags12(maxi);
% assign outputs
opt_delay = max_lag;
peak_corr = max_c12;

if pb
    figure; hold on;
    ph = plot(tvect, sig1); set(ph, 'Color', 'b');
    t2 = (1:length(sig2))*dt + tvect(1);
    ph = plot(t2, sig2); set(ph, 'Color', 'g');
    ph = plot(t2 + max_lag, sig2, 'g--');
    legend('Signal 1', 'Signal 2'); 
    figure; hold on;
    plot([0 0], [0 1], 'k');
    plot(t_lags12,c12, 'r');
    plot(t_lags11, c11, 'Color', [.5 .5 .5]);
    plot(t_lags22, c22, 'Color', [.5 .5 .5], 'LineStyle', ':'); 
   
    legend('Zero Lag', 'X Corr', 'Auto Corr 1', 'Auto Corr 2');
    xlabel('\tau (ms)', 'FontSize', 13);
    ylabel('Correlation Coefficient', 'FontSize', 13);
end
