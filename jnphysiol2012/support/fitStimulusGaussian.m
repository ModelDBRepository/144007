function [fit_mu, fit_sigma, fit_speed, fit_trace, stim_t] = fitStimulusGaussian(stim_t, stim_vals, extra_t, extra_stim, varargin)
%  [fit_mu, fit_sigma, fit_speed] = fitStimulusGaussian(stim_t, stim_vals, extra_t, extra_stim, varargin)
%
% This function does a least squares gaussian fit on the stimulus luminance sequence, like the linear fit function
% findStimSlope2, but this fits a cumulative gaussian function in order to get an estimate of speed based on the ratios
% of standard deviations of the spatial RF of photoreceptors and the temporal integrated gaussian. 
%
% Varargin: 1 - the width of the receptive field STD in deg. 

sigma_x = 2/3; %the light adapted width of a photoreceptor RF, in deg (Wilson 1974).
if(length(varargin) >= 1)
    sigma_x = varargin{1}; 
end
if(length(varargin) >= 2)
    pb=1;
    plot_ah = varargin{2};
end

tt = cat(1, stim_vals(:), extra_stim(:));
stim_t = cat(1, stim_t(:), extra_t(:));
[stim_t, sorti] = sort(stim_t);
stim_vals = tt(sorti);
if (stim_vals(end) >= stim_vals(1)) % This is for luminance increases
     % for the normcdf function, params(1) is the mean, and params(2) is the std.
     anon_gausscdf = @(params, x) normcdf(x,params(1), params(2)); %normal cumulative gaussian function
else % and for decreases
     anon_gausscdf = @(params, x) -1*normcdf(x,params(1), params(2)) + 1; %flipped cumulative gaussian function
end

mu0 = (max(stim_t) - min(stim_t))/2 + min(stim_t);
sigma0 = (max(stim_t) - min(stim_t))/3;
lb = [-4000 0]; % mean value range wide enough for at least l/v = 80ms looming luminance changes.
ub = [2000 1000];
opts =  optimset('MaxFunEvals', 1000, 'TolFun', 1e-7, 'Display', 'off');


% fit the later section of the curve
between = stim_vals < max(stim_vals) & stim_vals > min(stim_vals); %take the actual luminance change portion
bl = length(between); 
% take the later 2/3 of the luminance change because for many of the slower luminance changes the first 1/3 are 
% much slower than the rest - the final 2/3 yields a better qualitative fit. 
range = (1:ceil(bl*2/3)) + floor(bl/3) + between(1); 
[fit_params,resnorm,resid,ef,output] = lsqcurvefit(anon_gausscdf, [mu0, sigma0], stim_t(range), stim_vals(range), lb, ub, opts);

fit_mu = fit_params(1);
fit_sigma = fit_params(2);
fit_speed = sigma_x/fit_sigma * 1000; %convert to deg/sec from deg/msec
fit_trace = anon_gausscdf(fit_params, stim_t);

% plots the fits
if(exist('plot_ah'))
    line(stim_t, fit_trace, 'Parent', plot_ah, 'Color', 'r');
end
