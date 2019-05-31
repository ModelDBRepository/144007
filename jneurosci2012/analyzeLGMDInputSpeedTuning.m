function analyzeLGMDInputSpeedTuning(sim_results_inh, sim_results_noinh)
% This is a function/script to run the analysis of LGMD speed tuning
% written by Peter W Jones in 2011.

% These are the names of the simulation files that need to be generated in
% order to run this analysis. Generate them using generate_variability_sims.m, 
% run them using NEURON, process using analyze_model.m, and that will produce these
% mat files. Name the results what you want, put those names here. In the
% function, these are given as inputs, but if you comment out the function
% line, use these to define your filenames
 %sim_results_noinh = 'noinh/loom_snr_5_jit_6_exc_0.6_inh_0.0095_analyzed.mat';
 %sim_results_inh = '0319/loom_snr_5_jit_6_exc_0.6_inh_0.0095_analyzed.mat';

% These files hold the saved synaptic input timings from simulation generation
 synapticInputFile = 'synapticInputs.mat';
 avg_synapticInputFile = 'averageSynInput.mat';

speedTuningModelSetup; %script that just offloads some setup tasks.

% The linear fit parameters for the peak timing of the LGMD response to looming,
% which peaks a fixed time after the looming stimulus reaches and angular threshold.
% The values used here will determine what the exact angular threshold and fixed delay are. 
% Those things in turn determine the stimulus velocities at which the angular threshold
% is reached...a value that is plotted extensively as a reference throughout and used as
% a maximum for response fitting.  What we have done is analyzed the model output beforehand
% to get the fit parameters, then entered them here. The LGMD model itself contains noise, thus
% there is no guarantee that future model runs will result in exactly the same fit parameters.
% lgmd_peak_time_line = [4.4 18.7]; %my data from Jones and Gabbiani (2010), for reference
lgmd_peak_time_line = [4.8 14.7]; %the LGMD model's peak firing rate fit

bin_width = 2; %in ms

%this defines the gaussian filter to use to filter timecourses
filt_width = 8/bin_width; % 8ms std in samples
filtx = -3*filt_width:1:3*filt_width;
filty = my_normpdf(filtx,0,filt_width);
filty = filty/sum(filty);

pre = 2000; %ms - time before the stimulus
post = 300; %ms - time after collision or stim reaches max size
options = optimset('TolFun', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'off');
powfun = @(p,x)p(2).*(rectify(x).^p(1)); %power law function for spike threshold, p(1):exponent, p(2):scale
peak_loverv_fun = @(p,x) polyval([p(1) -p(2)], x); % linear function but the sign of the offset is switched for the peak time (alpha*l/v - delta)

for i = 1:length(mparams)
     
    % Going to make histograms of the input times
    edges = (min(mparams(i).mov_t)-pre):bin_width:(max(mparams(i).mov_t)+post);
    centers = edges(1:end-1)+(bin_width/2); 
    
    % Get the response latencies and magnitudes for each facet based on the luminance change
    % information present in mparams
    lats = transition2responseLatency(mparams(i).transition_durations, 'lmc', 'monitor');
    lmc_times = mparams(i).transition_start_times' + lats;
    lmc_mags = speed2lgmdResponse(mparams(i).RFspeeds, 'lmc');
    lats = transition2responseLatency(mparams(i).transition_durations, 'vc', 'monitor');
    vc_times = mparams(i).transition_start_times' + lats;
    vc_mags = speed2lgmdResponse(mparams(i).RFspeeds, 'vc');
    lats = transition2responseLatency(mparams(i).transition_durations, 'cc', 'monitor');
    cc_times = mparams(i).transition_start_times' + lats;
    cc_mags = speed2lgmdResponse(mparams(i).RFspeeds, 'cc');

    % makes timecourses based on the weight and timings of the inputs
    lmc_weighted = makeWeightedHist(edges, lmc_times, lmc_mags);
    vc_weighted = makeWeightedHist(edges, vc_times, vc_mags);
    cc_weighted = makeWeightedHist(edges, cc_times, cc_mags);
    % filter them to smooth some.
    lmc_weighted = conv2(lmc_weighted(:), filty(:), 'same');
    vc_weighted = conv2(vc_weighted(:), filty(:), 'same');
    cc_weighted = conv2(cc_weighted(:), filty(:), 'same');
    
    min_t = -400;
    max_t = 100;

    peak_time = -1*peak_loverv_fun(lgmd_peak_time_line, abs(mparams(i).loverv)); %time of the LGMD firing peak
    thresh_time(i) = peak_time - lgmd_peak_time_line(2); %time of corresponding stim ang thresh 
    max_val = max([max(lmc_weighted), max(vc_weighted), max(cc_weighted)]);
    %stimulus timecourse, in terms of angle and velocity
    t = centers; %using the centers of the time bins
    pz = (t<=0); %logical array to select pre-collision timepoints
    max_theta = 41; %maximum half-angle
    theta = -atan(mparams(i).loverv./t)*180/pi; %this is the half angle.
    stoppedi = find(theta >= max_theta, 1, 'first');
    theta(stoppedi:end) = max_theta;
    vel_func = @(loverv, t) -loverv./(t.^2+loverv^2)*180/pi .* 1000; %gives deg/s
    vel = abs(vel_func(mparams(i).loverv, t));%deg/s
    vel(stoppedi:end) = 0; %velocity is zero when stimulus is stopped
    vel_at_lgmd_peak(i) = abs(vel_func(mparams(i).loverv, peak_time)); %velocity at the lgmd peak
    vel_at_ang_thresh(i) = abs(vel_func(mparams(i).loverv, thresh_time(i))); %velocity at the ang threshold
    pz = (t<=0 & theta < max_theta); %logical array to select pre-collision timepoints
    
    %%%%%%%%%%%%%%% Begin plotting %%%%%%%%%%%%%%%%%%%%%%%%%
    % This figure pops up anew for each l/v value, and looks at the responses through the pathway for each.
    
    % This first subplot plots the relative weights of individual inputs through time.
    fh5(i) = figure;
    subplot(3,2,5); %third row, plot the CC responses
    plot(cc_times, cc_mags, '.k'); hold on;
    xlim([min_t, max_t]); ylim([0 1]); 
    title('Current Clamp - LGMD');
    xlabel('time (ms)');
    ylabel('Single facet Input strength');
    set(gca, 'TickDir', 'out'); 
    %plot the VC responses, indicative of medulla input
    plot( vc_times,  vc_mags, '.c');
    xlim([min_t, max_t]); ylim([0 1]);
    title('Voltage Clamp - Medulla');
    set(gca, 'TickDir', 'out');
    % Now lmc responses
    plot(lmc_times, lmc_mags, '.r'); hold on;
    xlim([min_t, max_t]); ylim([0 1]);
    title(['LMC - Response Strength over Time, l/v=' num2str(mparams(i).loverv)]);
    set(gca, 'TickDir', 'out');
  
    % This plots the overall response at each stage through time.
    subplot(3,2,4); hold on;
    plot(centers, lmc_weighted(:), 'r-', 'Linewidth', 1);
    plot(centers, vc_weighted(:), 'c-', 'Linewidth', 1);
    plot(centers, cc_weighted(:), 'k-', 'Linewidth', 1);
    plot([thresh_time(i) thresh_time(i)], [0 200], 'k--');
    ymax = max([max(lmc_weighted(:)), max(vc_weighted(:)), max(cc_weighted(:))]);
    xlim([min_t, max_t]); ylim([0 ymax]);
    xlabel('Time (ms)');
    ylabel('Input Strength');
    set(gca, 'TickDir', 'out');
    
    % Plotting the stimulus velocity
    subplot(3,2, 2); hold on;
    xlim([min_t, max_t]);
    plot(t(pz),vel(pz),'k', 'LineWidth', 2);
    plot([thresh_time(i) thresh_time(i)], [0 vel_at_ang_thresh(i)],'k--'); plot([min_t, thresh_time(i)], [vel_at_ang_thresh(i) vel_at_ang_thresh(i)], 'k--');
    xlabel('Time (ms)'); ylabel('Ang Velocity (deg/s');
    ylim([0 vel_at_ang_thresh(i) + 100]);
    ylim([0 500]);
    setPresentationDefaults(fh5(i), 0);
    
    %save across l/v so that we can compare them between
    resp(i).vel = vel(:);
    resp(i).lmc_weighted = lmc_weighted(:);
    resp(i).vc_weighted = vc_weighted(:);
    resp(i).cc_weighted = cc_weighted(:);
    resp(i).theta = theta(:);
    resp(i).mov_t = centers(:);

end

%% Loading versions of the LGMD model results without inhibition
% Need to do this now in order to use the stimulus time vector information, etc.
load(sim_results_noinh); %model results
noinh_stim = stim;
clear stim;

%%  This generates some proper stimulus kinematic vectors for the model results.
theta_fun = @(loverv, t) 2 * atan(-loverv./t)*180/pi; %theta full-size
max_theta = 82; %this is the max because it is limited by the screen size we have in the lab
figure; hold ;
for ii= 1:length(noinh_stim)   
    % make a velocity vector that matches the timebase of the synaptic/model time
    v = vel_func(-noinh_stim(ii).loverv, noinh_stim(ii).tvec);
    ang_size = theta_fun(noinh_stim(ii).loverv,  noinh_stim(ii).tvec);
    ang_size(ang_size > max_theta | ang_size < 0) = max_theta;
    endt = (noinh_stim(ii).tvec >= 0 | ang_size == max_theta);  
    v(endt) = 0;
    theta_max_t = noinh_stim(ii).tvec(find(endt == 1, 1, 'first'));
    noinh_stim(ii).vel =  v;
    noinh_stim(ii).theta = ang_size;
    plot(noinh_stim(ii).tvec, noinh_stim(ii).vel, 'Color', colors{ii}); hold on;
end
% make sure that we have the proper mean membrane potential
% We don't have the mean proximal Vm that we need, so let's build it.
for ii = 1:length(noinh_stim)    
    vm_prox_m = zeros(length(noinh_stim(ii).mu_vmfilt),length(noinh_stim(ii).trial));
    conv_inst_freq_m = zeros(length(noinh_stim(ii).mu_vmfilt),length(noinh_stim(ii).trial));
    for jj = 1:length(noinh_stim(ii).trial)
        vm_prox_m(:,jj) = noinh_stim(ii).trial(jj).vm_filt_prox;
        conv_inst_freq_m(:,jj) = noinh_stim(ii).trial(jj).conv_inst_freq;
    end
    noinh_stim(ii).mu_vmprox_filt = mean(vm_prox_m,2);
    noinh_stim(ii).a_vmprox_filt = vm_prox_m;
    noinh_stim(ii).a_conv_inst_freq = conv_inst_freq_m; 
end

%%%%%%%%%%%%%%% TIMING COMPENSATION SCHEME  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So, in order to look at the speed tuning in a way that makes sense, we want to compensate for neural delays
% between processing stages in order to get the relationship between the population signal and the stimulus 
% velocity that triggered it.  We are doing this by offsetting the signals by the amount which produces a maximum correlation 
% between the two.
ifr_dt = mean(diff(noinh_stim(i).tvec));
input_dt = mean(diff(resp(i).mov_t));
% The goal here is to compute the optimal delay directly by finding where it correlates best with the velocity signal. 
% In practice, the result of this is that the offsets found are similar to those introduced by the delays from the data 
% (single facet linear relationship with luminance change duration).
for i = 1:length(noinh_stim)
    vel = resp(i).vel;
    % LGMD Vm data
    [best_off, peak_corr] = findOptimalDelay(resp(i).cc_weighted, vel, [], input_dt, 0)
    offset_t(i).cc = best_off;
    % LGMD Im data
    [best_off, peak_corr] = findOptimalDelay(resp(i).vc_weighted, vel, [], input_dt, 0)
    offset_t(i).vc = best_off;
    % LMC signals
    [best_off, peak_corr] = findOptimalDelay(resp(i).lmc_weighted, vel, [], input_dt, 0)
    offset_t(i).lmc = best_off;
    % LGMD model IFR
    [best_off, peak_corr] = findOptimalDelay(noinh_stim(i).mu_conv_inst_freq, noinh_stim(i).vel, [], ifr_dt, 0, noinh_stim(i).tvec)
    offset_t(i).ifr = best_off;
    
    % calculate the offset in samples
    offset(i).cc = round(mean(offset_t(i).cc)/input_dt);
    offset(i).vc = round(mean(offset_t(i).vc)/input_dt);
    offset(i).lmc = round(mean(offset_t(i).lmc)/input_dt);
    offset(i).ifr = round(mean(offset_t(i).ifr)/ifr_dt);
end
%these are the offset values in index offsets, means across l/v.
mean_offset.cc = round(mean([offset_t.cc])/input_dt);
mean_offset.vc = round(mean([offset_t.vc])/input_dt);
mean_offset.lmc = round(mean([offset_t.lmc])/input_dt);
mean_offset.ifr = round(mean([offset_t.ifr])/ifr_dt);
mean_offset_t.cc = (mean([offset_t.cc]));
mean_offset_t.vc = (mean([offset_t.vc]));
mean_offset_t.lmc = (mean([offset_t.lmc]));
mean_offset_t.ifr = (mean([offset_t.ifr]));

disp(mean_offset_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Revisiting the past figures to fill in the speed tuning curves 
for i=1:length(resp)
    figure(fh5(i));
    subplot(3,2,6); hold on;
    pz = find(resp(i).mov_t<=0 & resp(i).theta < max_theta & resp(i).vel > 2);
    nfacets = length(mparams(i).transition_start_times);
    plot(resp(i).vel(pz), resp(i).lmc_weighted(pz+mean_offset.lmc), 'r-', 'Linewidth', 1);
    plot(resp(i).vel(pz), resp(i).vc_weighted(pz+mean_offset.vc), 'c-', 'Linewidth', 1);
    plot(resp(i).vel(pz), resp(i).cc_weighted(pz+mean_offset.cc), 'k-', 'Linewidth', 1);
    plot([vel_at_ang_thresh(i) vel_at_ang_thresh(i)], [0 20], 'k--');
    ylim([0 max([resp(i).lmc_weighted(pz); resp(i).vc_weighted(pz); resp(i).cc_weighted(pz)])]);
    ylim([0 5]);
    xlim([0 vel_at_ang_thresh(i) + 40]);
    ylabel('Population response strength')
    xlabel('Stimulus Angular Speed (deg/s)');
    legend('LMC', 'LGMD Im', 'LGMD Vm');
    setPresentationDefaults(gcf, 0);
end

%% want to normalize by a value that makes sense for each stage, so that we can compare across stages while still looking across l/v values
if (1 == 1)
    for i = 1:length(mparams)
        normi = find(resp(i).vel >= vel_at_ang_thresh(i), 1, 'first');
        lmc_normval(i) = resp(i).lmc_weighted(normi + mean_offset.lmc);
        vc_normval(i) = resp(i).vc_weighted(normi + mean_offset.vc);
        cc_normval(i) = resp(i).cc_weighted(normi + mean_offset.cc);
    end
    for i = 1:length(mparams)
        resp(i).lmc_weighted_norm = resp(i).lmc_weighted ./ lmc_normval(i);
        resp(i).vc_weighted_norm = resp(i).vc_weighted ./ vc_normval(i);
        resp(i).cc_weighted_norm = resp(i).cc_weighted ./ cc_normval(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots Fig. 5C,D - The angular velocity tuning curves for the populations giving input to the LGMD.
% It does this plotting different l/v values on the same axis. Also, we fit the power-law functions
% in this loop
figure;
colors = {[0 1 0], [1 0 0], plotblue, [0 0 0]};
lmc_fitp = []; vc_fitp = []; cc_fitp = [];
max_y = 50;

for i=1:length(mparams)
    disp(['l/|v| value: ' num2str(mparams(i).loverv)]);
    color_light{i} = (colors{i} + [1 1 1])./2;
    edges = resp(i).mov_t - (bin_width/2);
    centers = resp(i).mov_t;
    t = centers;
    theta = atan(mparams(i).loverv./t)*180/pi;
    theta = resp(i).theta;
    nonzero = resp(i).lmc_weighted_norm > 0; %the bins where there is response
    pz = find(theta < max_theta & t < 0 & nonzero); % just a set of indices to plot that doesn't overrun the useful bounds
    % we are getting indices to define the bounds of plotting and fitting 
    threshi = find(resp(i).vel >= vel_at_ang_thresh(i), 1, 'first');
    maxi = find(resp(i).vel >= vel_at_ang_thresh(i), 1, 'first');
    mini = find(resp(i).vel >= 2, 1, 'first');
    fiti = mini:threshi;
    ploti = mini:maxi;
    
    %%%%%%%%%%%%%%  LMC section  %%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,1);  hold on;
    [max_val, maxi] = nanmax2(resp(i).lmc_weighted);
    ploti = 1:(maxi-mean_offset.lmc);
    % actual tuning function
    plot(resp(i).vel(ploti), resp(i).lmc_weighted(ploti+mean_offset.lmc), 'Color', colors{i}, 'linestyle', '-', 'LineWidth',1);
    plot([vel_at_ang_thresh(i) vel_at_ang_thresh(i)], [0 100], '-', 'Color', colors{i}); %vertical line
    resp_at_ang_thresh(i).lmc = resp(i).lmc_weighted(threshi+mean_offset.lmc);
    % Now going to fit a power law function for the section till the speed that the angular threshold was reached
    xd = resp(i).vel(fiti); %vectors to fit to
    yd = resp(i).lmc_weighted(fiti+mean_offset.lmc);
    lmc_fitp(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e3 1e3], options);
    yfit2 = powfun(lmc_fitp(i,:), xd);
    lmc_pow_r2(i) = calcR2(yd, yfit2)
    plot(xd, yfit2, '--', 'Color', color_light{i}, 'LineWidth',2);
    title('LMC - Speed tuning to instantaneous stimulus velocity');
    ylim([0 max_y]);  xlim([0 max(vel_at_ang_thresh)+20]);
    
    %%%%%%%%%%%%%% LGMD Voltage Clamp (Im derived) section %%%%%%%%%%%%%%%%%%%%
    subplot(3,1,2);  hold on;
    [max_val, maxi] = nanmax2(resp(i).vc_weighted);
    ploti = 1:(maxi-mean_offset.vc);
    di = mean_offset.vc;
    plot(resp(i).vel(ploti), resp(i).vc_weighted(ploti+di), 'Color', colors{i}, 'linestyle', '-', 'LineWidth',1);
    plot([vel_at_ang_thresh(i) vel_at_ang_thresh(i)], [0 100], '-', 'Color', colors{i});
    resp_at_ang_thresh(i).vc = resp(i).vc_weighted(threshi+mean_offset.vc);
    % Now going to fit a power law function for the section till the speed that the angular threshold was reached
    xd = resp(i).vel(fiti); %vectors to fit to
    yd = resp(i).vc_weighted(fiti+mean_offset.vc);
    yd = resp(i).vc_weighted(fiti+di);
    vc_fitp(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e3 1e3], options);
    yfit2 = powfun(vc_fitp(i,:), xd);
    vc_pow_r2(i) = calcR2(yd, yfit2)
    title('LGMD Im');
    ylim([0 max_y]);  xlim([0 max(vel_at_ang_thresh)+20]);
    plot(xd, yfit2, '--', 'Color', color_light{i}, 'LineWidth',2);
    
    %%%%%%%%%%%%%%%%% LGMD Current clamp (Vm derived) section %%%%%%%%%%%%%%%%%%
    subplot(3,1,3); hold on;
    [max_val, maxi] = nanmax2(resp(i).cc_weighted);
    ploti = 1:(maxi-mean_offset.cc);
    plot(resp(i).vel(ploti), resp(i).cc_weighted(ploti+mean_offset.cc), 'Color', colors{i}, 'linestyle', '-', 'LineWidth',1);
    plot([vel_at_ang_thresh(i) vel_at_ang_thresh(i)], [0 100], '-', 'Color', colors{i});
    resp_at_ang_thresh(i).cc = resp(i).cc_weighted(threshi+mean_offset.cc);
    % Now going to fit a power law function for the section till the speed that the angular threshold was reached
    xd = resp(i).vel(fiti); %vectors to fit to
    yd = resp(i).cc_weighted(fiti+mean_offset.cc);
    cc_fitp(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e3 1e3], options);
    yfit2 = powfun(cc_fitp(i,:), xd);
    cc_pow_r2(i) = calcR2(yd, yfit2)
    plot(xd, yfit2, '--', 'Color', color_light{i}, 'LineWidth',1);
    
    ylim([0 max_y]); xlim([0 max(vel_at_ang_thresh)+20]);
    title('LGMD Vm');
    xlabel('Stimulus Velocity (deg/s)');
    ylabel('Popultation Response Strength');
    legendstr{3*(i-1)+1} = ['l/v= ' num2str(mparams(i).loverv) ' ms'];
    legendstr{3*(i-1)+2} = ['Latency Adjusted'];
    legendstr{3*(i-1)+3} = ['LGMD peak'];
    legend(legendstr);
    setPresentationDefaults(gcf, 0);
    
    %set the y axis limits
    stages = {'lmc', 'vc', 'cc'}; %for LMC, LGMD Im, and LGMD Vm
    for j = 1:3
        subplot(3,1,j);
        slowi = find(resp(i).vel(pz) <= 120);
        yv = eval(['resp(i).' stages{j} '_weighted(pz+mean_offset.' stages{j} ');']);
        ymax = max(yv(slowi));
        ylim([0 ymax]);
        ylim([0 50]);
    end
end
%% print the power law coeffecients for the paper - just transform the fit ones. The power law function in the 
% paper is scaled by the max speed that it's fit to, so we need to transform the fit scaling parameter  
for i=1:length(resp)
    lmc_fit_a(i) = vel_at_ang_thresh(i)^lmc_fitp(i,1) * lmc_fitp(i,2);
    vc_fit_a(i) = vel_at_ang_thresh(i)^vc_fitp(i,1) * vc_fitp(i,2);
    cc_fit_a(i) = vel_at_ang_thresh(i)^cc_fitp(i,1) * cc_fitp(i,2);
end
lmc_fit_a
vc_fit_a
cc_fit_a

%% Let's do the same data in subplots, but plot stages on top of each other - see how they compare
% These plots are not in the paper.
figure; hold on;
lmc_fitp3_norm = []; vc_fitp3_norm= []; cc_fitp3_norm= [];
for i = 1:length(resp)
    centers = resp(i).mov_t;
    t = centers;
    theta = atan(mparams(i).loverv./t)*180/pi;
    nonzero = resp(i).lmc_weighted_norm > 0; %the bins where there is response
    pz = find(theta < max_theta & t < 0 & nonzero); % just a set of indices to plot that doesn't overrun the useful bounds
    max_offset = max([mean_offset.lmc, mean_offset.vc, mean_offset.cc]);
    threshi = find(resp(i).vel >= vel_at_ang_thresh(i), 1, 'first');
    maxi = find(resp(i).vel >= max(120, vel_at_ang_thresh(i)), 1, 'first');
    mini = find(resp(i).vel >= 2, 1, 'first');
    fiti = mini:threshi;
    
    subplot(length(resp),1,i); hold on;
    title(['l/v = ' num2str(mparams(i).loverv)]);
    
    % LMC
    [maxv, maxi] = nanmax2(resp(i).lmc_weighted_norm); pz = 1:(maxi-max_offset);
    ph(1) = plot(resp(i).vel(pz), resp(i).lmc_weighted_norm(pz+mean_offset.lmc), 'Color', 'r', 'linestyle', '-', 'LineWidth',2);
    xd = resp(i).vel(fiti);  yd = resp(i).lmc_weighted_norm(fiti+mean_offset.lmc);  %vectors to fit to
    lmc_fitp3_norm(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e4 1e4], options);
    yfit = powfun(lmc_fitp3_norm(i,:), xd);
    plot(xd, yfit, 'Color','r', 'linestyle', '--', 'Linewidth', 2);
    % LGMD Im (VC)
    [maxv, maxi] = nanmax2(resp(i).vc_weighted_norm); pz = 1:(maxi-max_offset);
    ph(2) = plot(resp(i).vel(pz), resp(i).vc_weighted_norm(pz+mean_offset.vc), 'Color', 'c', 'linestyle', '-', 'LineWidth',2);
    yd = resp(i).vc_weighted_norm(fiti+mean_offset.vc);  %vectors to fit to
    vc_fitp3_norm(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e4 1e4], options);
    yfit = powfun(vc_fitp3_norm(i,:), xd);
    plot(xd,yfit, 'Color','c', 'linestyle', '--', 'Linewidth', 2);
    % LGMD Vm (CC)
    [maxv, maxi] = nanmax2(resp(i).cc_weighted_norm); pz = 1:(maxi-max_offset);
    ph(3) = plot(resp(i).vel(pz), resp(i).cc_weighted_norm(pz+mean_offset.cc), 'Color', 'k', 'linestyle', '-', 'LineWidth',2);
    yd = resp(i).cc_weighted_norm(fiti+mean_offset.cc);  %vectors to fit to
    cc_fitp3_norm(i,:) = lsqcurvefit(powfun, [1 1], xd, yd, [0 0], [1e4 1e4], options);
    yfit = powfun(cc_fitp3_norm(i,:), xd);
    plot(xd, yfit, 'Color','k', 'linestyle', '--', 'Linewidth', 2);
   
    %set the limits
    xlim([0 max(vel_at_ang_thresh)+20]);
    slowi = find(resp(i).vel(pz) <= 120);
    ccy = resp(i).cc_weighted_norm(pz+mean_offset.cc); ymax = max(ccy(slowi));
    ylim([0 maxv]);
    plot([vel_at_ang_thresh(i) vel_at_ang_thresh(i)], [0 max_val], '-', 'Color', 'k');
end
setPresentationDefaults(gcf, 0);

%% Loading versions of the LGMD model results run with and without inhibition
% we already loaded the  no inhibition file above
[noinh_stim, noinh_vtune] = loadModelResults(noinh_stim, mean_offset_t); %does some supplementary processing

load(sim_results_inh); %model w/ inhibition results
[reg_stim, reg_vtune] = loadModelResults(stim, mean_offset_t);
clear stim; % save memory

% Also, read the synaptic inputs from many trials and average them.
% This function reads a file written when the simulations are generated that saves the individual synaptic
% inputs.  This is necessary because it is cumbersome to go back and parse all of the hoc files in order to 
% get the synaptic input.  Much easier to read in a mat file.
if ~exist('mean_exc_input','var')
    if exist(avg_synapticInputFile, 'file')
        load(avg_synapticInputFile);
    else
        [mean_exc_input, mean_inh_input, synInputTime] = averageSynapticInputs(synapticInputFile, [],[],0); %#ok<ASGLU>
        save(avg_synapticInputFile, 'mean_exc_input', 'mean_inh_input', 'synInputTime');
    end
end
exc_input_max = max(max(mean_exc_input));


%% Figure 6: Plotting speed tuning functions for model Vm and IFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options for function fitting
options = optimset('TolFun', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'off', 'algorithm',  'trust-region-reflective');
powfun = @(p,x)p(2).*(rectify(x).^p(1)); %power law function for spike threshold, p(1):exponent, p(2):scale
vel_params_ifr3 = zeros(length(reg_stim),3); vel_params_vm3 = zeros(length(reg_stim),3);
figure; hold on;
for ii=1:length(reg_stim)
    lb = [0 0 0]; %for fitting Weibull functions
    ub = [1e1 1e4 1e4];
    %finding the range over which to fit - to the velocity at the angular threshold
    maxi = find(noinh_vtune.vm.vel_m(:,ii) >= vel_at_ang_thresh(ii), 1, 'first');
    if isempty(maxi)
        maxi = length(noinh_vtune.vm.vel_m(:,ii));
    end
    mini = find(noinh_vtune.vm.vel_m(:,ii) >= 2, 1, 'first');
    lini = mini:maxi;
    nz = noinh_vtune.vm.vel_m(:,ii) > 0;
    
    %%%%%%%%%%%%%%%%%    Vm    %%%%%%%%%%%%%%%%%%
    subplot(2,1,1); hold on;
    ph2(ii) = plot(noinh_vtune.vm.vel_m(nz,ii), noinh_vtune.vm.norm_resp_m(nz,ii), ...
                    'color', colors{ii}, 'LineWidth',2, 'Linestyle', '-');
    % 3 parameter weibull function fitting
    wbl_p = lsqcurvefit(@weibull_func3, [1 20 1], noinh_vtune.vm.vel_m(lini,ii), ...
                        noinh_vtune.vm.norm_resp_m(lini,ii), lb, ub, options);
    vel_params_vm3(ii, :) = wbl_p;
    vmfit_y = weibull_func3(wbl_p, noinh_vtune.vm.vel_m(lini,ii));
    plot(noinh_vtune.vm.vel_m(lini,ii), vmfit_y, 'Color', color_light{ii}, 'LineWidth', 3, 'LineStyle','--');
    vm_r2(ii) = calcR2(noinh_vtune.vm.norm_resp_m(lini,ii), vmfit_y) %calculate a goodness of fit measure
    ts = ['\alpha = ' num2str(wbl_p(1)) '  \lambda = ' num2str(wbl_p(2)) '  \kappa = ' num2str(wbl_p(3))]
    th= text(60, .1*ii, ts); %parameter printing 
    set(th, 'Color', colors{ii}, 'FontSize', 12);
    plot(1:450, powfun(cc_fitp3_norm(ii,:), 1:450), 'LineStyle', ':', 'LineWidth',2, 'Color', color_light{ii});
    % vertical line giving peak velocity
    plot([vel_at_ang_thresh(ii) vel_at_ang_thresh(ii)], [0 10], '-', 'Color', colors{ii}); 
    xlim([0 max(vel_at_ang_thresh)+20]);
    
    %%%%%%%%%%%%%  Firing Rates - IFR   %%%%%%%%%%%%%%%%%
    subplot(2,1,2); hold on; %top plot for Vm, and the bottom for the IFR
    % We need to find the range again since the indices are different
    maxi = find(noinh_vtune.ifr.vel_m(:,ii) >= vel_at_ang_thresh(ii), 1, 'first');
    if isempty(maxi)
        maxi = length(noinh_vtune.ifr.vel_m(:,ii));
    end
    mini = find(noinh_vtune.ifr.vel_m(:,ii) >= 2, 1, 'first');
    lini = mini:maxi;
    nz = noinh_vtune.ifr.vel_m(:,ii) > 0;
    
    plot(noinh_vtune.ifr.vel_m(nz,ii), noinh_vtune.ifr.norm_resp_m(nz,ii), 'color', colors{ii}, 'LineWidth',2, 'Linestyle', '-');
    % 3 parameter weibull function fitting
    wbl_p = lsqcurvefit(@weibull_func3, [1 50 2], noinh_vtune.ifr.vel_m(lini,ii), noinh_vtune.ifr.norm_resp_m(lini,ii), lb, ub, options);
    vel_params_ifr3(ii, :) = wbl_p;
    fit_y = weibull_func3(wbl_p, noinh_vtune.ifr.vel_m(lini,ii));
    plot(noinh_vtune.ifr.vel_m(lini,ii), fit_y, 'Color', color_light{ii}, 'LineWidth', 3, 'LineStyle','--');
    ifr_r2(ii) = calcR2(noinh_vtune.ifr.norm_resp_m(lini,ii), fit_y)
    ts = ['\alpha = ' num2str(wbl_p(1)) '  \lambda = ' num2str(wbl_p(2)) '  \kappa = ' num2str(wbl_p(3))]
    th= text(60, .1*ii, ts); %parameter printing 
    set(th, 'Color', colors{ii}, 'FontSize', 12);
    plot([vel_at_ang_thresh(ii) vel_at_ang_thresh(ii)], [0 10], '-', 'Color', colors{ii}); % vertical line giving peak velocity
    xlim([0 max(vel_at_ang_thresh)+20]);
    
end
% Just cleaning up the plots a little
subplot(2,1,1);
xlim([0 max(vel_at_ang_thresh)+20]);
ylim([0 1]);
xlabel('Stimulus velocity (deg/sec)', 'FontSize', 14);
ylabel('Normalized Model Vm', 'FontSize', 14);
leg_str = {'l/v = 10', 'l/v=40', 'l/v=80'};
legend(ph2, leg_str, 'FontSize', 12);

subplot(2,1,2);
xlim([0 max(vel_at_ang_thresh)+20]);
ylim([0 1]);
xlabel('Stimulus velocity (deg/sec)', 'FontSize', 14);
ylabel('Normalized Model IFR', 'FontSize', 14);
setPresentationDefaults(gcf, 0);

%% Figure 7 - What we want to do is lay out methodically the transformation happening from stimulus variables to LGMD output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Want to plot the membrane Vm the near SIZ along with the aggregate synaptic input and the dendritic Vm. 
% What this actually shows is that the input current near the SIZ rises much less sharply than the linear sum
% of the inputs. This is Figure 7C.
figure; hold on;
vm_adj = zeros(length(noinh_stim(1).mu_vmfilt), length(noinh_stim));
vm_dend_adj = zeros(length(noinh_stim(1).mu_vmfilt), length(noinh_stim));
% This is near the SIZ - adjusted for the resting Vm.
for ii = 1:length(noinh_stim)
    vm_adj(:,ii) = noinh_stim(ii).mu_vmprox_filt - mean(noinh_stim(ii).mu_vmprox_filt(200:500));
end
% This is the dendritic Vm halfway along one of the central branches
for ii = 1:length(noinh_stim)
    vm_dend_adj(:,ii) = noinh_stim(ii).mu_vm_dend_horiz - mean(noinh_stim(ii).mu_vm_dend_horiz(200:500));
end
max_vm = max(max(vm_adj)); %first find the max 
max_vm_dend = max(max(vm_dend_adj));
% plotting model outputs of various kinds
for ii=1:length(noinh_stim)
    plot(thresh_time(ii)*ones(1,2), [-.4 1], 'Color', colors{ii});
    ph = plot(noinh_stim(ii).tvec, vm_adj(:,ii)./max_vm, 'Color', colors{ii}, 'LineWidth', 2); hold on
    ph = plot(noinh_stim(ii).tvec, vm_dend_adj(:,ii)./max_vm_dend, 'Color', color_light{ii}, 'LineWidth', 1, 'LineStyle', '-'); hold on
end
% Now onto the model inputs
max_cc = 0;
for ii=1:length(resp)
    max_cc = max(max_cc, max(resp(ii).cc_weighted));
end
for ii = 1:length(resp)
    color_dark{ii} = (colors{ii})./1.5;
    plot(synInputTime, mean_exc_input(:,ii) ./ exc_input_max, 'Color', color_dark{ii}, 'LineWidth', 2, 'LineStyle', '--'); hold on;
end
xlim([-800 200]);
setPresentationDefaults(gcf, 0);
ylabel('Normalized Signal', 'FontSize', 16); xlabel('Time from Collision (ms)', 'FontSize', 16);

%% Now, the velocity and the inputs, Fig 7A
theta_fun = @(loverv, t) 2 * atan(-loverv./t)*180/pi; %theta full-size
max_theta = 82; %this is the max because it is limited by the screen size we have in the lab
timefig = figure;
hold on;
%subplot(2,1,1); hold on;
syn_input_v = zeros(length(synInputTime), length(reg_stim));
abs_max_in = max(max(mean_exc_input));
for ii=1:length(reg_stim)
    % make a velocity vector that matches the timebase of the synaptic time
    v = vel_func(-reg_stim(ii).loverv, synInputTime);
    ang_size = theta_fun(reg_stim(ii).loverv, synInputTime);
    ang_size(ang_size > max_theta | ang_size < 0) = max_theta;
    endt = (synInputTime >= 0 | ang_size == max_theta);  v(endt) = 0;
    %theta_max_t = reg_stim(ii).tvec(find(endt == 1, 1, 'first'))
    syn_input_v(:,ii) = v; 
end
abs_max_vel = max(max(syn_input_v));
for ii=1:length(reg_stim) %loop for plotting
    v = syn_input_v(:,ii);
    si = find(synInputTime >= -800, 1, 'first');
    [max_vel, mvi] = nanmax(v); [max_in, mii] = nanmax(mean_exc_input(:,ii));
    plot(synInputTime, v./max_vel, 'Color', colors{ii}, 'LineWidth', 1, 'Linestyle', '-'); hold on;
    ph = plot(synInputTime, mean_exc_input(:,ii)./max(mean_exc_input(:,ii)), 'Color', color_dark{ii}, 'LineWidth', 1, 'LineStyle', '--'); hold on
    % We also want to find the optimal delay between the velocity and excitatory synaptic input
    [optDelay(ii), optCorr] = findOptimalDelay(v(si:mii)./abs_max_vel, mean_exc_input(si:mii,ii)./abs_max_in, [], mean(diff(synInputTime))); 
end
% We now use the mean value of delay to shift all of the traces.
mean_optDelay = mean(optDelay);
for ii=1:length(reg_stim)
    plot(synInputTime-mean_optDelay, v./max_vel, 'Color', colors{ii}, 'LineWidth', 1, 'Linestyle', '-');
end
xlim([-400 100]);
xlabel('Time from Collision (ms)');
ylabel('Norm Syn Input / Velocity');

%% Here we are fitting the relationship between velocity and the excitatory synaptic input to the LGMD (Figure 7B)
transfig = figure;
vg_fitp = [];
options = optimset('TolFun', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'final', 'algorithm', 'trust-region-reflective');
for ii=1:length(reg_stim)
    dt = mean(diff(synInputTime)); di = round(optDelay(ii)/dt);
    di = round(mean_optDelay/dt);
   
    % make a velocity vector that matches the timebase of the synaptic time
    v = syn_input_v(:,ii);
    subplot(2,2,1); hold on;
    [max_vel, mvi] = nanmax(v); 
    [max_in, mii] = nanmax(mean_exc_input(:,ii));
    len = find(v >= min(max_vel,1000), 1,'first'); %fit only the speeds before the stimulus stops, or 1000 deg/s
    % Plotting input (shifted back in time) versus velocity
    plot(v(1:len),  mean_exc_input(1-di:len-di,ii), 'Color', colors{ii}, 'LineWidth', 1, 'Linestyle', '-'); hold on;
    vg_fitp(ii,:) = lsqcurvefit(powfun, [1 1], v(1:len),  mean_exc_input(1-di:len-di,ii), [-1e2 -1e2], [1e6 1e4], options);
    plot(v(1:len),  powfun(vg_fitp(ii,:), v(1:len)),  'Color', color_light{ii}, 'LineWidth', 1, 'Linestyle', '--'); hold on;
    vmax(ii) =  min(max_vel,1000)
    v_input_r2 = calcR2(mean_exc_input(1-di:len-di,ii),  powfun(vg_fitp(ii,:), v(1:len)))
    subplot(2,2,3);
    plot(synInputTime(1:len), v(1:len), 'Color', colors{ii}, 'LineWidth', 1, 'Linestyle', '-'); hold on;
    
end
vg_coeff = vg_fitp(:,2).* (vmax(:).^ vg_fitp(:,1)) %again, re-parameterize the fitted scale value as before
xlabel('Time (ms)');
ylabel('Angular speed (deg/s)');
subplot(2,2,1);
xlabel('Angular velocity (deg/s)');
ylabel('Total gsyn (mS/cm2)');
xlim([0 1000]);

%% plotting the transform between the inputs and proximal (near SIZ) Vm - Fig 7D
% Fitting options
logfun = @(p,x) rectify(p(1).* log(p(2).*x));
options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'final', 'algorithm', 'trust-region-reflective');
% Now resample the Exc synaptic vector to match the Vm
mean_exc_input_rs = resampleSignal(synInputTime, mean_exc_input, noinh_stim(1).tvec);
figure(transfig); subplot(2,2,2);
for ii=1:length(reg_stim)
    [~, maxi] = max(mean_exc_input_rs(:,ii)); % fit till the maximum of the trace
    ph = semilogx(mean_exc_input_rs(1:maxi,ii),  noinh_stim(ii).mu_delta_vmprox_filt(1:maxi), 'Color', colors{ii}, 'LineWidth', 1, 'LineStyle', '-'); hold on
    mei = mean_exc_input_rs(1:maxi,ii);
    mei_thresh = nanmax(mei) / 5000; %.2 percent of max 
    vmd =  noinh_stim(ii).mu_delta_vmprox_filt(1:maxi);
    fiti = mei >= mei_thresh;
    % The relationship is not evenly sampled, so in order to make the fitting procedure weight the curve evenly,
    % resample it.
    [eiX, Ymean, Ysd] = resampleVectorEvenly(mei(fiti), vmd(fiti), 'linear', 50);
    excg_vm_fitp(ii,:) = lsqcurvefit(logfun, [100 1], eiX, Ymean, [0 1e-2], [1e8 1e8], options);
    % plot the fit
    semilogx(mei(fiti), logfun(excg_vm_fitp(ii,:), mei(fiti)),'Color', color_light{ii}, 'LineWidth', 1, 'LineStyle', '--');
    r2 = calcR2(Ymean, logfun(excg_vm_fitp(ii,:), eiX))
end
semilogx(mei(fiti), logfun(mean(excg_vm_fitp), mei(fiti)),'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
xlim(16*[1e-3 20]);
ylim([0 45]);
ylabel('delta Vm (mV)'); xlabel('Total gsyn (mS/cm2)');
%set(gca, 'Xscale', 'linear'); xlim(15*[0 20]); %comment this line for semilogx plots, the inset
%set(gca, 'Xscale', 'log'); xlim(15*[1e-3 20]); 
setPresentationDefaults(gcf, 0);

%% Inhibition! - the script produces, amoungst other plots, Figures 7E,F,G.
compareModelsInhibition;

%%  Plot the spike transform of the model - Figure 7H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This calculates for each l/v value individually
% figure; axes; hold on;
% for ii = 1:length(reg_stim)
%     [datah, fith] = plotModelSpikeTransform(gca, reg_stim(ii).tvec, reg_stim(ii).a_vmprox_filt, reg_stim(ii).a_conv_inst_freq);
%     set(datah, 'Color', colors{ii}, 'LineWidth', 1);
%     set(fith, 'Color', colors{ii}, 'LineWidth', 1);
% end
% xlabel('Median Filtered Vm (mV)', 'FontSize', 16);
% ylabel('Firing Rate (Hz)', 'FontSize', 16);
% ylim([0 200]); 
% %xlim([-75 -20]);
% xlim([0 30]);
% setPresentationDefaults(gcf, 0);

%  Now, we are lumping the conditions (l/v values) to make a single spike threshold function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; axes; hold on;
all_vm_prox = [reg_stim.a_vmprox_filt];
all_inst_freq = [reg_stim.a_conv_inst_freq];
[datah, fith, spike_fitp, r2, spk_fun, spkthresh_data] = ...
                        plotModelSpikeTransform(gca, reg_stim(1).tvec, all_vm_prox, all_inst_freq);
% This converts the fitted coefficient to the one used in the paper
spkCoeff = spkthresh_data.vm(end)^spike_fitp(1) * spike_fitp(2) 
set(datah, 'Color', 'k', 'LineWidth', 1);
set(fith, 'Color', 'r', 'LineWidth', 1);
xlabel('Median Filtered Vm (mV)', 'FontSize', 16);
ylabel('Firing Rate (Hz)', 'FontSize', 16);
ylim([0 100]); xlim([0 20]);
setPresentationDefaults(gcf, 0);
