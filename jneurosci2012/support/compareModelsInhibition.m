
lin_thresh_fun = @(p, x) p(1) .* rectify(x - p(2));
theta_fun = @(loverv, t) 2 * atan(-loverv./t)*180/pi; %theta full-size
max_theta = 82; %this is the max because it is limited by the screen size we have in the lab

%% Subtraction section - effect of inhibition
% Here we want to do a subtraction between the two model versions in order to observe the effect of inhibition on the
% responses.  
options = optimset('TolFun', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'final', 'algorithm', 'trust-region-reflective'); %options for fitting
fh2 = figure;
stim_ah = axes('Parent', gcf, 'Position', [.1 .75 .8 .2], 'Visible', 'off'); hold on;
xlim([-800 200]);
vm_ah = axes('Parent', gcf, 'Position', [.1 .1 .8 .6]); hold on;
xlim([-800 200]);
mean_vm_diff = NaN*zeros(length(reg_stim(ii).tvec), length(reg_stim));
theta_inh = NaN* zeros(length(reg_stim(ii).tvec), length(reg_stim));
for ii=1:length(reg_stim)
    theta = theta_fun(reg_stim(ii).loverv, reg_stim(ii).tvec);
    theta(theta > max_theta) = max_theta;
    theta(theta<0) = max_theta;
    theta_inh(:,ii) = theta;
    plot(stim_ah, mparams(ii).mov_t, mparams(ii).theta, 'Color', colors{ii}); % plotting theta
    
    plot(thresh_time(ii)*ones(1,2), [0 1], 'Color', colors{ii});
    %plotting the traces
    mean_vm_diff(:,ii) = noinh_stim(ii).mu_vmprox_filt - reg_stim(ii).mu_vmprox_filt;
    off_dt(ii) = -mean_optDelay;
    %need to compare timecourses, so plotting the normalized versions on one another, Figure 7E
    norm_vm_diff = mean_vm_diff(:,ii) ./ max(max(mean_vm_diff));
    peak_vm_diff = max(mean_vm_diff(:,ii))
    norm_theta = mparams(ii).theta ./ max(mparams(ii).theta);
    ph1 = plot(vm_ah, reg_stim(ii).tvec, norm_vm_diff, 'Color', colors{ii}, 'LineWidth', 1);
    ph2 = plot(vm_ah, mparams(ii).mov_t + off_dt(ii), norm_theta, 'Color', colors{ii}, 'LineStyle', '--', 'LineWidth', 1);
end
ylabel('Inhibitory Influence on Vm');
xlabel('Time from Collision (ms)');
legend([ph1, ph2], {'Inhibitiory Influence', 'Stimulus \theta + toffset'});
ylim(vm_ah, [-.2 1]);
setPresentationDefaults(gcf, 0);

%% Plotting/fitting the angular size versus the inhibitory influence.
dt = mean(diff(reg_stim(1).tvec)); 
i_off = round(off_dt(ii)/dt); % the offset
theta_inh_linp = []; % initialize fit parameters to empty
% Since the model output and synaptic vectors are possibly on different timebases, we actually
% need to resample the synaptic input vector
mean_inh_input_rs = resampleSignal(synInputTime, mean_inh_input, reg_stim(ii).tvec);
% Fitting options
lb = [0 0 0]; ub = [100 500 1000];
options = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 4000, 'MaxIter', 1000, 'Display', 'notify', 'algorithm','trust-region-reflective');
figure;
abs_max_inh = max(max(mean_inh_input));
abs_max_vmdiff =  max(mean_vm_diff(:));
for ii=1:length(reg_stim)
   %make a theta vector on the same timebase
   theta = theta_fun(reg_stim(ii).loverv, reg_stim(ii).tvec);
   theta(theta > max_theta) = max_theta;
   theta(theta<0) = max_theta;
   theta_thresh =  2*atand(1/lgmd_peak_time_line(1));
   theta_inc = theta(theta<max_theta);
   timv = max(theta_inc)
   subplot(3,2,1);
   plot(reg_stim(ii).tvec + off_dt(ii), theta./max_theta, 'Color', color_light{ii}, 'LineWidth', 1); hold on;
   plot(reg_stim(ii).tvec, mean_vm_diff(:,ii)./ abs_max_vmdiff, 'Color', colors{ii}, 'LineWidth', 1, 'LineStyle', '--'); hold on;
   seli = (1:length(theta_inc));
   thetaX = theta_inc(seli);
   inhY = mean_vm_diff((1:length(seli))+i_off,ii);
   xlabel('Time (ms)');
   ylabel('Norm theta / Vm Inh Influence');
   %we are testing the segments over time
   subplot(3,2,3); 
   plot(reg_stim(ii).tvec(seli), thetaX, 'Color', colors{ii}, 'LineWidth', 1, 'LineStyle', '--'); hold on;
   plot(reg_stim(ii).tvec(seli), inhY, 'Color', color_dark{ii}, 'LineWidth', 1, 'LineStyle', '-'); hold on;
   [tX, Ymean, Ysd] = resampleVectorEvenly(thetaX, inhY, 'linear', max_theta+1);
    xlabel('Time (ms)');
   ylabel('Theta / Vm Inh Influence');
   %plotting the theta against the Inhibitory Influence (Vm) other, Figure 7F
   subplot(3,2,2); 
   plot(tX, Ymean, 'Color', colors{ii}, 'LineWidth', 1, 'LineStyle', '--'); hold on;
   plot(thetaX, inhY, 'Color', colors{ii}, 'LineWidth', 1); hold on;
   xlabel('Theta (deg)');
   ylabel('Vm Inh Influence');
   % Now fitting those relationships with Weibull functions
   theta_inhvm_fitp(ii,:) = lsqcurvefit(@weibull_func3, [1 50 2], tX(:), Ymean(:), lb, ub,options);
   theta_inh_linp(ii,:) =  theta_inhvm_fitp(ii,:);
   fitX = 0:max_theta;
   inh_weibull_r2 = calcR2(Ymean(:), weibull_func3(theta_inhvm_fitp(ii,:), tX(:)))
   plot(fitX, weibull_func3(theta_inhvm_fitp(ii,:), fitX), 'Color', color_light{ii}, 'LineStyle', '--', 'LineWidth', 1);
   ylim([0 40]); xlim([0 max_theta]);
   % Now we're going to look at the inhibitory Gsyn as a function of time
   subplot(3, 2, 4);
   plot(reg_stim(ii).tvec + off_dt(ii), theta./max_theta, 'Color', color_light{ii}, 'LineWidth', 1); hold on;
   plot(synInputTime, mean_inh_input(:,ii)./abs_max_inh, 'Color', color_dark{ii}, 'LineStyle', '--', 'LineWidth', 1); hold on;
   xlim([-800 200]); ylim([0 1]);
   xlabel('Time (ms)');
   ylabel('Norm theta / Total Ginh');
   % Now as a function of theta, Figure 7G
   subplot(3,2,5);
   sel_inh_input = mean_inh_input_rs((1:length(seli))+i_off,ii);
   plot(thetaX, sel_inh_input, 'Color', colors{ii}, 'LineWidth', 1); hold on;
   theta_inhg_fitp(ii,:) = lsqcurvefit(powfun, [1 1], thetaX, sel_inh_input, [0 0], [1e2 1e3], options);
   plot(thetaX, powfun(theta_inhg_fitp(ii,:), thetaX), 'Color', color_light{ii}, 'LineStyle', '--', 'LineWidth', 1);
   theta_inhsyn_r2(ii) = calcR2(sel_inh_input, powfun(theta_inhg_fitp(ii,:), thetaX))
   xlim([0 max_theta]); ylim([0 Inf]);
   xlabel('Theta (deg)');
   ylabel('Total Ginh');
end
% This is again reparameterizing the power function fit as before.
thetag_coeff =  theta_inhg_fitp(:,2).* (max_theta.^  theta_inhg_fitp(:,1))

setPresentationDefaults(gcf, 0);


