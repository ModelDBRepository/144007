function [stim_out, vel_tuning] = loadModelResults(stim, mean_offset_t)
% function [stim_out, vel_tuning] = loadModelResults(stim, mean_offset_t)
%
% just makes sure that the correct Vm measures are present in the model data and 
% then normalizes them, creating new "_norm" fields in the stim data structure.
spont_range = 200:500;
if(~isfield(stim, 'mu_delta_vmprox_filt'))
    for ii = 1:length(stim)    
        vm_prox_m = zeros(length(stim(ii).mu_vmfilt),length(stim(ii).trial));
        ifr_m = zeros(length(stim(ii).mu_vmfilt),length(stim(ii).trial));
        for jj = 1:length(stim(ii).trial)
            vm_prox_m(:,jj) = stim(ii).trial(jj).vm_filt_prox;
            ifr_m(:,jj) = stim(ii).trial(jj).conv_inst_freq;
        end
        stim(ii).mu_vmprox_filt = mean(vm_prox_m,2);
        stim(ii).a_vmprox_filt = vm_prox_m;
        stim(ii).mu_delta_vmprox_filt = stim(ii).mu_vmprox_filt - mean(stim(ii).mu_vmprox_filt(spont_range));
        stim(ii).a_conv_inst_freq = ifr_m;
    end
end

% This is the Vm value to use
vm_measure = 'mu_vmprox_filt';

% Now, let's plot the same types of plots for the LGMD model responses
types = {'vm', 'ifr'};
for j = 1:2
    type = types{j};
    for ii=1:length(stim)
         % Firing Rates
        switch j 
            case 1
                [vel, resp, vel_t] = computeLoomingVelocityTuning(stim(ii).tvec, stim(ii).(vm_measure), stim(ii).loverv, mean_offset_t.cc);
            case 2
                [vel, resp, vel_t] = computeLoomingVelocityTuning(stim(ii).tvec, stim(ii).mu_conv_inst_freq, stim(ii).loverv, mean_offset_t.ifr);
        end
        if(ii == 1)
            vel_m = NaN*zeros(length(vel), length(stim));
            resp_m = NaN*zeros(length(vel), length(stim));
            vel_t_m = NaN*zeros(length(vel), length(stim));
        end
        vel_m(1:length(vel),ii) = vel;
        resp_m(1:length(vel),ii) = resp;
        vel_t_m(1:length(vel),ii) = vel_t;  
    end
    %get normalized values
    delta_resp_m = resp_m - mean(resp_m(1,:));
    norm_resp_m = delta_resp_m ./ max(max(delta_resp_m));

    vel_tuning.(type).vel_m = vel_m;
    vel_tuning.(type).resp_m = resp_m;
    vel_tuning.(type).vel_t_m = vel_t_m;
    vel_tuning.(type).delta_resp_m =  delta_resp_m;
    vel_tuning.(type).norm_resp_m = norm_resp_m;
end

stim_out = stim;


function [vel, resp, t_out] = computeLoomingVelocityTuning(t_in, response, loverv, tshift)
% Just plots the instantaneous dependence of the response on the velocity
% tshift ms in the past.  This shift is to compensate for the latency.
% Time and response must be the same size.
pb=0;
dt = mean(diff(t_in));
ishift = round(tshift/dt);
front_t = (-ishift:-1)* dt + t_in(1); nfront = length(front_t);
end_t = (1:-ishift)*dt + t_in(end); nend = length(end_t); %for negative time shifts - to the future
vel_t = [front_t(:); t_in(:); end_t(:)]; % we want to calculate velocities for comparison before because of shifting
max_theta = 41; %maximum half-angle
theta = -atan(loverv./vel_t)*180/pi; %this is the half angle of the stimulus
stoppedi = find(theta >= max_theta, 1, 'first');
theta(stoppedi:end) = max_theta;
vel_func = @(loverv, t) -loverv./(t.^2+loverv^2)*180/pi .* 1000; %gives deg/s
shifted_vel =  abs(vel_func(loverv, vel_t));%deg/s
shifted_vel(stoppedi:end) = 0; %velocity is zero when stimulus is stopped
% need to repeat for the non-extended arrays
short_theta = atan(loverv./t_in)*180/pi; %this is the half angle of the stimulus
short_vel = abs(vel_func(loverv, t_in));

[~,maxi] = nanmax2(response);
colli = find(t_in <= 0, 1, 'last');
endi = max(maxi, colli);
pz = 1:endi;
t_out = t_in(pz);
resp = response(pz);
short_vel = short_vel(pz);

sel = (1:length(resp)) + nfront - ishift;
vel = shifted_vel(sel);
vel_t_orig = vel_t;
vel_t = vel_t(sel);

if (pb)
    figure;
    subplot(2,1,1); hold on;
    plot(t_out, short_vel, 'k');
    plot(t_out, vel, 'r');
    legend('Normal', 'Shifted');
    subplot(2,1,2); hold on;
    plot(vel, resp, '-r');
    plot(short_vel, resp, '-k');
    legend('Shifted', 'Not Shifted');
    xlabel('Velocity (deg/sec)');
    ylabel('Response');
end
