% Making a chart showing the relative contribution of different factors to the overall variability/SNR
if exist('compareLoomingVariability.mat', 'file')  %check for a saved data file...delete to regenerate
    load('compareLoomingVariability.mat');
end

% These define the files in which the full sets of simulations are saved.  They are omitted from the 
% DB submitted model due to space constraints, therefore if the below 'if ~exist..." were true, then this 
% script will error.  
fnroot = 'results';
filenames = {'averagedFull', 'averagedExcGconst', ...
    'averagedExcTconst',  'averagedNoexcvar', ...
    'averagedNoinhvar', 'averagedSpontonly'}; 

snr_tags = {'F - Exc Gmax', 'F - Exc Jitter', 'F - Exc Var', 'F - Inh Var'};
sd_tags = {'F - Exc Gmax', 'F - Exc Jitter', 'F - Exc Var', 'F - Inh Var', 'Spontaneous Only'};
full_tags = {'Full', 'F - Exc Gmax', 'F - Exc Jitter', 'F- Exc Var', 'F - Inh Var', 'Spontaneous Only'};
ext = '_analyzed.mat';

if ~exist('rsnr_fmax','var')
    for i = 1:length(filenames)
        load([fnroot '/' filenames{i} ext]);
        snr_fmax(i,:) = [stim.snr_fmax];
        snr_vmpeak(i,:) = [stim.snr_vmpeak];
        mu_fmax(i,:) = [stim.mu_fmax];
        mu_vmpeak(i,:) = [stim.mu_vmpeak];
        mu_fmax_t(i,:) = [stim.mu_fmax_t];
        mu_vmpeak_t(i,:) = [stim.mu_vmpeak_t];
        sd_fmax_t(i,:) = [stim.sd_fmax_t];
        sd_vmpeak_t(i,:) = [stim.sd_vmpeak_t];
        
        clear stim stim_str;
    end
end  


% Comparisons 
for i = 1:(length(filenames)-1)
    rsnr_fmax(i) = mean((snr_fmax(1,:)- snr_fmax(i,:))./(snr_fmax(1,:) - snr_fmax(end,:)));
    rsnr_vmpeak(i) =  mean((snr_vmpeak(1,:) - snr_vmpeak(i,:))./(snr_vmpeak(1,:) - snr_vmpeak(end,:)));
    rmu_fmax(i) = mean(mu_fmax(i,:)./mu_fmax(end,:));
    rmu_vmpeak(i) = mean(mu_vmpeak(i,:)./mu_vmpeak(end,:));
    rmu_fmax_t(i) = mean(mu_fmax_t(i,:)./mu_fmax_t(end,:));
    rmu_vmpeak_t(i) = mean(mu_vmpeak_t(i,:)./mu_vmpeak_t(end,:)); 
    
    % What we want is the fraction relative to full variability simulations
    cond= 1:3;
    rsd_fmax_t(i) = mean((sd_fmax_t(1,cond)-sd_fmax_t(i+1,cond))./sd_fmax_t(1,cond));
    rsd_vmpeak_t(i) = mean((sd_vmpeak_t(1,cond)-sd_vmpeak_t(i+1,cond))./sd_vmpeak_t(1,cond));
    rsd_fmax_t_lv(i,:) = (sd_fmax_t(1,:)-sd_fmax_t(i+1,:))./ sd_fmax_t(1,:);
end

% Do the plotting
figure;
subplot(2,2, 1);
xp =( 1:length(rsnr_fmax)) + .15;
bar(xp(2:end), rsnr_fmax(2:end), 'facecolor','k', 'barwidth', .25); hold on;
xp =( 1:length(rsnr_fmax)) - .15;
bar(xp(2:end), rsnr_vmpeak(2:end),'facecolor','r', 'barwidth', .25);
set(gca, 'xtick', 2:length(xp), 'xticklabel', snr_tags);
title('SNR_{peak}'); 
ylabel('SNR_{peak}, (Cond-Spont)/(Full - Spont)');
set(gca, 'tickdir', 'out');
%ylim([0 1.2]);
subplot(2,2,2);
xp =( 1:length(rsd_fmax_t)) + .15;
bar(xp, rsd_fmax_t, 'facecolor','k', 'barwidth', .25); hold on;
xp =( 1:length(rsd_vmpeak_t)) - .15;
bar(xp, rsd_vmpeak_t,'facecolor','r', 'barwidth', .25);
set(gca, 'xtick', 1:length(xp), 'xticklabel', sd_tags);
title('Peak Timing SD');
ylabel('SD_{peak time}, Fractional Difference'); 
set(gca, 'tickdir', 'out');
%set(gca, 'xticklabel', tags{2:end});

subplot(2,2,4);
lv = [10 40 80];
%plot(lv, rsd_fmax_t_lv);
plot(lv, sd_vmpeak_t);
legend(full_tags);
ylabel('IFR SD_{peak time} (ms)');
xlabel('l/|v| (ms)');

subplot(2,2,3);
plot(lv, snr_fmax);
ylabel('IFR SNR_{peak}');
