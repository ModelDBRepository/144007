% Analyze_model.m
% 
% Originally, Speron Oct. 7 2008
% Modified/extended by pwjones, 2010-2012.
%
function analyze_model(datadir, file_base)
  % --- preliminaries
  if ~exist('datadir')
	datadir = ['simout' filesep];
  else
    datadir = [datadir filesep];
  end
  
  % use a global variable, so need to clear that 
  clear global totalTrialN;
	  
  if (~exist('file_base'))
      % these are filenames for variability simulations looking at looming stimuli and response variability
      tags = {'loom_excTconst_snr_5_jit_6_exc_0.54_inh_0.0075'};
      tags = {'loom_excGconst_snr_0_jit_6_exc_0.54_inh_0.0075'};
      tags = {'loom_snr_10_jit_6_exc_0.54_inh_0.0075'};
  else
      tags = file_base;
  end

  % main loop
  for t=1:length(tags)

            % --- analyze the data ...
        analyze_tag(datadir, tags{t});
        
           % ---plot some of the data
        if strncmp(tags{t}, 'singlefacet', length('singlefacet'))
            plotSingleFacetData(datadir, tags{t});
        elseif strncmp(tags{t}, 'loom',4)
            if (t==1) fh = []; end;
            fh = plotLoomData(datadir, tags{t}, [1 1 0 1 1 1], fh);
        elseif strcmp(tags{t}, 'iclamp')
            plotIClampData(datadir, tags{t});
        end
  end
  
  % need to compare the variability for the different SNR values - breaking the flow of the code by putting things here, but this is quickest
  %compareLoomResps(datadir, tags);
  
% ----------------------------------------------------
function analyze_tag(datadir, tag)
% ---------------------------------------------------
% This will analyze a set of data - specified as an entry in TAGS, above


    outfile = [datadir tag '_analyzed.mat'];
	if (exist(outfile) ~= 0) % skip
	  disp([outfile ' already exists.  No analysis.']);
    else
      stim_str = tag;
      % These are the strings from stimuli that I'm using
      stim_str = {[tag '_gsnr_2_jit_6'], [tag '_gsnr_4_jit_6'], [tag '_gsnr_6_jit_6'], ... %suffixes for single facet stimuli
        [tag '_gsnr_8_jit_6'], [tag '_gsnr_10_jit_6'], [tag '_gsnr_15_jit_6'], ...
        [tag '_lv_10'], [tag '_lv_40'], [tag '_lv_80']}; %suffixes for looming
      
      inj_vals = [-3:0 2:10,12,15,20];% set the current injection values used, in nA
      for j = 1:length(inj_vals)
        stim_str{1+length(stim_str)} = sprintf('%s_na_%d_', tag, inj_vals(j));
      end
      % -- single facet
      for s=1:6
        flist = dir([datadir filesep stim_str{s} '*']); for f=1:length(flist) ; fnames{f} = flist(f).name ; end
          if (exist('fnames') & ~isempty(fnames))
              tstim = analyze_single_dataset(datadir, fnames, [100 125]);
              if (isstruct(tstim)); stim(s) = tstim ; end %#ok<AGROW>
          end
          fnames = {};
      end
      % -- Looming stimuli
      for s=7:9
          flist = dir([datadir filesep stim_str{s} '*']); for f=1:length(flist) ; fnames{f} = flist(f).name ; end
          if (exist('fnames') & ~isempty(fnames))
              tstim = analyze_single_dataset(datadir, fnames, [1900 2000]);
              if (isstruct(tstim)); stim(s) = tstim ; end %#ok<AGROW>
          end
          fnames = {};
      end
      % -- Current Injections
      for s=10:(9+length(inj_vals))
          flist = dir([datadir filesep stim_str{s} '*']); for f=1:length(flist) ; fnames{f} = flist(f).name ; end
          if (exist('fnames') & ~isempty(fnames))
              tstim = analyze_single_dataset(datadir, fnames, [200 350]);
              if (isstruct(tstim)); stim(s) = tstim ; end %#ok<AGROW>
          end
          fnames = {};
      end
      % --- save ...
		save(outfile, 'stim', 'stim_str', '-v7.3');
    end

		
% --------------------------------------------------------------
function stim = analyze_single_dataset(datadir, outfiles, t_fss)
% ---------------------------------------------------------------
% Analyzes a single type of data - all trials here are treated as the same stimulus
% 
  
  % flag for the type of firing rate filtering - gaussian smoothing (0) or ISI based (1).
  % Looming stimuli data are processed with gaussian smoothed rates, so use 0 for those.
  % For current injections, use the ISI based measures for comparability to those in Peron & Gabbiani (Nature Neurosci 2009)
  doISIfilt = 0; 
  if (strncmp(outfiles, 'iclamp', 6)) % if these are current injection trials, set the analysis method to be ISI based
          doISIfilt = 1; 
  end
  
  if (length(outfiles) == 0) ; stim = -1 ; return ; end

    % loop thru files
  nfiles = length(outfiles);
  for o=1:length(outfiles)
      disp(['Processing ' datadir '/' outfiles{o}]);

      % --- load data
      d = load([datadir '/' outfiles{o}]);
      
      % --- process
      spike_idx = get_spikes(-30, d(:,4));
      
      stim.tvec = d(:,1);
      if (strncmp(outfiles, 'loom', 4)) % if these are looming trials, adjust the time to be relative to collision
          stim.tvec = stim.tvec - (1900-50);
      end
      spike_vec = zeros(size(stim.tvec)); spike_vec(spike_idx) = 1; %spike vector for convolution
      
      sponti = find(d(:,1) > 20 & d(:,1) < 50);
      
      if (o==1) %set up guassian filters
        dt = mean(diff(stim.tvec));
        stdg = ceil(20/dt); % 20 ms in number of samples
        %stdg = ceil(.5/dt); % short gaussian useful for getting spike threshold of single spikes
        filtx = -3*stdg:1:3*stdg; 
        % this is the gaussian to convolve the spike train with.
        filty = my_normpdf(filtx,0,stdg);  filty = filty/sum(filty);
        % Simon's gaussian - if dt = .1ms, this is a 1 ms SD filter.
        gauss = dt*normpdf(-100:dt:100,0,2*dt); 
      end
      
      if (doISIfilt) % the ISI based method of computing the spike rates
        % gaussian convolve ... and normalize so that integral equals number of spks
        inst_freq = 1000*get_inst_freq(d(:,1), spike_idx); %uses spike timing for instantaneous freq
        gauss_if = conv(inst_freq, gauss);
        gauss_if = gauss_if((length(gauss)-1)/2:length(inst_freq) + (length(gauss)-1)/2 -1);
        denom = sum(gauss_if);
        if (denom == 0) %there are no spikes, and dividing by it would give NaNs
            denom=1;
        end
        gauss_if = gauss_if * (sum(inst_freq)/denom);
        gauss_if2 = gauss_if;
      else % the binary spike vector/gaussian filter method
        gauss_if2 = conv2(spike_vec(:),filty(:), 'same');
        %normalize to the number of spikes in the trial as in Gabbiani et al. 1999
        if (length(spike_idx) ~= 0) %if no spikes, then divide by 0 gives NaNs for data where 0s should be.
              gauss_if2 = gauss_if2 ./(nansum2(gauss_if2)*dt/1000);
              gauss_if2 = gauss_if2 * length(spike_idx);
        end
      end
      
      if (1 == 0) %optional plotting showing the difference between the two spike rate analysis methods
        figure; hold on;
        plot(stim.tvec, gauss_if, 'k', 'linewidth', 2);
        plot(stim.tvec, gauss_if2, 'b', 'linewidth', 2);
      end
      vm_axon = d(:,2);
      vm_prox = d(:,3);
      vm_siz = d(:,4);
      vm_dist = d(:,5);
      ca_siz = d(:,6);
      % median filter the dendritic Vm to remove spikes
      vm_filt = median_filt(d(:,5), 8, mean(diff(stim.tvec))); %8 msec width median filter
      vm_filt_prox = median_filt(d(:,3), 8, mean(diff(stim.tvec))); 
      vm_filt_siz = median_filt(vm_siz, 8, dt); 
      vm_rest(o) = mean(vm_filt(sponti));
      nsamp = length(vm_prox);
      
      % --- store relevant info
      stim.trial(o).conv_inst_freq = gauss_if2(:);
      stim.trial(o).spike_times = stim.tvec(spike_idx);
      stim.trial(o).spike_vec = spike_vec; 
      stim.trial(o).vm_filt = vm_filt;
      stim.trial(o).vm_filt_prox = vm_filt_prox;
      stim.trial(o).vm_filt_siz = vm_filt_siz;
      stim.trial(o).vm_axon =  d(:,2);
      stim.trial(o).vm_prox =  d(:,3);
      stim.trial(o).vm_siz =  d(:,4);
      stim.trial(o).vm_dist = d(:,5);
      stim.trial(o).vm_prox_std = std(vm_prox(sponti));
      stim.trial(o).vm_rest = vm_rest(o);
      stim.trial(o).ca_siz = ca_siz;
      stim.trial(o).gkca = d(:,7);
      stim.trial(o).ina = d(:,8);
      stim.trial(o).vm_dend_vent = median_filt(d(:,10), 2, dt);
      stim.trial(o).vm_dend_horiz = median_filt(d(:,11), 2, dt);
      stim.trial(o).vm_dend_dors = median_filt(d(:,12), 2, dt);
      
      nspks(o) = length(spike_idx);
      if (~isempty(t_fss))
          itemp = find(d(:,1) > t_fss(1), 1, 'first');
          if ~isempty(itemp)
              i_fss(1) = itemp; else
              i_fss(1) = 1;
          end
          itemp = find(d(:,1) >= t_fss(2), 1, 'first');
          if ~isempty(itemp) 
              i_fss(2) = itemp; else
              i_fss(2) = length(gauss_if2);
          end;
      else
          i_fss(1) = 1;
          i_fss(2) = length(gauss_if2);
      end
      [fmax(o), maxi] = max(gauss_if2(1:i_fss(1)));
      fss(o) = mean(gauss_if2(i_fss(1):i_fss(2)));
      vmss(o) = mean(vm_filt(i_fss(1):i_fss(2)));
      
      stim.trial(o).fmax = fmax(o);
      stim.trial(o).fmax_t = stim.tvec(maxi);
      stim.trial(o).fss = fss(o);
      stim.trial(o).vmss = vmss(o);
      stim.trial(o).nspks = nspks(o);
      [stim.trial(o).vmpeak, peaki] = nanmax2(vm_filt_prox - mean(vm_prox(sponti)));
      stim.trial(o).vmpeak_t = stim.tvec(peaki);
      %compute the widths of the responses.  Useful for some sim conditions.
      stim.trial(o).ifr_width = computePeakWidth(stim.tvec, gauss_if2, stim.trial(o).fmax_t, stim.trial(o).fmax, 1, mean(gauss_if2(1:100))); 
      stim.trial(o).vmpeak_width = computePeakWidth(stim.tvec, vm_filt, stim.trial(o).vmpeak_t, stim.trial(o).vmpeak, 1, mean(vm_filt(1:100))); 
      stim.trial(o).ca_ss = mean(ca_siz(i_fss(1):i_fss(2)));
      
      if (o==1) % allocate matrices
          a_conv_inst_freq = zeros(nsamp, nfiles); a_vm_filt = zeros(nsamp, nfiles);
          a_vm_filt_prox = zeros(nsamp, nfiles); a_vm_prox = zeros(nsamp, nfiles);
          a_vm_dist = zeros(nsamp, nfiles); a_vm_axon = zeros(nsamp, nfiles);
          a_spike_vec = zeros(nsamp, nfiles); a_ca_siz = zeros(nsamp, nfiles);
          a_m = zeros(nsamp, nfiles); a_h = zeros(nsamp, nfiles); 
          a_n = zeros(nsamp, nfiles);
          a_vm_dend_vent = zeros(nsamp, nfiles);
          a_vm_dend_horiz= zeros(nsamp, nfiles);
          a_vm_dend_dors= zeros(nsamp, nfiles);
          a_vm_siz=  zeros(nsamp, nfiles);
          a_vm_filt_siz = zeros(nsamp, nfiles);
      end   
      a_conv_inst_freq(:,o) = gauss_if2;
      a_vm_filt(:,o) = vm_filt;
      a_vm_filt_prox(:,o) = vm_filt_prox;
      a_vm_prox(:,o) =  stim.trial(o).vm_prox;
      a_vm_dist(:,o) =  stim.trial(o).vm_dist;
      a_vm_axon(:,o) = vm_axon;
      a_vm_siz(:,o) = vm_siz;
      a_vm_filt_siz(:,o) = vm_filt_siz;
      a_spike_vec(:,o) = stim.trial(o).spike_vec;
      a_ca_siz(:,o) = ca_siz;
      a_vm_filt_siz(:,o) = vm_filt_siz;
      a_spike_vec(:,o) = stim.trial(o).spike_vec;
      a_vm_dend_vent(:,o) = stim.trial(o).vm_dend_vent;
      a_vm_dend_horiz(:,o) = stim.trial(o).vm_dend_horiz;
      a_vm_dend_dors(:,o) = stim.trial(o).vm_dend_dors;
  end
  
  % get 'average' data
  nt = length(outfiles);
  stim.num_trials = nt;
  stim.mu_nspks = mean(nspks);
  stim.sd_nspks = std(nspks);
  stim.se_nspks = std(nspks)/sqrt(nt);
  stim.mu_fss = mean(fss);
  stim.sd_fss = std(fss);
  stim.se_fss = std(fss)/sqrt(nt);
  stim.mu_fmax = mean(fmax);
  stim.sd_fmax = std(fmax);
  stim.se_fmax = std(fmax)/sqrt(nt);
  stim.snr_fmax = stim.mu_fmax/stim.sd_fmax;
  stim.sd_fmax_t = std([stim.trial.fmax_t]);
  stim.mu_fmax_t = mean([stim.trial.fmax_t]);
  stim.mu_ifr_width = mean([stim.trial.ifr_width]);
  stim.sd_ifr_width = std([stim.trial.ifr_width]);
  
  stim.mu_conv_inst_freq = nanmean2(a_conv_inst_freq,2);
  stim.sd_conv_inst_freq = nanstd2(a_conv_inst_freq');
  stim.se_conv_inst_freq = nanstd2(a_conv_inst_freq')./sqrt(nt);
  stim.mu_vmrest = mean(vm_rest);
  stim.mu_vmfilt = mean(a_vm_filt,2);
  stim.sd_vmfilt = std(a_vm_filt');
  stim.mu_vmss = mean(vmss);
  stim.sd_vmss = std(vmss);
  stim.se_vmss = std(vmss)/sqrt(nt);
  stim.mu_vmpeak = mean([stim.trial.vmpeak]);
  stim.sd_vmpeak = std([stim.trial.vmpeak]);
  stim.snr_vmpeak = stim.mu_vmpeak / stim.sd_vmpeak;
  stim.sd_vmpeak_t = std([stim.trial.vmpeak_t]);
  stim.mu_vmpeak_t = mean([stim.trial.vmpeak_t]);
  stim.mu_vmpeak_width = mean([stim.trial.vmpeak_width]);
  stim.sd_vmpeak_width = std([stim.trial.vmpeak_width]);
  stim.a_vm_filt = a_vm_filt;
  stim.mu_vm_prox_std = mean([stim.trial.vm_prox_std]);
  stim.a_ca_siz = a_ca_siz;
  stim.a_vm_siz = a_vm_siz;
  stim.a_vm_filt_siz = a_vm_filt_siz;
  stim.mu_vmfilt_siz = mean(a_vm_filt_siz,2);
  stim.sd_vmfilt_siz = std(a_vm_filt_siz');
  stim.mu_ca_ss = mean([stim.trial.ca_ss]);
  stim.mu_vm_dend_vent = mean(a_vm_dend_vent, 2);
  stim.mu_vm_dend_horiz = mean(a_vm_dend_horiz, 2);
  stim.mu_vm_dend_dors = mean(a_vm_dend_dors, 2);
  
  %spikes not in axon are at a slight delay, compensate, then measure spike height
  tr = [50 200];
  prox_off = round(.29/dt);

    %plotting to visualize the model responses per trial type 
  if ( 1 == 1) 
      simple = 1; %how many subplots you wanna see 
      if simple 
          nr = 2; nc = 3;%rows, columns    
      else
          nr = 4; nc = 3;
      end
      plotn = 5; % number of individual traces to plot
      xl = [min(stim.tvec) max(stim.tvec)]; %limits
      
      figure; subplot(nr,3,1); hold on;
      title(outfiles{o});
      plotSpikeRasters(gca, stim.tvec, a_spike_vec', 0);
      ylim([0 20]); xlim(xl);
      
      %look at how the Vm itself looks
      subplot(nr,3,2);
      if (~simple)
          plot(stim.tvec, a_vm_siz(:,1:plotn), 'Color', 'k'); hold on;
          plot(stim.tvec, a_vm_axon(:,1:plotn), 'Color', 'c');
      end
      plot(stim.tvec, a_vm_prox(:,1:plotn), 'Color', 'r'); hold on;
      plot(stim.tvec, a_vm_dist(:,1:plotn), 'Color', 'b');
      ylabel('Vm (mV)');  xlim(xl);
      
      %look at the filtered Vm and the firing rate
      subplot(nr,3,3);
      plot(stim.tvec, a_vm_filt(:,1:plotn), 'Color', [.5 .5 .5]); hold on;
      plot(stim.tvec, stim.mu_vmfilt, 'linewidth', 2);
      xlabel('time (ms)'); ylabel('Vm (mV)');  xlim(xl);
      subplot(nr,3,5); plot(stim.tvec, stim.sd_vmfilt);
      xlabel('time (ms)'); ylabel('Vm \sigma (mV)'); xlim(xl);
      subplot(nr,3,4); % these will be for the firing rates
      plot(stim.tvec, a_conv_inst_freq(:, 1:plotn), 'Color', [.5 .5 .5]); hold on;
      plot(stim.tvec, stim.mu_conv_inst_freq, 'linewidth', 2);
      xlabel('time (ms)'); ylabel('IFR (Hz)'); xlim(xl);
      subplot(nr,3,6); plot(stim.tvec, stim.sd_conv_inst_freq);
      xlabel('time (ms)'); ylabel('IFR \sigma (Hz)'); xlim(xl);
      if (~simple)
          subplot(nr,3,7); plot(stim.tvec, stim.a_ca_siz(:,1:plotn), 'Color', 'g');
          xlabel('time (ms)'); ylabel('[Ca]'); xlim(xl);
          subplot(nr,3,8); plot(stim.tvec, stim.mu_gkca, 'k-'); ylabel('I_{KCa}'); xlim(xl);
          subplot(nr,3,9); plot(stim.tvec, stim.mu_ina, 'k-'); ylabel('I_{Na}'); xlim(xl);
          subplot(nr,3,10); 
          plot(stim.tvec, stim.a_m(:,1:plotn), 'b', stim.tvec, stim.a_h(:,1:plotn), 'r', stim.tvec, stim.a_n(:,1:plotn), 'g')  
          %graph the [Ca] versus firing rate.  This should be approximately linear, but better fit by a square root
          %relationship, as in Wang (1998)
          stimt = [50 200]; stimi = stim.tvec >= stimt(1) & stim.tvec <= stimt(2);
          subplot(nr,3,11);
          meanCa = mean(stim.a_ca_siz(stimi,:), 2) .* 1e3;
          plot(meanCa, stim.mu_conv_inst_freq(stimi), 'k.'); hold on;
          
      end
  end

  
  
% ------------------------------------------
function plotSingleFacetData(datadir, tag)
% ------------------------------------------
% This function plots the response characteristics of single facet responses as a function of their 
% input parameters.
  anfile = [datadir tag '_analyzed.mat'];
  load(anfile);
  bootn = 5000;
  
  for i=1:length(stim_str)
     val = textscan(stim_str{i}, 'singlefacet_gain_%d_gsnr_%d_jit_%d');
     try
         gain(i) = val{1}; gsnr(i) = val{2}; t_jitter(i) = val{3};
     end
  end
  
  % Calculate and print the "spontaneous" variability of the more distal membrane potential
  sponti = find(stim(1).tvec > 50 & stim(1).tvec < 100);
  spont_vm = zeros(length(sponti), length(stim(1).trial), length(stim));
   snrfun = @(x)mean(x)./std(x);
  for i = 1:length(stim)
        temp = [stim(i).trial.vm_dist];
        spont_vm(:,:,i) = temp(sponti,:);
        ci = bootci(bootn, snrfun, [stim(i).trial.vmpeak]);
        snr_ci(i) = stim(i).snr_vmpeak - ci(1); %take the lower diff from the mean - assumes symmetry
        ci = bootci(bootn, @std, [stim(i).trial.vmpeak_t]);
        sd_vmpeak_t_ci(i) = stim(i).sd_vmpeak_t - ci(1);
  end
  disp(sprintf('Spontaneous membrane potential: mean = %d  ,  sd = %d', mean(spont_vm(:)), std(spont_vm(:))));
 
  figure; 
  % plot the single facet response height
  subplot(2, 2, 1);
  plot(gsnr, [stim.mu_vmpeak], 'k-o');
  addErrorBars(gca, gsnr, [stim.mu_vmpeak], [stim.sd_vmpeak], 'k');
  xlabel('G_{syn} SNR'); ylabel('Vm_{peak}');
  % plot the SNR of the peak response
  subplot(2,2,2);
  plot(gsnr, [stim.snr_vmpeak], 'k-o');
  addErrorBars(gca, gsnr, [stim.snr_vmpeak], snr_ci, 'k');
  xlabel('G_{syn} SNR'); ylabel('Vm_{peak} SNR');
  subplot(2,2,3);
  plot(gsnr, [stim.sd_vmpeak_t], 'k-o');
  addErrorBars(gca, gsnr, [stim.sd_vmpeak_t], sd_vmpeak_t_ci, 'k');
  xlabel('G_{syn} SNR'); ylabel('Vm_{peak} time \sigma (ms)');
  subplot(2,2,4);
  plot(gsnr, [stim.mu_vmpeak_width], 'k-o');
  addErrorBars(gca, gsnr, [stim.mu_vmpeak_width],[stim.sd_vmpeak_width], 'k');
  xlabel('G_{syn} SNR'); ylabel('Vm_{peak} width (ms)');
  
% -----------------------------------
function plotIClampData(datadir, tag)
%------------------------------------- 
% Analyzes the current injection data

% This function plots the response characteristics of single facet responses as a function of their 
% input parameters.
  anfile = [datadir tag '_analyzed.mat'];
  load(anfile);
  for i=1:length(stim) injstims(i) = ~isempty(stim(i).tvec); end %create logical vect for where data exists.
  stim = stim(injstims); stim_str = stim_str(injstims); %just select the subset that are current injections
  for i=1:length(stim_str)
     val = textscan(stim_str{i}, 'iclamp_na_%f');
     try
         na(i) = val{1}; %get the current injection values
     end
  end

% plot the Vm responses to current injections and the input resistance
figure;
subplot(2, 2, 1);
plot(na, [stim.mu_vmpeak]+[stim.mu_vmrest], 'k-o'); hold on;
plot(na, [stim.mu_vmss], 'r-o');
xlabel('Current Injection (nA)'); ylabel('Vm');
legend(gca, {'Peak', 'SS'});
addErrorBars(gca, na, [stim.mu_vmpeak]+[stim.mu_vmrest], [stim.sd_vmpeak], 'k',.25);
addErrorBars(gca, na, [stim.mu_vmss], [stim.sd_vmss], 'r',.25);
vmss = [stim.mu_vmss];
negi = na < 0;
if sum(negi)>=2
    fitp = polyfit(na(negi), vmss(negi), 1)
    text(min(na), stim(1).mu_vmrest-3, ['R_{in} =' num2str(fitp(1)) ' M\Omega']);
    disp(['Input resistance is ' num2str(fitp(1)) ' MOhms']);
    plot(na,polyval(fitp, na), 'r--');
end

% plot the firing rates in response to the current injections
subplot(2,2,2); hold on;
plot(na, [stim.mu_fss], 'k-o');
plot(na, [stim.mu_fmax], 'r-x');
xlabel('Current Injection (nA)'); ylabel('SS IFR');
legend(gca, {'Steady State', 'Peak'});
subplot(2,2,3);
plot(na, [stim.sd_vmpeak_t], 'k-o');
xlabel('Current Injection (nA)'); ylabel('Vm_{peak} time \sigma (ms)');
subplot(2,2,4);
plot(stim(1).tvec, [stim.mu_conv_inst_freq]');

%plot the steady state [Ca] versus firing rate - this should be roughly linear
figure; axes; hold on;
plot([stim.mu_ca_ss], [stim.mu_fss], 'ko-');
xlabel('[Ca]'); ylabel('F_{ss} (Hz)');

 % Calculate and print the "spontaneous" variability of the more distal membrane potential
 sponti = find(stim(1).tvec > 30 & stim(1).tvec < 100);
 spont_vm = zeros(length(sponti), length(stim(1).trial), length(stim));
 for i = 1:length(stim)
     temp = [stim(i).trial.vm_dist];
     spont_vm(:,:,i) = temp(sponti,:);
 end
 disp(sprintf('Spontaneous membrane potential: mean = %d  ,  sd = %d', mean(spont_vm(:)), std(spont_vm(:))));
 dvm=.1;
 vm_range = [-inf, -70:dvm:-58, inf];
 bin_centers = [vm_range(2)-dvm*3/2, (vm_range(2:end-1) + dvm/2), vm_range(end-1)+dvm*3/2];
 counts = histc(spont_vm(:), vm_range);
 figure; bar(bin_centers, counts);

% ------------------------------------------------
function fh = plotLoomData(datadir, tag, pb, fh)
% -----------------------------------------------
% This function plots the response characteristics of single facet responses as a function of their 
% input parameters. pb is a vector containing boolean values about which plots to make.

plot_ind = 0;
anfile = [datadir tag '_analyzed.mat'];
load(anfile);
for i=1:length(stim) loomstims(i) = ~isempty(stim(i).tvec); end %create logical vect for where data exists.
stim = stim(loomstims); 
stim_str = stim_str(loomstims); %just select the subset that are the looming
%extract the SNR value for the dataset from the file name
val = {9}; %assigning an snr of 9?
snr = double(val{1}); %and extracting it
color_scale = (10-snr)/10; %since we know it goes to 10
load([datadir '/loomingParameters.mat']);
dt = mean(diff(stim(1).tvec));
pc = {[0 1 0], [1 0 0], [0.3 0.6 1], [0 0 0]};
time_range = [-1500 150];

% Setting up the axes for plotting
if pb(1) % The Firing Rate plots
    if (length(fh)>=1 && fh(1) ~= 0)
        figure(fh(1));
        mu_ifr_ah = findobj(fh(1), 'tag', 'frah'); 
        mu_stim_ah = findobj(fh(1), 'tag', 'frstimah');
        mu_spk_ah = findobj(fh(1), 'tag', 'frspkah');
    else
        fh(1) = figure; mu_ifr_ah = axes('parent', fh(1), 'Position', [.1 .1 .8 .4], 'tag', 'frah'); hold on; xlim(time_range);
        mu_stim_ah = axes('parent', gcf, 'Position', [.1 .8 .8 .15], 'Visible', 'off', 'tag', 'frstimah'); hold on; xlim(time_range);
        mu_spk_ah = axes('parent', gcf, 'Position', [.1 .51 .8 .29], 'Visible', 'off', 'tag', 'frspkah'); hold on;  xlim(time_range);
    end
    ylabel(mu_ifr_ah, 'IFR (Hz)'); xlabel(mu_ifr_ah, 'Time (ms)');
end
if pb(2) % Vm plots 
    if (length(fh)>=2 && fh(2) ~= 0)
        figure(fh(2));
        mu_vm_ah = findobj(fh(2), 'tag', 'vmah');
        sd_vm_ah = findobj(fh(2), 'tag', 'vmsdah');
        vm_stim_ah = findobj(fh(2), 'tag', 'stimah');
    else
        fh(2) = figure;
        mu_vm_ah = axes('parent', fh(2), 'Position', [.1 .1 .8 .3], 'tag','vmah'); hold on; xlim(time_range);
        sd_vm_ah = axes('parent', fh(2), 'Position', [.1 .45 .8 .25], 'tag','vmsdah'); hold on; xlim(time_range);
        vm_stim_ah = axes('parent', fh(2), 'Position', [.1 .75 .8 .2], 'Visible', 'off', 'tag', 'stimah'); hold on; xlim(time_range);
    end
    ylabel(sd_vm_ah, 'Vm SD (mV)'); xlabel(mu_vm_ah, 'Time (ms)'); set(sd_vm_ah, 'TickDir', 'out');
end
if pb(6) %  %Figure for the IFR SD over time
    if (length(fh)>=6 && fh(6) ~= 0)
        figure(fh(6));
        ifr_sd_ah = findobj(fh(6), 'tag', 'ifr_sd_ah');
        ifrsd_stim_ah = findobj(fh(6), 'tag', 'ifrsd_stim_ah');
    else
        fh(6) = figure;
        ifr_sd_ah = axes('parent', fh(6), 'Position', [.1 .1 .8 .6], 'tag','ifr_sd_ah'); hold on; xlim(time_range);
        ifrsd_stim_ah = axes('parent', fh(6), 'Position', [.1 .75 .8 .2], 'tag','ifrsd_stim_ah'); hold on; xlim(time_range);
    end
    ylabel(ifr_sd_ah, 'IFR SD (Hz)');
end

global totalTrialN;
if isempty(totalTrialN) 
    ti = 0; 
else
    ti = totalTrialN; 
end

for i=1:length(stim)
    % find the l/v value from the string
    val = textscan(stim_str{i}, [tag '_lv_%d']);
    stim(i).loverv = double(val{1});
    
    if pb(1) % spiking figure
        % plot the mean line
        plot(mu_ifr_ah, stim(i).tvec, stim(i).mu_conv_inst_freq, 'Color', pc{i}, 'LineWidth', 2);
        errc = (pc{i} + [1 1 1])/2; %define the color for error envelope, halfway between the mean line and white
        % plot error envelope
        plot_err_poly(mu_ifr_ah, stim(i).tvec, stim(i).mu_conv_inst_freq, stim(i).sd_conv_inst_freq, pc{i}, errc, 1, 10, time_range);
        line(loomStimParameters(i).mov_t, abs(loomStimParameters(i).theta), 'Color', pc{i}, 'Linewidth', 2,'Parent',mu_stim_ah);
        % Now make and plot the rasters
        spk_m = zeros(length(stim(i).trial(1).spike_vec), length(stim(i).trial));
        for j=1:length(stim(i).trial) % make the spike matrix
            spk_m(:, j) = stim(i).trial(j).spike_vec;
        end
        rh = plotSpikeRasters(mu_spk_ah, stim(i).tvec, spk_m', ti);
        set(mu_spk_ah, 'xlim', time_range); set(rh, 'Color', pc{i});
        ti = ti + length(stim(i).trial);
    end
    if pb(6) %Figure for the IFR SD over time
        line(loomStimParameters(i).mov_t, abs(loomStimParameters(i).theta), 'Color', pc{i}, 'Linewidth', 2,'Parent', ifrsd_stim_ah);
        line(stim(i).tvec, stim(i).mu_conv_inst_freq(:)./stim(i).sd_conv_inst_freq(:), 'Parent',ifr_sd_ah,'Color', pc{i}, 'LineWidth', 1);
    end    
    if pb(2) % Vm figure
        figure(fh(2));
        % Now plot the filtered Vm for each condition
        plot(mu_vm_ah, stim(i).tvec, stim(i).mu_vmfilt, 'Color', pc{i}, 'LineWidth', 2);
        set(mu_vm_ah, 'TickDir', 'out');
        ylabel(mu_vm_ah, 'V_m (mV)');
        plot(vm_stim_ah, loomStimParameters(i).mov_t, abs(loomStimParameters(i).theta), 'Color', pc{i}, 'Linewidth', 2); hold on;
        xlabel('Time (ms)'); ylabel('Filtered Vm (mV)'); set(mu_ifr_ah, 'Tickdir', 'out');
        % Also, plot the SD of the filtered Vm in another panel.
        plot(sd_vm_ah, stim(i).tvec, stim(i).sd_vmfilt, 'Color', pc{i}, 'LineWidth', 2); 
    end
    
    if pb(3) % If we want to plot individual trials to check them out
        fh(3) = figure;
        ifr_ah = axes('parent', gcf, 'position', [.1 .1 .8 .6]); hold on; xlim(time_range);
        spk_ah = axes('parent', gcf, 'position', [.1 .72 .8 .2], 'visible', 'off'); xlim(time_range);
        for j=1:length(stim(i).trial)
            line('Parent', ifr_ah, 'xdata', stim(i).tvec, 'ydata', stim(i).trial(j).conv_inst_freq, 'color', [.5 .5 .5]);
        end
        plotSpikeRasters(spk_ah, stim(i).tvec, spk_m', 0);
        set(spk_ah, 'xlim', time_range); set(ifr_ah, 'xlim', time_range);
    else
        fh(3) = 0;
    end
    
    %let's print a few things out
    disp(sprintf('l/v = %d', stim(i).loverv));
    disp(sprintf('IFR: mean = %g   SD = %g  SNRpeak=%g SDpeak_time = %g', ...
        stim(i).mu_fmax, stim(i).sd_fmax, stim(i).snr_fmax, stim(i).sd_fmax_t));
    disp(sprintf('Vm: mean = %g   SD = %g  SNRpeak=%g SDpeak_time = %g', ...
        stim(i).mu_vmpeak, stim(i).sd_vmpeak, stim(i).snr_vmpeak, stim(i).sd_vmpeak_t));
    
end

totalTrialN = ti; %set a global variable for the total number of trials.
clear totalTrialN;

if pb(4) % plot the SNR values for the responses
    if (length(fh)>=4 && fh(4) ~= 0)
        figure(fh(4));
    else
        fh(4) = figure;
    end
    subplot(1,2,1); hold on;
    black = [0 0 0]; red = [1 0 0];
    
    snr_max = 21;
    bc = black + ([1 1 1] - black)*(snr_max-snr)/snr_max - [.1 .1 .1];
    rc = red + ([1 1 1] - red)*(snr_max-snr)/snr_max - [0 .1 .1];
    plot([stim.loverv], [stim.snr_fmax], 'ko-', 'MarkerSize', 10, 'LineWidth', 2, 'Color', bc);
    plot([stim.loverv], [stim.snr_vmpeak], 'r-o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', rc);
    legend_str = {'Peak IFR', 'Peak Vm'}; legend(legend_str);
    xlabel('l/v (ms)'); ylabel('SNR');
    set(gca, 'Xtick', [stim.loverv], 'Xticklabel', cellstr(num2str([stim.loverv]')), 'Tickdir', 'out');
    subplot(1,2,2); hold on;
    plot([stim.loverv], [stim.sd_fmax_t], '-ok', 'MarkerSize', 10, 'LineWidth', 2, 'Color', bc);
    plot([stim.loverv], [stim.sd_vmpeak_t], '-or', 'MarkerSize', 10, 'LineWidth', 2, 'Color', rc);
    xlabel('l/v (ms)'); ylabel('Peak Timing \sigma (ms)');
    set(gca, 'Xtick', [stim.loverv], 'Xticklabel', cellstr(num2str([stim.loverv]')), 'Tickdir', 'out');
    legend(legend_str);
    % compute bootstrap confidence intervals
    for i = 1:length(stim)
         fm = [stim(i).trial.fmax];
         vmp = [stim(i).trial.vmpeak];
         fm_t = [stim(i).trial.fmax_t];
         vmp_t = [stim(i).trial.vmpeak_t];
         fm_ci = bootci(5000, @snrfun, fm); %bootstrapped 95% confidence intervals for snr
         vmp_ci = bootci(5000, @snrfun, vmp);
         fm_t_ci = bootci(5000, @std, fm_t); 
         vmp_t_ci = bootci(5000, @std, vmp_t);
         subplot(1,2,1);
         addErrorBars(gca, stim(i).loverv, vmp_ci(1) + diff(vmp_ci)/2, diff(vmp_ci)/2, rc);
         addErrorBars(gca, stim(i).loverv, fm_ci(1) + diff(fm_ci)/2, diff(fm_ci)/2, bc);
         subplot(1,2,2);
         addErrorBars(gca, stim(i).loverv, vmp_t_ci(1) + diff(vmp_t_ci)/2, diff(vmp_t_ci)/2, rc);
         addErrorBars(gca, stim(i).loverv, fm_t_ci(1) + diff(fm_t_ci)/2, diff(fm_t_ci)/2, bc);
         %save these values back to the structure
         stim(i).snr_fmax_ci = fm_ci;
         stim(i).snr_vmpeak_ci = vmp_ci;
         stim(i).sd_fmax_t_ci = fm_t_ci;
         stim(i).sd_vmpeak_t_ci = vmp_t_ci;
    end
    
end

if pb(5)
    % Also, in order to check that responses to looming are like those that
    % we see in vivo, would like to check the timing of the peak, and whether
    % it is linear.
    fh(5) = figure;
    for i=1:2
        subplot(1,2,i);
        if i==1
            plot([stim.loverv], -([stim.mu_vmpeak_t]), 'o-b', 'linewidth', 2, 'MarkerSize', 10); hold on;
            addErrorBars(gca, [stim.loverv], -([stim.mu_vmpeak_t]), [stim.sd_vmpeak_t], 'b');
            title('Vm Peak');
            %plot the individual peak times
            x = zeros(max([stim.num_trials]), length(stim)); y = zeros(max([stim.num_trials]), length(stim));
            [fit_p, ~, ~, ci] = linear_fit_stats(x(:), y(:))
            plot([stim.loverv], polyval(fit_p, [stim.loverv]), 'b--', 'LineWidth', 1);
        else
            plot([stim.loverv], -([stim.mu_fmax_t]), 'o-b', 'linewidth', 2, 'MarkerSize', 10); hold on;
            addErrorBars(gca, [stim.loverv], -([stim.mu_fmax_t]), [stim.sd_fmax_t], 'b');
            x = zeros(max([stim.num_trials]), length(stim)); y = zeros(max([stim.num_trials]), length(stim));
            for j=1:length(stim)
                x(:,j) = stim(j).loverv*ones(stim(j).num_trials,1);
                y(:,j) = -[stim(j).trial.fmax_t]';
            end
            [fit_p, ~, ~, ci] = linear_fit_stats(x(:), y(:))
            plot([stim.loverv], polyval(fit_p, [stim.loverv]), 'b--', 'LineWidth', 1);
            fitstr=sprintf('fit: slope = %.2g, int= %.2f, thresh= %.2f', fit_p(1), fit_p(2), 2*atand(1/fit_p(1)));
            title(fitstr);
            
        end
        xlabel('l/v (ms)'); ylabel('Peak Time before Collision (ms)');
        hold on;
        
        % Note that since using polyval, the intercept here has the opposite sign as in the other functions.
        plot(10:10:80,  polyval([4.39 -18], 10:10:80, 1), 'r'); %My looming data
        %plot(10:10:80,  polyval([2.9 26], 10:10:80, 1), 'r'); %haleh DCMD peak (Fotowat and Gabbiani 2007)
        %plot(10:10:80,  polyval([4.7 -27], 10:10:80, 1), 'r'); %Fabrizio DCMD peak (Gabbiani et al 1999)
    end
end

% Let's save a summary that will be quickly loaded.
save(anfile, 'stim', 'stim_str', '-v7.3');


function snr = snrfun(vals)
meanval = mean(vals);
sdval = std(vals);
if (sdval == 0)
    snr = 0;
else
    snr = meanval./sdval;
end
	 
% --------------------------------------------
function spike_idx = get_spikes(thresh, data)
% --------------------------------------------
% Extracts spike idx, based on threshold + peak method
% 
greater = find(data > thresh);
diffs = diff(data);
greater = greater(find(greater < length(diffs)));

spike_idx = [];
greater = greater(find(greater < length(diffs)));

for g=1:length(greater)-1
    if (diffs(greater(g)) > 0 & diffs(greater(g)+1) < 0)
        spike_idx(length(spike_idx)+1) = greater(g)+1;
    end
end
	
 
  
% --------------------------------------------------------------
function [inst_freq_vals] = get_inst_freq(time_vals, spike_idx)
% --------------------------------------------------------------
% Computes instantaneous firing frequency and returns it along with SD
%

  % The heart of the matter
	inst_freq_vals=zeros(1,length(time_vals));
	spike_times = time_vals(spike_idx);
	if (length(spike_times) > 1)
		for i=1:length(inst_freq_vals)
			pre_spike = find(spike_times < time_vals(i));
			pre_spike = max(pre_spike);
			post_spike = find(spike_times > time_vals(i));
			post_spike = min(post_spike);

			% Are we AT a spike?
			if (length(find(spike_times == time_vals(i))) == 1)
				if (length(pre_spike) > 0 & length(post_spike) > 0)
					inst_freq_vals(i) = 0.5/(time_vals(i)-spike_times(pre_spike)) ...
					+ 0.5/(spike_times(post_spike)-time_vals(i)); 
				% At first spike?
				elseif (length(pre_spike) == 0)
					inst_freq_vals(i) = 1/(spike_times(post_spike)-time_vals(i)); 
				% At last spike?
				elseif (length(post_spike) == 0)
					inst_freq_vals(i) = 1/(time_vals(i)-spike_times(pre_spike)); 
				end
			% Not at a spike
			else
				% Before first spike or after last?
				if (length(pre_spike) == 0 | length(post_spike) == 0) 
					inst_freq_vals(i) = 0; 
				else
					inst_freq_vals(i) = 1/(spike_times(post_spike)-spike_times(pre_spike)); 
				end
			end
		end
	end       


