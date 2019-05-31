% modified by PWJ on 12/30/10-2012 from generate_dirsel_sims to generate the simulations for single facets and looming
% modified by FG on 12/22/10
%
% Function generate_variability_sims
%
%  This is the main generation function -- it generates a large number of simulations 
%  that employ the rake_final morphology.  Output consists of NEURON files that store
%  data at several compartments (see below).
%
function generate_variability_sims

  % --- Preliminaries
  template_path = 'active_template.t';
  hocfile_rootpath = 'hocfiles'; % where hocfiles will go and from where simulations are run
  outfile_path = 'simout'; % directory name where simulation output is saved
  % number of cpus (cores) to use for simulations vector notation 
  % spreads files over sum(n_cpus) with n_cpus(i) in separate batch 
  % files, runnable on separate machines
  %n_cpus = [6 7 10 5]; 
  n_cpus = [14];
  %where the 'special' program (the compiled NEURON binary) created by nrnivmodl is located for each machine 
  %special_dir = { '../x86_64', '../umac', '../umac',  '../umac'}; 
  special_dir = { '../i386'};
  dt = .01; % simulation timestep in ms
  %randn('state',sum(100*clock)); % sets the seed to be different every time matlab is run, randn underlies normrand.
  
  % --- Global simulation settings - paths, etc.

  % the params(xxx) structure is passed to generate_hocfiles ; <%params(x).idstr%> in the
  %  template file (template_path above) is replaced with params(x).value.

  % the morphology file
  params(1).idstr = 'morphoFile';
  params(1).value = 'rake_final';
  % where the hoc files, once run, output their RESULTS (in text files)
  params(2).idstr = 'outfilePath';
  params(2).value = ['../' outfile_path]; 

  % the format of the output text file; the seg[X].v(.5) write the Vm over time of a single compartment in the 
  % model. The segments listed correspond to [460]:axonal, [475]: proximal dendrite near SIZ, [468]: SIZ, [485]: distal end of dendritic trunk
  params(3).idstr = 'saveParamsLine';
  params(3).value = ['fprint(" %f %f %f %f %f %f %f %f %f %f %f %f\n", t, seg[460].v(.5), seg[475].v(.5), seg[468].v(.5), seg[485].v(.5),'...
      'seg[470].cai(.5), seg[470].cm_KCa(.5), seg[485].i_cap(.5), seg[485].i_pas(.5),seg[110].v(.5),seg[210].v(.5), seg[310].v(.5))'];

  % Copy the hoc file to the simulation directory in case it has been altered in the interim
  %copyfile([params(1).value '.hoc'], [hocfile_rootpath '/' params(1).value '.hoc']);  %use this in most conditions
  % above line throws an error due to using 'cp -p' on a linux nfs mount of a zfs storage pool
  unix(['cp ' params(1).value '.hoc ' hocfile_rootpath '/' params(1).value '.hoc']);

  % --- Neuron-wide settings
  rm = 10350; % membrane resistance, Ohms/cm^2
  rm_nospont = 4900; % this is the membrane resistance to give an input resistance of ~5Mohm at seg[485] w/o spont activity
  ri = 60; % intracellular resistivity Ohm x cm

  % --- Settings for the active conductances

  % HH channels for spiking
  gna =  0.217; %% sodium conductance unitary, S/cm^2, SIZ
  gkdr = 0.3532; % % potassium conductance, S/cm^2, SIZ
  gna_axon = gna/4;% axon conductances
  gkdr_axon =gkdr/4;
  gna_dend = 0; gkdr_dend = 0; %dendritic conductances are not currently used
  
  gna = [gna gna_axon gna_dend]; %for reducing parameters to pass below
  gkdr = [gkdr gkdr_axon gkdr_dend];
  
  % Ca parameters - adjusted to get exponential decay from peak firing rates during current injection and a ~150Hz firing rate at steady state 
  gcal = 0.005; % calcium channel; allows calcium entry, which contributes minimally to voltage
  gkca = 0.2; % the real effect of calcium is that it opens these K channels, calcium-sensitive potassium conductance (S/cm^2)
  %gcal = 0; gkca = 0; %Uncomment to disable Ca and Spike Frequency Adaptation
  
  % --- Spontaneous synapses
  % excitatory, mimic nicotinic acetylcholine receptors
  spont_exc.freq = .4/1000; % firing freq of EACH syn, in
                           % units of 1/ms = kHz, 1000 ms
                           % between events or 1 event per s
  spont_exc.gmax = .012; % max conductance of spontaneous
                          % synapses in micro Siemens, was .001
                         
  spont_exc.erev = 0; % reversal potential in mV
  spont_exc.tau_syn = 0.3; % tau parameter for alpha synapse --
                           % time to peak in ms
  spont_exc.nsyns = 5000; % how many events will you see?
                         % freq*nsyns/ms; 
                         %timing drawn from a normal distribution 
  spont_exc.cmpt = 1:400;% which compartments will they be
                         % located in?

  % inhibitory, mimic gaba_A receptors ; fields as above
  spont_inh.freq = 3.344/1000; %1.75/1000 works for single facet stuff
  spont_inh.gmax = .0054;
  spont_inh.erev = -75;
  spont_inh.tau_syn = 3;
  spont_inh.nsyns = 600;
  spont_inh.cmpt = 476:480;
  
  % --- Visual synapses 
  % excitatory, mimic nicotinic acetylcholine receptors
  vis_exc.gmax = .6; % peak conductance in microSiemens; 1.2 with nspf 8 was a good pair
  vis_exc.erev = 0; % reversal potential 
  vis_exc.tau_syn = 0.3; % rise time constant
  vis_exc.nspf = 6; % number of synapses per "facet"
  vis_exc.gstd = vis_exc.gmax/5; %variability on the gmax, std of a normal distribution (gmax being the mean)
  vis_exc.tjitter = 6; %ms, the timing jitter of visual excitatory synapses.  Times drawn from a normal distribution with tjitter STD
  vis_exc.cmpt = 1:400; % compartment numbers for the synapses
    
  % inhibitory
  vis_inh.gmax = .0095; %was .1 - at low values I was using .0025
  vis_inh.erev = -75;
  vis_inh.tau_syn = 3;
  vis_inh.nspf = 4; % number of synapses per "facet"
  vis_inh.gstd = vis_inh.gmax/20; %variability on the gmax, std of a normal distribution (gmax being the mean)
  vis_inh.tjitter = 10; %ms, the timing jitter of visual inhibitory synapses.  Times drawn from a normal distribution with tjitter STD
  vis_inh.delay = 70; %ms, latency delay relative to the luminance change at a facet 
  vis_inh.cmpt = 476:480; % compartment numbers for the synapses
  
  % to generate no-inhibition simulations, these are substituted (no conductance)
  vis_noinh.gmax = 0;
  vis_noinh.erev = 0;
  vis_noinh.tau_syn = 3;
  vis_noinh.nspf = 1; % number of synapses per "facet"
  vis_noinh.gstd = 0;
  vis_noinh.tjitter = 15; %ms
  vis_noinh.delay = 50; %ms 
  vis_noinh.cmpt = 476:480;
  
  % to generate simulations without spontaneous activity, use these parameters
  no_spont_exc.freq = 1/1000; 
  no_spont_exc.gmax = 0.002; 
  no_spont_exc.erev = 0;
  no_spont_exc.tau_syn = 0.3; 
  no_spont_exc.nsyns = 1; 
  no_spont_exc.cmpt = 1:400;
  % inhibitory, mimic gaba_A receptors ; fields as above
  no_spont_inh.freq = 1/1000; 
  no_spont_inh.gmax = 0.002; 
  no_spont_inh.erev = -70;
  no_spont_inh.tau_syn = 3;
  no_spont_inh.nsyns = 1;
  no_spont_inh.cmpt = 476:480;

  % --- for scaling synapses for facet density
  map_sf = [1 5]; % min/max scaling factor ; thus, at highest
                  % density, conductance is upped 5x
                  % it is at least 1x

  % --- the heart of the matter -- call generate_vis_stims
  % repetitively to generate a particular set of parameter
  % combinations each call will generate many hocfiles.  Note
  % that since hoc_idx is incremented internally, it is
  % perfectly safe to comment out as many of the calls as you
  % want to exclude a particular parameter set.  See
  % "generate_vis_stims" definition below for more details of
  % what these mean. 

  hoc_idx = 1; % this is incremented at every call so that we keep hocfile numbering consistent
 
  % 1) Single facet simulation 
  retmap = 1; % use retinotopic map to place synapses
  gain = 1;
  for j=1:length(gain)  %looops through different gain values. 
        sf_spont_exc = spont_exc; %sf_spont_exc.gmax = spont_exc.gmax*gain(j); 
        sf_spont_inh = spont_inh; %sf_spont_inh.gmax = spont_inh.gmax*gain(j);
        sf_vis_exc = vis_exc; 
        sf_vis_exc.gmax = vis_exc.gmax*gain(j); sf_vis_exc.gstd = vis_exc.gstd*gain(j);
%        hoc_idx = generate_vis_stims (params, rm, ri, map_sf, [gkdr gkdr_axon], [gna gna_axon], gkca, gcal, ...
%              sf_spont_exc, sf_spont_inh, sf_vis_exc, vis_noinh, retmap, ['singlefacet_gain_' num2str(gain(j))], ...
%              hoc_idx, template_path, hocfile_rootpath, outfile_path,  dt, 1, 10);
  end
  
  % 2) Looming simulation [each block: 3 l/v values * 1 snr value= 3 sims] 
  retmap = 1; % use retinotopic map to place synapses
       
  snr = 5; %[2 5 10 15] if you want to explore using different snr values for gmax as in jones & gabbiani (2012, Jneurophysiol fig 7b,c)
  for k = 1:length(snr) % loop through snr values
        if (snr(k) ~= 0) %set the excitatory SNR
            vis_exc.gstd = vis_exc.gmax./snr(k);
        else
            vis_exc.gstd = 0;
        end
       % get the SNR: useful mostly if we aren't already setting it here.
       if (vis_exc.gstd ~= 0) exc_snr = vis_exc.gmax/vis_exc.gstd; else exc_snr = 0; end 
       hoc_idx = generate_vis_stims(params, rm, ri, map_sf, gkdr, gna, gkca, gcal, ...
         spont_exc, spont_inh, vis_exc, vis_inh, retmap,  sprintf('loom_snr_%g_jit_%g_exc_%g_inh_%g', ...
         exc_snr, vis_exc.tjitter, vis_exc.gmax, vis_inh.gmax), hoc_idx, template_path, hocfile_rootpath, outfile_path, dt, 2, 10);
  end     
 
  
  % 3) Current Injections (set up to do multiple currents, with the current values given within generate_vis_stims)
  retmap = 1; % use retinotopic map to place synapses
 
  spont_gain = 1;
  spont_inh.gmax = spont_inh.gmax * spont_gain;
  spont_exc.gmax = spont_exc.gmax * spont_gain;
  
  % Current injections with ongoing synaptic input - 'noisy' ones 
%    hoc_idx = generate_vis_stims (params, rm, ri, map_sf, [gkdr gkdr_axon], [gna gna_axon], gkca, gcal, ...
%              spont_exc, spont_inh, vis_exc, vis_inh, retmap, 'iclamp', ...
%              hoc_idx, template_path, hocfile_rootpath, outfile_path,  dt, 3, 10);
  
  % Current injections without ongoing synaptic input ('noise') - be sure to set the Rm value to the alternative 
  % if running these
  % rm = rm_nospont; % this is the adjusted membrane resistance, as above
  % hoc_idx = generate_vis_stims (params, rm, ri, map_sf, [gkdr gkdr_axon], [gna gna_axon], gkca, gcal, ...
  %           no_spont_exc, no_spont_inh, vis_exc, vis_inh, retmap, 'iclamp', ...
  %           hoc_idx, template_path, hocfile_rootpath, outfile_path,  dt, 3, 1);
                
  %NOTE: you may generate additional conditions by using other
  %combinations of the flags controlling facet density,
  %retinoptopic placement of synapses and pre-synaptic direction
  %selectivity
  
  % generate sh files. This makes hocfile running MUCH easier -- we have 8 CPUs but leave 1 open
  % We call with hoc_idx-1 because the variable keeps the index of the next hoc file. 
  generate_sh_files(n_cpus, 1, hoc_idx-1, hocfile_rootpath, special_dir);

% ----------------------------------
% Function generate_vis_stims
% ---------------------------------
%  Generates simulation of visual stimulation, specifically by setting up unique
%  patterns of synaptic activity.  This will make repeated calls to
%  a synaptic block generator followed by a hocfile generator. Parameters:
%
%  params - things to assign in the hocfile that are global
%  rm, ri - passive properties (membrane resistance, intracellular resistivity)
%  gkdr, gna, gkca, gcal - peak conductance for various channels (HH K, HH Na, calcium-dep K, calcium)
%  spont_exc, spont_inh - the spontaneous synapses' properties
%  vis_exc, vis_inh - the visual synapses' properties
%  use_retmap - set to 1 and will assign exc vis syns topographically; otherwise, 
%  random position so effectively no electrotonic structure
%  tag - root tag for out simulation files
%  hoc_idx - idx for hoc files to start at; returns last one
%  template_path - template file to use
%  hocfile_rootpath - where hocfiles will end up
%
%  This function works by building up the params structure, then generating 
%  hocfiles using it. pi is incremented as it is the index in the params 
%  structure (it has nothing to do with circles or the number pi!)
%
function hoc_idx = generate_vis_stims (params, rm, ri, map_sf, gkdr, gna, gkca, gcal, ...
      spont_exc, spont_inh, vis_exc, vis_inh, use_retmap, tag, hoc_idx, ...
			template_path, hocfile_rootpath, outfile_path,  dt, mode, N)

  testmode = 0; % set to 1 to test; will then generate just a few simulations

  % VERY IMPORTANT: this vector determines which simulations are
  % run. A value of 1 means the corresponding simulations are
  % generated; 0 means no. Convenient if you want a subset or are
  % debugging. ORDER: single-facet, looming, current injection
  enabled = [0 0 0];
  enabled(mode) = 1;
 

  % --------- pre-preliminaries --------------------------------
  if (isempty(N))
    N = 100; % The number of iterations per simulation type that we will do 
  end
  if (testmode == 1)
        disp(['WARNING : running generator in test mode; I will only generate a FEW ' tag]);
        N = ceil(N/10);
  end

  % dt is determined by the argument given.  The model readings are saved at a lower rate.
  % The quality of read firing rates degrades with samping dts > .03 ms.
  samp_dt = dt*3;
  if (samp_dt > .03) samp_dt = .03; end
  
  % --------- preliminaries ------------------------------------
  % first, assign the passed parameters to the params
  % structure. Each parameter is converted to a string that will be
  % inserted at the appropriate location in the template hoc file
  % to generate the final hoc file 

  % The next 3 parameters will have to be assigned INDIVIDUALLY in
  % each section below, These parameters are the tag of the
  % outfile, the synaptic stimulation block and the simulation length
  pi = length(params) + 1;
  params(pi).idstr = 'outfileTag'; otpi = pi; pi = pi+1; % the output file that simulation writes to; note index; increment
  params(pi).idstr = 'synBlock'; sbpi = pi; pi = pi+1; % synaptic block
  params(pi).idstr = 'simLen'; slpi = pi ; pi = pi+1; % duration of simulation

  params(pi).idstr = 'dt';
  params(pi).value = num2str(dt); pi = pi+1;
  params(pi).idstr = 'samp_dt';
  params(pi).value = num2str(samp_dt); pi = pi+1;
  
  % passive props - global
  params(pi).idstr = 'Rm';
  params(pi).value = num2str(rm); pi = pi+1;
  params(pi).idstr = 'Ri';
  params(pi).value = num2str(ri); pi = pi+1;

  % conductanges - global
  params(pi).idstr = 'gKDr';
  params(pi).value = num2str(gkdr(1)); pi = pi+1;
  params(pi).idstr = 'gKDrAxon';
  params(pi).value = num2str(gkdr(2)); pi = pi+1;
  params(pi).idstr = 'gKDrDend';
  params(pi).value = num2str(gkdr(3)); pi = pi+1;
  params(pi).idstr = 'gNa';
  params(pi).value = num2str(gna(1)); pi = pi+1;
  params(pi).idstr = 'gNaAxon';
  params(pi).value = num2str(gna(2)); pi = pi+1;
  params(pi).idstr = 'gNaDend';
  params(pi).value = num2str(gna(3)); pi = pi+1;
  params(pi).idstr = 'gKCa';
  params(pi).value = num2str(gkca); pi = pi+1;
  params(pi).idstr = 'gCaL';
  params(pi).value = num2str(gcal); pi = pi+1;

  for n=1:N %number of repetitions
    
    % ----------- Single facet stimuli -----------------------------
    if enabled(1)
        simlen = 250;
        params(slpi).value = num2str(simlen);
        
        gstd_vals = [vis_exc.gmax/2 vis_exc.gmax/4 vis_exc.gmax/6 vis_exc.gmax/8 vis_exc.gmax/10 vis_exc.gmax/15];
        jitter_vals = [4 6 7 8 10 12 16];
        j = 2;
        for k=1:length(gstd_vals)
            
            % synaptic block - spontaneous
            syn_idx = 1;
            [synblock_spont syn_idx] = get_spontaneous_synapses(spont_exc, spont_inh, simlen, syn_idx);
            
            %mode parameter: 4 = single facet stimuli
            mode = 1;
            
            % inhibitory and excitatory visual synaptic blocks
            % here, we pass the stimulus parameters gvs_params, which
            % will be used to determine synapse positioning, see
            % function get_rake_box_cmpts immediately below
            gvs_params = [];
            vis_exc_sf = vis_exc; %just so that I don't mess with the default one
            vis_exc_sf.gstd = gstd_vals(k); vis_exc_sf.tjitter = jitter_vals(j);
            [synblock_vis syn_idx] = get_visual_synapses(vis_exc_sf, vis_inh, simlen, syn_idx, mode, ...
                gvs_params, use_retmap, map_sf);
            
            %update synaptic block value
            params(sbpi).value = [synblock_spont synblock_vis]; % SYNBLOCK HERE!
            
            %subtag added to the tag passed below
            subtag = sprintf('_gsnr_%d_jit_%d', vis_exc_sf.gmax/vis_exc_sf.gstd, vis_exc_sf.tjitter);
            params(otpi).value = [tag subtag];
            
            % The actual generation ...
            hoc_idx = generate_hocfiles(template_path, hocfile_rootpath, outfile_path,  params, hoc_idx, otpi, n);
        end
                                 
    end %if enabled(1)
    % --------------- Looming stimuli --------------------------------
    if enabled(2)
        % The generateLoomingSynanpticInput function uses the global LOOMSTIMPARAMS to create trials with identical
        % synaptic timing. Clearing this variable causes it to regenerate the synaptic timings, and we do this every
        % 50 trials to have sets of trials with identical timings. If we aren't trying to do identical timings, clearing
        % the variable doesn't hurt anything
        if(~mod(n-1, 50))
            clear global LOOMSTIMPARAMS; %this clears the LOOMSTIMPARAMS variable holding facet timing information - forces a new random vector
        end
        if (1 == 1) %if we want to plot the synaptic strengths over time.
            sfh = findobj('Tag', 'synapse_fig');
            if (~isempty(sfh)) 
                sah(1) = findobj('Tag', 'synapse_ax1'); %we are finding persistant axes
                sah(2) = findobj('Tag', 'synapse_ax2');
            else
                sfh = figure('Tag', 'synapse_fig');
                sah(1) = axes('Parent', sfh, 'Position', [.1 .55 .8 .4], 'Tag', 'synapse_ax1');
                sah(2) = axes('Parent', sfh, 'Position', [.1 .05 .8 .4], 'Tag', 'synapse_ax2');
            end
        else
            sah = [];
        end
            
        simlen = 2000; % This is not the full loom for some of the longer ones, but certainly most of it.
        params(slpi).value = num2str(simlen);
        loverv = [10 40 80]; %looming parameter
        loom_l = 400; %looming object half-size
        loom_v = loom_l./loverv;  %approach velocity
        for i=1:length(loom_v)
            % synaptic block - spontaneous
            syn_idx = 1;
            [synblock_spont syn_idx] = get_spontaneous_synapses(spont_exc, spont_inh, simlen, syn_idx);
            
            %mode parameter: 2 = looming
            mode = 2;
            gvs_params = [loom_v(i) loom_l sah];
            [synblock_vis syn_idx] = get_visual_synapses(vis_exc, vis_inh, simlen, syn_idx, mode, ...
                gvs_params, use_retmap, map_sf);
            
            %update synaptic block value
            params(sbpi).value = [synblock_spont synblock_vis]; % SYNBLOCK HERE!
            
            %subtag added to the tag passed below
            subtag = sprintf('_lv_%d', loom_l/loom_v(i));
            params(otpi).value = [tag subtag];
            
            % The actual generation ...
            hoc_idx = generate_hocfiles(template_path, hocfile_rootpath, outfile_path,  params, hoc_idx, otpi, n);
        end
    end
    % -------------------- Current injections ------------------------------
    if enabled(3) % current injections 
        % these will still have the normal spontaneous inputs, which are variable trial to trial, but
        % there is nothing random to the 'stimulus' in these.  Just straight current injections.
        simlen = 500;
        params(slpi).value = num2str(simlen);
        mode = 3;
         % synaptic block - spontaneous
        syn_idx = 1;
        [synblock_spont syn_idx] = get_spontaneous_synapses(spont_exc, spont_inh, simlen, syn_idx);
        %injected current in nA, small negative to get input resistance, and larger positive for the i/o curves
        amps = [-3:0 2:10 12 15 20]; 
        inj_cmpt = 485.5; %485.5 is base of dend - where we would be recording in vivo
        for i=1:length(amps) %generate trials for each current injection value
            clampstr = generate_iclampstr(inj_cmpt, 100, 350, amps(i)); %location, onset time, duration, amplitude
            
            %update synaptic block value
            params(sbpi).value = [synblock_spont clampstr]; % SYNBLOCK HERE!
            
            subtag = sprintf('_na_%d', amps(i));
            params(otpi).value = [tag subtag];
            
            % The actual generation ...
            hoc_idx = generate_hocfiles(template_path, hocfile_rootpath, outfile_path,  params, hoc_idx, otpi, n);
        end    
    end
  end %for n=1:N
	

% ---------------------------------------------------------------
function [synblock syn_idx] = get_visual_synapses(exc, inh, simlen, syn_idx, mode, visparams, use_retmap, map_sf)
% -------------------------------------------------------
%  returns a synblock for visual activity.  This is perhaps the
%  most essential function. Parameters:
%
%  exc - excitatory visual synapse parameters
%  inh - inhibitory visual synapse parameters
%  simlen - length of simulation in ms
%  syn_idx - index of current synaptic block; returned with
%  synaptic block
%  mode: 1 = single facet stimuli
%        2 = looming stimuli
%  visparams: gives extra stimulus info 
%  use_retmap: 0 - should excitatory synapses be mapped retinotopically across the dendrites
%  or just be distributed RANDOMLY in tree 
% ------------------------------------------------------
  
  % --- initialize synaptic variables
  inh_cmpt = inh.cmpt;
  exc_cmpt = exc.cmpt;

  %initialize the synblock and compartment, timing synaptic vectors
  synblock = '';
  inh_cmpt_vec = [];
  exc_cmpt_vec = [];
  inh_timing_vec = [];
  exc_timing_vec = [];

  pre_time = 50; %in ms
  post_time = 50; %in ms
  
  %correct for pre- and post- time.
  simlen = simlen-(pre_time + post_time);

  if (mode == 1) 
      % single facet stimulation: a small group of synapses are placed in a single compartment and
      % activated nearly simultaneously (with some jitter).
      cmpt_base = 210; %along the middle of the AP axis, along the middle of the length of the dendrite.
      nspikes = 3;
      
      exc_cmpt_vec = cmpt_base + rand(1,exc.nspf); % all in the same compartment
      exc_gsyn_vec = rectify(normrnd(exc.gmax, exc.gstd, 1, exc.nspf))/nspikes; %jitter in gsyn
      exc_timing_vec = normrnd(100, exc.tjitter, 1, exc.nspf); %timing jitter
      % So, the single synapse responses are too brief to comprise single facet responses.  
      % Thus, we implement trains of presynaptic stimulation for each input, as postulated in Jones and Gabbiani (2010).
      for i = 1:(nspikes-1) %short burst of n spikes
        exc_cmpt_vec = cat(1, exc_cmpt_vec, exc_cmpt_vec(1,:));
        exc_gsyn_vec = cat(1, exc_gsyn_vec, exc_gsyn_vec(1,:));
        exc_timing_vec = cat(1, exc_timing_vec, exc_timing_vec(1,:)+5*i); % 200 Hz 
      end
      exc_cmpt_vec = exc_cmpt_vec(:)'; exc_gsyn_vec = exc_gsyn_vec(:)'; exc_timing_vec = exc_timing_vec(:)';
      
      inh_cmpt_vec = inh_cmpt(1);
      inh_timing_vec = 0;
      inh_gsyn_vec = inh.gmax;
     
  elseif (mode==2) 
        % This is the stimulation pattern for looming. All arguments are just copied from the function args.
     [exc_gsyn_vec,exc_cmpt_vec,exc_timing_vec,inh_gsyn_vec,inh_cmpt_vec,inh_timing_vec] = ...
        generateLoomingSynapticInput(exc, inh, simlen, syn_idx, mode, visparams,  ...
                                    use_retmap, map_sf);
  end % ending stimulus mode specific synaptic vector generation
  
  
  % shuffle synapses randomly across tree?
  if (use_retmap == 0)
    exc_cmpt_vec = ceil(max(exc_cmpt)*rand(1,length(exc_cmpt_vec)));
  end
  
  
  % plot distribution of synaptic timings and weights - set pb=1 to see this. It is sometimes nice.
  pb=0;
  if (pb)
    % Histograms of input counts  
    figure; subplot(2,2,1);
    hist(inh_timing_vec,100);
    h = findobj(gca,'Type','patch'); set(h,'FaceColor','r','EdgeColor','k'); %red for stop
    subplot(2,2,3);
    hist(exc_timing_vec,100);
    length(exc_timing_vec);
    
    % now plot the strengths of the individual inputs
    subplot(2,2,2);
    plot(inh_timing_vec, inh_gsyn_vec, 'r.');
    subplot(2,2,4);
    plot(exc_timing_vec, exc_gsyn_vec, 'b.'); hold on;
    plot(inh_timing_vec, inh_gsyn_vec, 'r.'); %give a scale to the inhibitory input.
    
    %Now show the summed gsyn over time during the stimulus
    bin_len = 5; edges = 0:bin_len:simlen; 
    % make time binned vectors that are summed synaptic weights
    exc_weighted = zeros(1, length(edges)-1); inh_weighted = zeros(1, length(edges)-1);
    for j=1:length(exc_timing_vec)
        ii = find(edges >= exc_timing_vec(j) , 1, 'first'); %find the time bin it's in
        if (~isempty(ii) && ii > 1 && ii <= length(exc_weighted)) exc_weighted(ii-1) = exc_weighted(ii-1) + exc_gsyn_vec(j); end %add the weight
    end
    for j = 1:length(inh_timing_vec)
         ii = find(edges >= inh_timing_vec(j) , 1, 'first'); %find the time bin it's in
         if (~isempty(ii) && ii > 1 && ii <= length(inh_weighted)) inh_weighted(ii-1) = inh_weighted(ii-1) + inh_gsyn_vec(j); end %add the weight
    end
    figure; subplot(2,1,1);
    weight_t = edges(1:end-1) + bin_len/2;
    plot( weight_t, exc_weighted, 'b', weight_t, inh_weighted, 'r');
    subplot(2,1,2);
    plot(weight_t, exc_weighted./inh_weighted);
  end

  % generate synblocks ...
  synblock = [synblock generate_synblock(exc_timing_vec+pre_time, exc_cmpt_vec, exc_gsyn_vec, exc.tau_syn, exc.erev, 0, syn_idx, [])];
  syn_idx = syn_idx + length(exc_timing_vec);
  if (exist('inh_gsyn_vec')) %sometimes there is just a fixed inhibitory gsyn, sometimes it's variable
      synblock = [synblock generate_synblock(inh_timing_vec+pre_time, inh_cmpt_vec, inh_gsyn_vec, inh.tau_syn, inh.erev, 0, syn_idx, [])];
  else % inh.gmax is constant
      synblock = [synblock generate_synblock(inh_timing_vec+pre_time, inh_cmpt_vec, inh.gmax, inh.tau_syn, inh.erev, 0, syn_idx, [])];
  end
  syn_idx = syn_idx + length(inh_timing_vec);

  
% -----------------------------------------------------------------------------------------------
%                       generateLoomingSynapticInput
% -----------------------------------------------------------------------------------------------
function [exc_gsyn_vec,exc_cmpt_vec,exc_timing_vec,inh_gsyn_vec,inh_cmpt_vec,inh_timing_vec] = ...
            generateLoomingSynapticInput(exc, inh, simlen, syn_idx, mode, visparams, ...
                                         use_retmap, map_sf)
% ------------------------------------------------------------------------------------------------
% 
% The general idea for this set of synaptic generation is to use a facet front-end that calculates the luminance
% transition caused by each stimulus.  It then outputs a set of luminance change durations and times that can be used to
% trigger the synaptic inputs based on the observed timing and magnitudes of single facet responses. It has a mountain of
% arguments and return values, but is just meant to be called in one place, within 'get_visual_synapses' of 
% 'generate_variability_sims'.
%
% It takes in all of the iputs to 'get_visual_synapses' and outputs all of the vectors necessary for making synapses
% onto the model - excitatory and inhibitory strengths, compartments, and timings.

% Flags used to recreate data in fig
% These are flags to save and repeat the timings for each stimulus.  In order to separate the variabilities due to synaptic strength
% versus timing, the idea is to randomly jitter the timing once then use that same set for each trial. The timings can be saved with 
% loomStimParams. Generally, this flag is zero, and new random timings are generated each trial/run.
repeatExcTimings = 0;
repeatInhTimings = 0;
% Another flag for running simulations with excitatory synaptic jitter constant, rather than varying with stimulus speed
constExcJitter = 0;

loverv = visparams(2)/visparams(1); % l and v
if (length(visparams)>2)
    syn_ah = visparams(3);
    syn_ah2 = visparams(4);
    pb=1;
else 
    pb=0;
end


% We are going to use a global variable for the structure that holds the stimulus' resulting input times since it takes
% a while to generate and this function is called for every stimulus presentation in the model.
global LOOMSTIMPARAMS;
if (isempty(LOOMSTIMPARAMS))
    if (exist('looming_stim_params.mat', 'file')) %check if there is a saved file we can load
        load('looming_stim_params.mat');
    else
        LOOMSTIMPARAMS = photArrayLoomOutput(loverv); %start a new structure
        stimParams = LOOMSTIMPARAMS;
    end
end
si = find([LOOMSTIMPARAMS.loverv] == loverv,1,'first');
if (isempty(si)) %the l/v value is not in the struct, run the stimulus generation
    stimParams = photArrayLoomOutput(loverv);
    if (isfield(LOOMSTIMPARAMS, 'exc_timing_vec')) %need to keep the struct array fields consistent
        stimParams.exc_timing_vec = []; stimParams.inh_timing_vec = [];
    end
    LOOMSTIMPARAMS(length(LOOMSTIMPARAMS)+1) = stimParams;
    si = length(LOOMSTIMPARAMS);
else % use the one we have
    stimParams = LOOMSTIMPARAMS(si);
end   
  
% Get the per facet parameters
lats = transition2responseLatency(stimParams.transition_durations, 'vc', 'monitor'); %gets latencies as a function of luminance change duration
facet_mid_times = (stimParams.transition_start_times(:) + stimParams.transition_durations(:)./2)'; 
facet_times = (stimParams.transition_start_times(:) + lats(:))';
% Put in a single facet timing jitter that is dependent on speed of luminance change, value is based on linear fit 
% of jitter to LGMD VC single facet data 
if constExcJitter
    facet_time_std = zeros(size(facet_times)) + exc.tjitter; %alternative, constant jitter for excitatory inputs
else %the jitter depends on the transition duration - this is the relationship that fits the data
    facet_time_std = stimParams.transition_durations .* .192 + exc.tjitter;
end
facet_synweights = speed2lgmdResponse(stimParams.RFspeeds, 'vc'); %the per facet weighting based on speed
facet_cmpts = getCmptsFromDegrees(stimParams.RFpos(:,1), stimParams.RFpos(:,2)); %the location (compartment) of each synapse

%check for cmpts that are out of the range of angles mapped on the dendrites and discard
nn = ~isnan(facet_cmpts); 
if (~nn) disp('Eliminating some facet inputs that fall outside of the range mapped onto dendrites'); end
facet_cmpts = facet_cmpts(nn); 
facet_synweights = facet_synweights(nn); 
facet_times = facet_times(nn); facet_time_std = facet_time_std(nn);
facet_mid_times = facet_mid_times(nn);

% now get the per synapse parameters - expand into matrices by the number of synapses per facet
[exc_synweights, ~] = meshgrid(facet_synweights, ones(exc.nspf,1));
[exc_cmpt_vec, ~] = meshgrid(facet_cmpts,ones(exc.nspf,1));
[exc_timing_vec, ~] = meshgrid(facet_times, ones(exc.nspf,1));
[exc_timing_std, ~] = meshgrid(facet_time_std, ones(exc.nspf,1));

%generates synapse random weight first, then scales, for constant SNR 
exc_gsyn_vec = exc_synweights .* normrnd(exc.gmax*ones(size(exc_synweights)), exc.gstd); 
% linearize the matrices and make sure there is no negative synaptic input, adjust times to simulation time
exc_gsyn_vec = rectify(exc_gsyn_vec(:)); 
exc_cmpt_vec = exc_cmpt_vec(:); 
exc_timing_vec = exc_timing_vec(:) + (simlen - 50);
facet_mid_times = facet_mid_times(:) + (simlen - 50);

% Generate the synapse timings - use options to use the same synaptic times across simulations or to generate new each time
if (~repeatExcTimings) %if we don't repeat, just generate the random vectors
    exc_timing_vec = normrnd(exc_timing_vec(:), exc_timing_std(:)); % and times
else 
    if (isfield(stimParams, 'exc_timing_vec') && ~isempty(stimParams.exc_timing_vec)) %then we already have the timings
        exc_timing_vec = stimParams.exc_timing_vec;
    else %generate the timings and save them
        exc_timing_vec = normrnd(exc_timing_vec(:), exc_timing_std(:));
        LOOMSTIMPARAMS(si).exc_timing_vec = exc_timing_vec;
    end
end

% For synaptic inhibition, what we do to generate their timings is to use the midpoint time of the 
% luminance change at each facet.  We then just give them a constant delay (with some jitter) and 
% a constant gmax (also with some jitter).  Both of these are constant, unlike excitation, and result
% in an overall inhibitory time course that is proportional to the area of the stimulus. 
if (~repeatInhTimings) %if we don't repeat, just generate the random vectors
    [inh_timing_vec, ~] = meshgrid(facet_mid_times+inh.delay, ones(inh.nspf,1));
    inh_timing_vec = normrnd(inh_timing_vec(:), inh.tjitter);
else 
    if (isfield(stimParams, 'inh_timing_vec') && ~isempty(stimParams.inh_timing_vec)) %then we already have the timings
        inh_timing_vec = stimParams.inh_timing_vec;
    else %generate the timings and save them
        [inh_timing_vec, ~] = meshgrid(facet_mid_times+inh.delay, ones(inh.nspf,1));
        inh_timing_vec = normrnd(inh_timing_vec(:), inh.tjitter);
        LOOMSTIMPARAMS(si).inh_timing_vec = inh_timing_vec;
    end
end
inh_gsyn_vec = rectify(normrnd(inh.gmax*ones(size(inh_timing_vec)), inh.gstd));
inh_cmpt_vec = linspace(inh.cmpt(1), inh.cmpt(end), length(inh_timing_vec)); %even inhibitory spacing over proximal dendritic trunk
inh_cmpt_vec = inh_cmpt_vec(randperm(length(inh_cmpt_vec))); %randomize across the inputs

%get back the relative to collision time rather than sim time
exc_times = exc_timing_vec - (simlen - 50); 
inh_times = inh_timing_vec - (simlen - 50);
pb=0; %boolean for plotting timecourses
if (pb) %plot if wanted    
    hold on; line(exc_times, exc_gsyn_vec, 'Parent', syn_ah, 'Marker','.', 'Color', 'k', 'LineStyle', 'none');
    hold on; line(inh_times, inh_gsyn_vec, 'Parent', syn_ah, 'Marker','.', 'Color', 'r', 'LineStyle', 'none');
    % a little fancier, do the sum over time
    bin_size = 5; %ms
    edges = (min(exc_times)-1):bin_size:(max(exc_times)+1);
    centers = edges(1:end-1) + bin_size/2;
    exc_hist = makeWeightedHist(edges, exc_times, exc_gsyn_vec);
    inh_hist = makeWeightedHist(edges, inh_times, inh_gsyn_vec);
    hold on; line(centers, exc_hist/max(exc_hist),'Parent',syn_ah2,'Color','k','Marker','none');
    hold on; line(centers, inh_hist/max(inh_hist),'Parent',syn_ah2,'Color','r','Marker','none'); 
end

% we have the hoc files for the neuron simulations, but for the analysis,
% it's nice to have a parameter file that is more easily deciphered by matlab
loomStimParameters = LOOMSTIMPARAMS;
save('simout/loomingParameters.mat', 'loomStimParameters');
clear loomStimParameters; 

%also, let's save the synaptic input sequence in order to look at it later.
global synapticInputs;
inputs.loverv = loverv;
inputs.exc_times = exc_times;
inputs.inh_times = inh_times;
inputs.exc_gsyn_vec = exc_gsyn_vec;
inputs.inh_gsyn_vec = inh_gsyn_vec;
if ~isempty(synapticInputs)
    synapticInputs(length(synapticInputs)+1) = inputs;
else
    synapticInputs = inputs;
end
save('synapticInputs.mat', 'synapticInputs');
  
% -----------------------------------------
function cmpts = getCmptsFromDegrees(az, el)
% ------------------------------------------
% Translates visual space locations in degrees to corresponding compartments
% on the dendrites of the model.  Compartments are numbered 1:400, and 
% visual space sampled is -50-50 deg elevation and 40-140 in azimuth.  Calling for 
% compartment numbers outside of this space will cause the function to return NaNs.
  
  % compartments 5-15 along each dendrite correspond to azimuths 40:140
  % This corresponds to 10 deg per compartment
  idx_x_offs = (5 + 10*(az-40)/100);
  outofrange = (idx_x_offs < 5 | idx_x_offs > 15);
  idx_x_offs(outofrange) = NaN;

  % elevation has a 100 deg span too, spread amoung 20 dendrites, or 5 deg  
  dend_ang = -50 + (0:19)*5; %the minimum angle for each dendrite
  dend_ang_max = dend_ang + 5; %the max
  
  for j=1:length(el) %find the proper dendritic branch of each elevation
      fi = find(el(j) >= dend_ang & el(j) <= dend_ang_max);
      if (~isempty(fi))
          dend_i(j) = fi(1);
      else
          dend_i(j) = NaN;
      end
  end
  
  % This computation should return numbers within 1-400, where each of twenty branches has 20 
  % compartments, and those numbered 5-15 along their length each take a 10x5 deg (az x el) space of the visual field.
  % In neuron, compartments are numbered continuously, so the idx_x_off fractional part gives the location within the comp.
  cmpts = (dend_i(:) -1)*20 + idx_x_offs(:);
  cmpts = cmpts';

  
% ---------------------------------------------------------------------
function iclampstr = generate_iclampstr(location, onset, dur, mag)
% ---------------------------------------------------------------------------
% This is produces the string to add a current clamp electrode for current injection for the model.

iclampstr = [sprintf('  objref stim\n'),  ...
            sprintf('  seg[%d] stim = new IClamp(%f)\n', floor(location), location-floor(location)), ...
            sprintf('  stim.del = %g\n',  onset), ....
            sprintf('  stim.dur = %g\n',  dur), ...
            sprintf('  stim.amp = %g\n',  mag), ...
            sprintf('  \n')];
            
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
        ii = find(bin_edges >= vals(j) , 1, 'first'); %find the time bin it's in
        if ~isempty(ii) hist_vect(ii-1) = hist_vect(ii-1) + weights(j); end %add the weight
end

