% load the saved file
load('looming_stim_params.mat');
mparams = LOOMSTIMPARAMS;
clear LOOMSTIMPARAMS;
%mparams = addLatencyParams(mparams); % The parameters necessary for the simple model
plotblue = [.3, .6, 1];
colors = {[0 1 0], [1 0 0], plotblue, [0 0 0]};

load('synapticInputs.mat'); %these are the synaptic inputs generated from 50 repetitions of each l/v value using the variability in the model

% Copying some of the static parameters over from the NEURON model generation scripts to here.
% --- Neuron-wide settings
cell.rm = 4900; % % membrane resistance, Ohms/cm^2
cell.ri = 60; % intracellular resistivity Ohm x cm
% cell.area = areaFromMorphologyFile * 1e-8; %this reads the morphology file and creates an array of compartment areas, cm^2

% excitatory, mimic nicotinic acetylcholine receptors
cell.vis_exc.gmax = .6;%.675; % peak conductance in microSiemens; 1.2 with nspf 8 was a good pair
cell.vis_exc.erev = 0; % reversal potential
cell.vis_exc.tau_syn = 0.3; % rise time constant
cell.vis_exc.nspf = 6; % number of synapses per "facet"
cell.vis_exc.gstd = cell.vis_exc.gmax/5; %cell.vis_exc.gmax/5; %cell.vis_exc.gmax/10; %variability on the gmax, std of a normal distribution (gmax being the mean)
cell.vis_exc.tjitter = 6; %6;%ms, the timing jitter of visual excitatory synapses.  Times drawn from a normal distribution with tjitter STD
cell.vis_exc.cmpt = 1:400; % compartment numbers for the synapses

% inhibitory
cell.vis_inh.gmax = .0095; %.007; %was .1 - at low values I was using .0025
cell.vis_inh.erev = -75;
cell.vis_inh.tau_syn = 3;
cell.vis_inh.nspf = 4; % number of synapses per "facet"
cell.vis_inh.gstd = cell.vis_inh.gmax/20; %cell.vis_inh.gmax/20; %variability on the gmax, std of a normal distribution (gmax being the mean)
cell.vis_inh.tjitter = 10;%10; %ms, the timing jitter of visual inhibitory synapses.  Times drawn from a normal distribution with tjitter STD
cell.vis_inh.delay = 70; %ms, latency delay relative to the luminance change at a facet
cell.vis_inh.cmpt = 476:480; % compartment numbers for the synapses

logfun = @(p,x) rectify(p(1).* log(p(2).*x));