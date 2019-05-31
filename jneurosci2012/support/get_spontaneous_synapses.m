%
%  returns a 'synblock' as well as the next synapse index for spontaeous synaptic activity
%
function [synblock syn_idx] = get_spontaneous_synapses(spont_exc, spont_inh, simlen, syn_idx)
  % --- ruleZ
	inh_cmpt = spont_inh.cmpt;
	exc_cmpt = spont_exc.cmpt;

	synblock = '';

	pre_time = 10;
	post_time = 10;

	simlen = simlen-(pre_time + post_time);

  % --- excitatory block
  n_exc = ceil(simlen*spont_exc.nsyns*spont_exc.freq);
	timing_vec = rand(1,n_exc)*simlen;
	cmpt_vec = ceil(rand(1,n_exc)*max(exc_cmpt));
  synblock = [synblock generate_synblock(timing_vec+pre_time, cmpt_vec, spont_exc.gmax, spont_exc.tau_syn, spont_exc.erev, 0, syn_idx, [])];
	syn_idx = syn_idx + n_exc;

	% --- inhibitory block
  n_inh = ceil(simlen*spont_inh.nsyns*spont_inh.freq);
	timing_vec = rand(1,n_inh)*simlen;

	cmpt_vec = linspace(inh_cmpt(1), inh_cmpt(end), n_inh);
	cmpt_vec = cmpt_vec(randperm(length(cmpt_vec)));
  synblock = [synblock generate_synblock(timing_vec+pre_time, cmpt_vec, spont_inh.gmax, spont_inh.tau_syn, spont_inh.erev, 0, syn_idx, [])];
	syn_idx = syn_idx + n_inh;
  
