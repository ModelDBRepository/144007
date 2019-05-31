% Modified by PWJ on 12/30/2010
%
% This will generate a synapse activation text block that can be plugged
% directly into a NEURON simulation given a set of synaptic positions
% and timings, and conductances (three vectors of equal size).  The other variables are the 
% proerties of the synapse.
%
% to call multiple blocks for a single simulation, pass the first_syn var
% so that synapse indexing is incremented
%
% cmpt_dens -- if assigned, will scale gmax by cmpt_dens(x) for compartment x
function synblock = generate_synblock(timing_vec, cmpt_vec, gmax, tau_syn, erev, jitter_max, first_syn, cmpt_dens)
	sidx = 1;

  disp(['Generating ' num2str(length(timing_vec)) ' synapses ...']);

	% The main loop
	for t=1:length(timing_vec)
        if (length(cmpt_dens) > 0)
            gsf = cmpt_dens(floor(cmpt_vec(t)));
        else
            gsf = 1;
        end
        if (length(gmax) > 1) %if gmax is a vector then use this string generation
            sblock{sidx} = [sprintf('  objref syn_%d\n', t+first_syn), ...
                sprintf('  seg[%d] syn_%d = new AlphaSynapse(%f)\n', floor(cmpt_vec(t)), t+first_syn, cmpt_vec(t)-floor(cmpt_vec(t))), ...
                sprintf('  syn_%d.onset = %g\n', t+first_syn, timing_vec(t)), ....
                sprintf('  syn_%d.tau = %g\n', t+first_syn, tau_syn), ...
                sprintf('  syn_%d.gmax = %g\n', t+first_syn, gmax(t)*gsf), ...
                sprintf('  syn_%d.e = %d\n', t+first_syn, erev), ...
                sprintf('  \n')];
        else % use one that doesn't index gmax
            sblock{sidx} = [sprintf('  objref syn_%d\n', t+first_syn), ...
                sprintf('  seg[%d] syn_%d = new AlphaSynapse(%f)\n', floor(cmpt_vec(t)), t+first_syn, cmpt_vec(t)-floor(cmpt_vec(t))), ...
                sprintf('  syn_%d.onset = %g\n', t+first_syn, timing_vec(t)), ....
                sprintf('  syn_%d.tau = %g\n', t+first_syn, tau_syn), ...
                sprintf('  syn_%d.gmax = %g\n', t+first_syn, gmax*gsf), ...
                sprintf('  syn_%d.e = %d\n', t+first_syn, erev), ...
                sprintf('  \n')];
        end
        sidx = sidx + 1;
	end

% Now generate file output ... 	
synblock = char(sblock);
a = synblock';
synblock = a(:)';
  
