================================  BASICS  =====================================

The model in this folder consists of 2 sections - one is the MATLAB
model description of the inputs to the LGMD and the second is the
generation of a set of NEURON simulation (.hoc) files that run
simulations of the responses of a model of the LGMD with simulated
synaptic inputs in response to looming stimuli.

The general approach is to run generate_variability_sims, which will
generate a large number of NEURON hocfiles, as well as shell scripts
to execute them. Run the shell scripts and make sure data appears
reasonable.

Then, use analyze_model.m and analyzeLGMDInputSpeedTuning.m to analyze.

NOTE: the following instructions are geared to unix-like operating
systems (Linux, Mac OS X). If you are running WINDOWS, you will have
to figure out a way to do adapt the following instructions (best is to
get a proper OS like Linux).

REFERENCES:

Jones PW, Gabbiani F. (2012) Logarithmic Compression of Sensory
    Signals within the Dendritic Tree of a Collision-Sensitive
    Neuron. J Neurosci, In press.
Jones PW, Gabbiani F. (2012B) Impact of neural noise on a
    sensory-motor pathway signaling impending collision. J
    Neurophysiol 107:1067–1079.

================================  STEP - BY - STEP ============================

1) install NEURON on your system & use it to compile the .mod
files. This last step is achieved by invoking from a shell the script
nrnivmodl, which is located in the NEURON install directory. This
script should be invoked from the directory where the .mod files are
located.  In the following we will call this directory SIMS_ROOT.

The call to nrnivmodl will create a directory containing an executable
version of neuron called special. The name of the directory is related
to the architecture of the system. For example on a linux machine the
name will most likely be i686 (but maybe x86-64). On a unix-mac, umac
(maybe i386). This name needs to put in the variable

    	special_dir= '../umac'; %wherever the special program created
        by %nrnivmodl is located

By default, it is set to umac. Note that this is RELATIVE TO THE
hocfiles_rootpath DIRECTORY.

2) make a directory where generated NEURON .hoc files will be put; in
generate_dirsel_sims, replace 'hocfiles' in the line:

	hocfile_rootpath = 'hocfiles'; % where hocfiles will go and
	from where simulations are run

with whatever your directory is.  Or just make a hocfiles subdirectory
in the directory where all the simulation code is (i.e., SIMS_ROOT).

3) make an output directory where the simulations will output ASCII
files to when run from NEURON.  In generate_dirsel_sims, replace
'../simout' with your directory:

		params(2).idstr = 'outfilePath';
		params(2).value = '../simout';

Note that this is RELATIVE TO THE hocfiles_rootpath DIRECTORY.  In the
default scheme, if hocfiles is in wherever all the hoc code is, then
in that same directory you would make a directory using

     	       cd SIMS_ROOT/hocfiles
     	       mkdir ../simout 

This will create a subdirectory of the parent to hocfiles_rootpath (by
default SIMS_ROOT).  Use full paths if you are not sure.

4) Configure generate_variability_sims as you wish ; definitely assign
the number of CPUs you have:

   	n_cpus = 7;

Here, 7 = (8-1) because this was done with 8 CPUs.  Use 1 instead of 7
on a 1 CPU system.  If you use multiple machines (with simulation
files on shared network drive), then you can use notation like

    n_cpus = [7 7];
    
where each number is the number of processors used on each
machine. Also note that the 'special_dir' folder becomes a cell array
of strings specifying the appropriate directory for each machine. Also
note that newer machines support hyper-threading, which lets you run 2
processes per core simultaneously. This lets more instances than
processor cores be run simultaneously.

This is only if you use the shell-generated run scripts that are
created by the call

	generate_sh_files(n_cpus, 1, hoc_idx, hocfile_rootpath, special_dir);  


5) start MATLAB; from within the directory that
generate_variability_sims.m resides, add the support directory (with
subdirectories) to the matlab path, and run generate_variability_sims.
	 
6) Go into the hocfile_rootpath directory with a shell.  Change the
*.sh file permissions to executable (e.g., 777).

	    cd SIMS_ROOT/hocfiles
	    chmod u+x *.sh
	 
7) Execute runbatchX.sh on each machine you want, with X being the
machine number:

   cd SIMS_ROOT/hocfiles
   ./runbatchX.sh

8) In analyze_model.m, point datadir to wherever your NEURON output
is, e.g.:

	datadir = ['simout' filesep];

The above is the default.

9) Run analyze_model.m, e.g.

 analyze_model('simout/', {'loom_snr_5_jit_6_exc_0.6_inh_0.0095'}); 

will analyze looming simulations (if they exist and are named such). 

10) To try a different parameter set. Go crazy.

For speed tuning analysis: 

11) In order to run the speed tuning analysis, sets of simulations
    both with and without inhibition are needed.  So, to run the
    without-inhibition simulations, change the parameter 'vis_inh' to
    'vis_noinh' on line 183 of generate_variability_sims.m. Run a
    decent set of looming simulations (>20 trials/condition, we used
    50/condition in the paper) and analyze them.
12) Take the .out files and .mat file and put them in a subdirectory
    of the simulation folder (let's call this NOINH_DIR). Also move
    the synaptic_inputs.mat file in SIMS_ROOT to NOINH_DIR.
13) Now, go and change the 'vis_noinh' parameter back to 'vis_inh',
    and generate a set of simulations of the same size with
    inhibition.  Rerun the model generation and NEURON simulations,
    analyze. Take those results and put them in another folder (let's
    say INH_DIR).
14) Call the matlab function
analyzeLGMDInputSpeedTuning('INH_DIR/loom_snr_5_jit_6_exc_0.6_inh_0.0095_analyzed.mat',
    'NOINH_DIR/loom_snr_5_jit_6_exc_0.6_inh_0.0095_analyzed.mat');
    This function should spit out all of the analysis in the paper.

======================= HOW TO GENERATE THE PAPER FIGURES ============

These instructions refer to the figures of:

Jones PW, Gabbiani F. (2012) Logarithmic Compression of Sensory
Signals within the Dendritic Tree of a Collision-Sensitive Neuron. J
Neurosci, In press.

Simulation times in NEURON for looming simulations are ~ 45 min per
simulation on a linux machine with a Xeon E5420 CPU @2.50GHz (circa
2009).  This makes large batches of looming simulations quite time
consuming to run on the above hardware, and we often spread these
simulations across 4 computers to make the overall batch simulation
times tractable. Single facet simulations are much shorter to run -
the default batch can be run in <15 min on that machine using 7 cores.

To generate most of the panels of figure 4(C,E,F) just run full
looming stimulations, using generate_variability_sims.m and NEURON,
and run the output through analyze_model.m - as detailed in the
instructions above.

Figure 5, 6, 7 - Run analyzeLGMDInputSpeedTuning as above.  This
should run fairly quickly, but may run into problems on machines where
loading the saved simulation files consume a significant portion of
the system's memory.  Each of the files with 50 simulations is ~500MB.

================================  FILE LIST ===================================

Simulation generation files:

    generate_variability_sims.m: the main file that calls most others
	  Running single facet simulations will produce plots like
	  those in Figure 6B of Jones and Gabbiani (2012B, J
	  Neurophysiol), and running looming simulations (which take
	  much longer, about ) will produce plots like Figure 4 (2012,
	  J Neurosci).
	  To produce Figure 7A (2012B) , the simulations must be run
	  with the 'constExcJitter' flag set to 1.
	  To make Figure 7B-C (2012B), make sure that the snr value on
          line 176 of generate_variability_sims is set to [2 5 10 15].
	  To make Figure 8 (2012B), use the other version of the model
	  geared towards analyzing variability and run
	  compareLoomingVariability.m
    get_spontaneous_syanpses.m: wrapper for generate_synblock that
	  deals with spontaneous synaptic activity.
    generate_synblock.m: generates the NEURON hoc notation-based
	  synapses based on synaptic timing.  Note that all the stuff
	  dealing with visual stimulation is handled in
	  generate_dirsel_sims ; this simply generates a very long
	  string of NEURON code that is then inserted into a hocfile
	  via generate_hocfiles.
    generate_hocfiles.m: this tool uses active_template.t to
	  generate a specific hoc file with specific parameters
    generate_sh_files.m: based on the # of cpus you expect to have,
	  this will generate a shell script and subordinate shell
	  scripts that run the individual NEURON hocfiles across
	  processors.
    parameter_comparison.rtf: This is a list of the major model
       parameters and how they have changed (or haven't) between this
       model and that of Peron et al (2009).

    'support' directory - contains a number of helper functions for
        generating the simulations and doing the analysis.

NEURON support files:

  LGMD_HH_Kdr_sf35.mod: K channel (Hodgkin-Huxley)
  LGMD_HH_Na_sf35.mod: Na channel (Hodgkin-Huxley)
  LGMD_KCa.mod: calcium-dependent potassium channel module
  LGMD_CaL.mod: calcium channel module
  LGMD_Ca.mod: calcium clearance module

  rake_final.hoc: the morphology file simulating LGMD as a simplified
	  "rake" (or "trident', for the more poetically inclined)
  active_template.t: template for 'active' model with conductances;
	  used by generate_hocfiles.m

Analysis files:

  analyze_model.m: the main plotting tool that reproduces the paper's
	  figure 4.  Mild discrepancies due to synaptic jitter
	  possible.
  analyzeLGMDInputSpeedTuning.m: the analysis program for the speed
        tuning analysis.  Generates the plots in figures 5,6,7 of the
        paper. Again, mild discrepencies with published results due to
        jitter possible in figures 6 and 7.
