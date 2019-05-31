================================  BASICS  =====================================

The general approach is to run generate_variability_sims, which will generate a
large number of NEURON hocfiles, as well as shell scripts to execute
them. Run the shell scripts and make sure data appears reasonable. 

Then, use analyse_model.m to analyze.

NOTE: the following instructions are geared to unix-like operating
systems (Linux, Mac OS X). If you are running WINDOWS, you will
have to figure out a way to do adapt the following instructions (best
is to get a proper OS like Linux). 
 
================================  STEP - BY - STEP ============================

1) install NEURON on your system & use it to compile the .mod
files. This last step is achieved by invoking from a shell the script nrnivmodl,
which is located in the NEURON install directory. This script should
be invoked from the directory where the .mod files are located.
In the following we will call this directory SIMS_ROOT.

The call to nrnivmodl will create a directory containing an executable
version of neuron called special. The name of the directory is related
to the architecture of the system. For example on a linux machine the
name will most likely be i686. On a unix-mac, umac. This name needs to put in 
the variable 

    	special_dir= '../umac'; %where the special program created by
                             	    	%nrnivmodl is located     

By default, it is set to umac. Note that this is RELATIVE TO THE
hocfiles_rootpath DIRECTORY. 

2) make a directory where generated NEURON .hoc files will be put; in 
generate_dirsel_sims, replace 'hocfiles' in the line:

	hocfile_rootpath = 'hocfiles'; % where hocfiles will go and from where simulations are run

with whatever your directory is.  Or just make a hocfiles subdirectory in
the directory where all the simulation code is (i.e., SIMS_ROOT).

3) make an output directory where the simulations will output ASCII files to
when run from NEURON.  In generate_dirsel_sims, replace '../simout'
with your directory:

		params(2).idstr = 'outfilePath';
		params(2).value = '../simout';

Note that this is RELATIVE TO THE hocfiles_rootpath DIRECTORY.  In the
default scheme, if hocfiles is in wherever all the hoc code is, then in that
same directory you would make a directory using 

     	       cd SIMS_ROOT/hocfiles
     	       mkdir ../simout 

This will create a subdirectory of the parent to hocfiles_rootpath (by
default SIMS_ROOT).  Use full paths if you are not sure.

4) Configure generate_variability_sims as you wish ; definitely assign the number of
CPUs you have:

   	n_cpus = 7;

Here, 7 = (8-1) because this was done with 8 CPUs.  Use 1 instead of 7 on a
1 CPU system.  If you use multiple machines (with simulation files on shared network drive), 
then you can use notation like 

    n_cpus = [7 7];
    
where each number is the number of processors used on each machine. Also note that 
the 'special_dir' folder becomes a cell array of strings specifying the appropriate 
directory for each machine.  

This is only if you use the shell-generated run scripts that are created by the call

	generate_sh_files(n_cpus, 1, hoc_idx, hocfile_rootpath, special_dir);  


5) start MATLAB; from within the directory that generate_variability_sims.m resides,
add the support directory (with subdirectories) to the matlab path, and run generate_variability_sims.
	 
6) Go into the hocfile_rootpath directory with a shell.  Change the *.sh file
permissions to executable (e.g., 777).  

	    cd SIMS_ROOT/hocfiles
	    chmod u+x *.sh
	 
7) Execute runbatchX.sh on each machine you want, with X being the machine number:

   cd SIMS_ROOT/hocfiles
   ./runbatchX.sh

8) In analyze_model.m, point datadir to wherever your NEURON output is:

	datadir = ['simout' filesep];

The above is the default.

9) Run analyze_model.m, e.g.

 analyze_model('simout/', {'singlefacet_gain_1'}); 

will analyze single facet simulations (if they exist and are named such). 

10) To try a different parameter set  

================================  FILE LIST ===================================

Simulation generation files:

    generate_variability_sims.m: the main file that calls all others -- really, this
	  is all you should need to look at.

    get_spontaneous_syanpses.m: wrapper for generate_synblock that deals with
	  spontaneous synaptic activity.

    generate_synblock.m: generates the NEURON hoc notation-based synapses based
	  on synaptic timing.  Note that all the stuff dealing with visual stimulation
		is handled in generate_dirsel_sims ; this simply generates a very long
		string of NEURON code that is then inserted into a hocfile via 
		generate_hocfiles.

	generate_hocfiles.m: this tool uses active_template.t to generate a specific
	  hoc file with specific parameters

    generate_sh_files.m: based on the # of cpus you expect to have, this will 
	  generate a shell script and subordinate shell scripts that run the
		individual NEURON hocfiles across processors.

    parameter_comparison.rtf: This is a list of the major model parameters and how they have 
       changed (or haven't) between this model and that of Peron et al (2009).

    'support' directory - contains a number of helper functions for generating the simulations and doing
        the analysis.

NEURON support files:

  LGMD_HH_Kdr_sf35.mod: K channel (Hodgkin-Huxley)
  LGMD_HH_Na_sf35.mod: Na channel (Hodgkin-Huxley)
  LGMD_KCa.mod: calcium-dependent potassium channel module
  LGMD_CaL.mod: calcium channel module
  LGMD_Ca.mod: calcium clearance module

	rake_final.hoc: the morphology file simulating LGMD as a simplified "rake" 
	  (or "trident', for the more poetically inclined)

	active_template.t: template for 'active' model with conductances; used by
	  generate_hocfiles.m

Analysis files:

  analyze_model.m: the main plotting tool that reproduces the paper's figures
	  6, and 7.  Mild discrepancies due to synaptic jitter possible.  

  compareLoomingVariability.m: This file generates the bar graphs of figure 8.  
        compareLoomingVariability.mat is the condensed data file for this function. It has the 
        averaged data (.mat) from 10 sets of 50 trials (over 3 l/v values) to estimate variability when 
        eliminating individual sources. This many trials takes a long time to run (even over 4 machines
        if you want to recreate these results - but computation speeds are always increasing). 

  averageModelResultStructs.m: This averages multiple .mat files produced by analyze_model.m and returns
        them in a format that when saved can be used by compareLoomingVariability.m
  
  





