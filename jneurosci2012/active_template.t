/**
 * This is a template file for synapse placement at a set of locations.
 * As the morphology is a free parameter, you insert it.
 *
 * These template files all work by replacing bracket-percent???percent-bracket
 *  enclosed segments with some value.  The replacement is done with a MATLAB 
 *  script, and the template is the basis for a NEURON simulation, hence the
 *  C-like commenting.  So, for instance, the first line below will replace the
 *  morphoFile string with whatever the user assigns to params(x).value where
 *  params(x).idstr is "morphoFile".  This is all done in generate_hocfiles.m,
 *  so look there if the above sounds confusing.
 *
 *     morpho: cell = <%morphoFile%>
 *  variables: Rm = <%Rm%>
 *             Ri = <%Ri%>
 *             Cm = 1.5
 *       synalpha = <%synAlpha%>
 *   sampInterval = .025 ms
 *         simLen = <%simLen%>
 *         output = <%outfilePath%>/<%outfileTag%>.out
 *  
 */ 

strdef morphologyFile, parameterFile
lineRead = 1

/* Definitions */
morphologyFile = "<%morphoFile%>.hoc"
simLen = <%simLen%> // length of simulation in ms
dt = <%dt%> // resolution of simulation
samfreq_den = <%samp_dt%>/dt // resolution of data -- 0.025/.005 means every 5th is 
                              // sampled, so 1/.025=40 kHz (remember, .025 is ms)

/* Load the cell file */
xopen (morphologyFile)

  fRm = <%Rm%> // membrane resistsance
  fCm = 1.5 // capactance
  fRi = <%Ri%>  // intracellular resistivity

  v_init = -68 // base voltage -- note that this is where the simulation STARTS, 
	             // not Erest, which will stabilize based on active conductances
							 // after some time
  tstop1 = simLen // when the simulation stops

  /* Membrane properties ... */
  Rm = fRm // Ohm cm2
  Cm = fCm    // micro farad / cm 2
  Rax = fRi    // ohm cm

  // Derive the file to write to's name
  strdef fil_name // output filename
  sprint(fil_name, "<%outfilePath%>/<%outfileTag%>.out")
  g_m = 1 / Rm // membrane condctance -- this is what NEURON uses so we convert

  /* Create the leak conductance (and hence Rm) */
  forall {
    Ra = Rax // axial resistivity = intracellular resistivity
    cm = Cm // membrane capacitiance
  }


  /* setup active conductances ; this uses the class variable you assign in 
   * the morphology file, allowing this to work independent of morpho file
   * as long as that variable is used; specifically, 
   *   1 = dendrite
   *   2 = SIZ (with KCa)
   *   3 = axon
   */
 
  // Hodgkin Huxley conductance settings (HH throughout)
  gna_axon = <%gNaAxon%>   //axonal, starting as the SIZ widens on axon side
  gkdr_axon = <%gKDrAxon%> 
  gna_siz = <%gNa%>  // concentration in the narrowing
  gkdr_siz = <%gKDr%> 
  gna_dend = <%gNaDend%>  // along the main trunk of the dendritic tree, for back-propogation, not initiation
  gkdr_dend = <%gKDrDend%> 
  

  // settings for KCa -- only in SIZ
  gcal = <%gCaL%> // conductance for calcium channel (L-type like)
  gkca = <%gKCa%> // conductance for ca-dependent K channel
	// These parameters govern calcium clearance; see Peron and Gabbiani 2009 
	//  nature neuroscience, and the accompanying biological cybernetics paper
	//  from same year
  alphaCa = 0.006 // mM * cm2 / ms / mA
  tauCa = 132 // FIXED from GAbbiani and Krapp 2006 -- extrusion time constant
  kD = .030 // mM


  // For different compartment classes (dendrite, spike-initiation-zone, axon)
	//  the conductance complement is different.  This is reflected below.
  for i=1, nseg {
    access seg[i]
    if (cmpt_class[i] == 1) {  // DEND (In future models, H-current)
        insert pas g_pas=g_m e_pas=v_init // this line creates passive properties
   } else if (cmpt_class[i] == 2) {  // SIZ (HH + KCa)
        insert pas g_pas=g_m e_pas=v_init // this line creates passive properties
  		insert HH_Na35 ena=70 gmax_HH_Na35=gna_siz
			insert HH_Kdr35 ek=-80  gmax_HH_Kdr35=gkdr_siz
			insert CaL eca=120 gmax_CaL=gcal
			insert CaInternal alpha_ca_CaInternal=alphaCa tau_ca_CaInternal=tauCa 
			insert KCa gmax_KCa = gkca kD_ca = kD ek=-80
			cai = 0
   } else if (cmpt_class[i] == 3) {  // AXON (HH only)
            insert pas g_pas=g_m e_pas=v_init // this line creates passive properties
			insert HH_Na35 ena=70 gmax_HH_Na35=gna_axon
			insert HH_Kdr35 ek=-80  gmax_HH_Kdr35=gkdr_axon
   } else if (cmpt_class[i] == 4) { // proximal dendritic trunk
           insert pas g_pas=g_m e_pas=v_init
           //insert dendNa ena=70 gmax_dendNa=gna_siz*1.3
           //insert HH_Kdr35 ek=-80  gmax_HH_Kdr35=gkdr_siz
           //insert CaL eca=120 gmax_CaL=gcal
		   //insert CaInternal alpha_ca_CaInternal=alphaCa tau_ca_CaInternal=tauCa 
	       //insert KCa gmax_KCa = gkca kD_ca = kD ek=-80
		   //cai = 0
    } else { // main dendritic trunk, create a Na channel gradient
           insert pas g_pas=g_m e_pas=v_init
           //insert dendNa ena=70 gmax_dendNa=gna_dend
           //insert HH_Kdr35 ek=-80  gmax_HH_Kdr35=gkdr_dend
           //insert HH_Na35 ena=70 gmax_HH_Na35=gna_siz/(1+1*(cmpt_class[i]-4))
           //insert HH_Kdr35 ek=-80  gmax_HH_Kdr35=gkdr_siz/(1+1*(cmpt_class[i]-4))
    }
  }  
 
  /* Synaptic block goes here -- this is generally quite complex so it is 
	 * generated by a MATLAB script and then inserted.  It generally constitutes the
	 * bulk of this script.  Synaptic properties are defined in the MATLAB code. */
	<%synBlock%>

  /* **** SETUP IS NOW DONE ; SIMULATION IS RUN BELOW **** */

  /* Prepare file for writing */
  wopen (fil_name)

  finitialize (v_init)
  fcurrent()
    
  // forall for (x,0) ri_dummy(x) = ri(x) // skip the 0 area nodes at 0 and 1
  
  wv = 0
  while (t<tstop1) {
    wv = wv + 1 
    fadvance()
    /* record the voltage at the middle of the compartment */
    if (wv >= samfreq_den) {
      wv = 0
      // Sample saveparams line: fprint(" %f %f %f %f %f %f %f\n", t, seg[133].v(.5), seg[463].v(.5), seg[320].v(.5), seg[1470].v(.5), seg[1466].v(.5), seg[1480].v(.5))
			<%saveParamsLine%>
    }
  }

  /* close the file */
  wopen()
  quit()

