function cmpts = getCmptsFromDegrees(az, el)
% ------------------------------------------
% Translates visual space locations in degrees to corresponding compartments
% on the dendrites of the model.  Compartments are numbered 1:400, and 
% visual space sampled is -50-50 deg elevation and 40-140 in azimuth.  Calling for 
% compartment numbers outside of this space will cause the function to return NaNs.
  
  % compartments 5-15 along each dendrite correspond to azimuths 40:140
  % This corresponds to 10 deg per compartment. Out of range positions get no compartment.
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
      else % if the elevation is out of range, assign nan
          dend_i(j) = NaN;
      end
  end
  
  % This computation should return numbers within 1-400, where each of twenty branches has 20 
  % compartments, and those numbered 5-15 along their length each take a 10x5 deg (az x el) space of the visual field.
  % In neuron, compartments are numbered continuously, so the idx_x_off fractional part gives the location within the comp.
  cmpts = (dend_i(:) -1)*20 + idx_x_offs(:);
  cmpts = cmpts';
