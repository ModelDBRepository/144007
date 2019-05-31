 %Also, show the summed gsyn over time during the stimulus
 edges = 0:5:simlen;
 exc_weighted = zeros(1, length(edges)-1);
 inh_weighted = zeros(1, length(edges)-1);
 for j=1:length(exc_timing_vec)
     ii = find(edges >= exc_timing_vec(j) , 1, 'first'); %find the time bin it's in
     if (~isempty(ii) && ii > 1 && ii <= length(exc_weighted)) exc_weighted(ii-1) = exc_weighted(ii-1) + exc_gsyn_vec(j); end %add the weight
 end
 for j = 1:length(inh_timing_vec)
     ii = find(edges >= inh_timing_vec(j) , 1, 'first'); %find the time bin it's in
     if (~isempty(ii) && ii > 1 && ii <= length(inh_weighted)) inh_weighted(ii-1) = inh_weighted(ii-1) + inh_gsyn_vec(j); end %add the weight
 end
 figure;
 subplot(2,1,1);
 weight_t = edges(1:end-1) + 2.5 - (simlen - 50);
 plot( weight_t, exc_weighted ./ max(exc_weighted), weight_t, inh_weighted ./ max(inh_weighted));
 subplot(2,1,2);
 plot(weight_t, exc_weighted./inh_weighted);