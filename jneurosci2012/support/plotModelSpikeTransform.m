function [dh, fith, p_fit, r2, spk_fun, outdata] = plotModelSpikeTransform(ah, t, vm, ifr)
% function [dh, fith, p_fit, r2, spk_fun, outdata] = plotModelSpikeTransform(ah, t, vm, ifr)
%
% function assumes that vm and ifr are taken at the same sampling rate
% uses a set of trials, so each vm and spk trial are row vectors where 
% each row is a trial.

% Option
linear_thresh = 0; % Makes it return a 3 parameter fit that is the linear threshold relationship, rather than a power law

if (isempty(ah))
    figure; axes;
    ah = gca;
end
xlabel('Cell Input','FontSize', 12); ylabel('Cell Output','FontSize', 12);

[~, maxi] = max(vm(:));
vm_rest = vm(200:2000, :);
vm_rest = mean(vm_rest(:));
vm = vm - vm_rest;
dt = mean(diff(t)); % the sampling rate of the vectors
bint = 5; %ms, the width of the spike counting bin
step = bint;
stepi = floor(step/dt);
bini = floor(bint/dt);

nbins = floor(length(t)/stepi); %round down to not overrun the vector
ntrials = size(vm,2);
ifr_mean = zeros(nbins, ntrials);
vm_mean = zeros(nbins, ntrials); 
vm_bin_w = .5;
%vm_bin_edges = (min(vm(:))-vm_bin_w/2):vm_bin_w:(max(vm(:))+vm_bin_w);
vm_bin_edges = 0:vm_bin_w:(max(vm(:))+vm_bin_w);
vm_bin_centers = vm_bin_edges(1:end-1) + vm_bin_w/2;
ifr_vm = zeros(length(vm_bin_centers),ntrials);
for jj = 1:ntrials %averaging over periods of time.
    for ii=1:nbins
        mini = (ii-1)*stepi + 1;
        maxi = mini + bini;
        maxi = min(length(t), maxi); %just in case there is a single element over
        ifr_mean(ii,jj) = mean(ifr(mini:maxi, jj));
        vm_mean(ii,jj) = mean(vm(mini:maxi, jj));
    end
    % after averaging a little, we now construct a regularly spaced ifr(vm) relationship for each trial
    [~, max_ifri] = max(vm_mean(:,jj));
    tifr = ifr_mean(1:max_ifri, jj); %select the rising faze of the Vm
    tvm = vm_mean(1:max_ifri, jj);
    ifr_vm_temp = NaN*zeros(length(vm_bin_centers),1);
    for ii=1:(length(vm_bin_edges)-1)
        vmi = find(tvm >= vm_bin_edges(ii) & tvm < vm_bin_edges(ii+1));
        if(~isempty(vmi))
            ifr_vm_temp(ii) = mean(tifr(vmi)); %use the mean of any ifr's in the bin
        end
    end
    %Now we must deal with any holes that may be in the vector by linearly interpreting the ifr
    nani = find(isnan(ifr_vm_temp));
    for ii = 1:length(nani)
        if (nani(ii) > 1 && nani(ii) < length(ifr_vm_temp)) % deal w/ edge cases separately
            ifr_vm_temp(nani(ii)) = (ifr_vm_temp(nani(ii)-1) + ifr_vm_temp(nani(ii)+1))./2; %average the two points on either side
        elseif (nani(ii) == 1)
            ifr_vm_temp(nani(ii)) = ifr_vm_temp(2) - diff(ifr_vm_temp(2:3));
        elseif (nani(ii) == length(ifr_vm_temp)) 
            ifr_vm_temp(end) = ifr_vm_temp(end-1) + diff(ifr_vm_temp((end-2):(end-1)));
        end
    end
    ifr_vm(:,jj) = ifr_vm_temp;
end

%dh = plot(ah, vm_bin_centers, ifr_vm);
mean_ifr_vm = nanmean2(ifr_vm,2);
sd_ifr_vm = nanstd2(ifr_vm, 0, 2);
nani = isnan(mean_ifr_vm);
mean_ifr_vm = mean_ifr_vm(~nani); % eliminate any remaining nans
sd_ifr_vm = sd_ifr_vm(~nani); sd_ifr_vm(sd_ifr_vm < 1) = 1; %set a 1Hz lower bound on SD.
vm_bin_centers = vm_bin_centers(~nani); 
%also, we only want to fit those points that are lower than the peak
[maxv, maxi] = max(mean_ifr_vm);
mean_ifr_vm = mean_ifr_vm(1:maxi);
sd_ifr_vm = sd_ifr_vm(1:maxi);
vm_bin_centers = vm_bin_centers(1:maxi);

dh = plot(ah, vm_bin_centers, mean_ifr_vm);
dh(2) = plot(ah, vm_bin_centers, sd_ifr_vm+mean_ifr_vm, ':');
dh(3) = plot(ah, vm_bin_centers, mean_ifr_vm-sd_ifr_vm, ':');

if linear_thresh %a linear threshold!
    % 3 fitted parameters
    powfitf = @(p,x)p(3).*(rectify(x(:)-p(2)).^p(1)); %parameters 1) exponent, 2) offset 3)scale
    p0 = [1 1 1];
else %power law relationship
    %2 fitted parameters
    powfitf = @(p,x)p(2).*(rectify(x(:)).^p(1)); %parameters 1) exponent, 2)scale
    p0 = [1 1];
end

powErr = @(p) sum((powfitf(p,vm_bin_centers) - mean_ifr_vm(:)).^2 ./ (sd_ifr_vm.^2)); 
options = optimset('TolFun', 1e-12,'TolX', 1e-12, 'TolCon', 1e-12, 'MaxFunEvals', 8000, 'MaxIter', 4000, 'Display', 'final', 'algorithm',  'interior-point');
if linear_thresh
    p_fit = fmincon(powErr, p0, [], [], [], [], [1, 0, 0], [1, 50, 100], [], options); %fitting using 3 params
else
    p_fit = fmincon(powErr, p0, [], [], [], [], [0, 0], [100, 100], [], options); % fitting for 2
end
yfit = powfitf(p_fit, vm_bin_centers);
fith = plot(ah, vm_bin_centers, yfit, '--k');
r2 = calcR2(mean_ifr_vm, yfit);
spk_fun = powfitf;


outdata.vm = vm_bin_centers;
outdata.ifr_mean = mean_ifr_vm;
outdata.ifr_std = sd_ifr_vm;


% -------------------------------------
function rect_vect = rectify(in_vect)
% function rect_vect = rectify(in_vect)
%
% Support function to rectify (eliminate negative values from) a vector
% ----------------------------------------
rect_vect = in_vect;
rect_vect(rect_vect < 0) = 0; 