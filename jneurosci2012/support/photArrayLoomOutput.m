function mparams = photArrayLoomOutput(loverv)
% mparams = photArrayLoomOutput(loverv)
%
% Generates the paramters about the luminance changes occuring across an array of facets on the compound eye
% in response to a looming stimulus with looming parameter 'loverv'.

pb = 1; %plotting boolean
dbg = 0; % debug flag - enables extra plotting that was useful for debugging
dt = 5; %ms, sample rate of the traces, refresh rate of the stimulus
t = -2500:dt:80; %ms, time range - 0 ms is collision time
degperpx = .16; % spatial resolution

% just need luminance in fractional intensities
output_lum = (0:63) ./ 63; % 6 bit values to correspond with experiments
[stim_movie, mov_t, theta, vel]= generateLoomMovie(t, loverv, degperpx, output_lum); % generate the loom movies
mov_size = size(stim_movie);
global_weights = zeros(size(stim_movie,1), size(stim_movie,2));

% photoreceptor RF parameters
acceptanceAngle = 3; %the full width at 95% response decay, deg
dphi = 1.6; %the separation of adjacent sampling axes, deg - also the same as a single px size on LG monitor
maxRFwidth = 5/degperpx; %if RF size is variable, this is the maximum

% xc = 90/degperpx;
% yc = 90/degperpx;
xc = floor(mov_size(1)/2); % center the array in the stimulus, in px
yc = floor(mov_size(2)/2);
max_theta_deg = 50; %maximum half-size in deg
min_theta_px = (90-max_theta_deg)/degperpx;
max_theta_px = (90+max_theta_deg)/degperpx; %maximum half-size in px

eye_center = [90 0];
acuteSamp = 1;
variableSizeRF = 0;
if(acuteSamp) %visual sampling according to the reconstructed eye in Krapp and Gabbiani (2004).
    sampling_path = './acute_synmap.mat';
    load(sampling_path);
    centers_deg = [synmap(:,3) synmap(:,4)]; %RF centers, in deg
    x_centers = (centers_deg(:,1)-eye_center(1))/degperpx + xc; %xcenters in px, offset by 90 deg since that is the center of stimulus
    y_centers = (centers_deg(:,2)-eye_center(2))/degperpx + yc; %ycenters in px
    %find the facets whose centers are in the specified range of the stimulus to eliminate some facets to cut memory needs - the sampling map is too big for stims
    x_inrange = find(x_centers > min_theta_px & x_centers < max_theta_px);
    y_inrange = find(y_centers > min_theta_px & y_centers < max_theta_px);
    inrange = intersect(x_inrange, y_inrange);
    x_centers = x_centers(inrange); y_centers = y_centers(inrange);
    centers_deg = centers_deg(inrange,:); 
    nRFs = length(x_centers);
    
else %uniform visual sampling
    nrows = 80;
    ncolumns = 80;
    nrows = 40;
    ncolumns = 1;
    nRFs = nrows*ncolumns;
    centers_deg = makeHexGrid(nrows, ncolumns, dphi); %in deg, but starts at zero
    x_centers = centers_deg(:,1)./degperpx + xc - dphi/degperpx*floor(nrows/2); %in px
    y_centers = centers_deg(:,2)./degperpx + yc - max(centers_deg(:,2))./degperpx/2;
    % do the same as above, cut out RFs outside the stimulus area
    x_inrange = find(x_centers > min_theta_px & x_centers < max_theta_px);
    y_inrange = find(y_centers > min_theta_px & y_centers < max_theta_px);
    inrange = intersect(x_inrange, y_inrange);
    x_centers = x_centers(inrange); y_centers = y_centers(inrange);
    centers_deg = [x_centers, y_centers-yc]*degperpx;
    nRFs = length(x_centers);
end

if pb
    t_subsamp = 1; %temporal subsampling factor of luminance traces
    rf_subsamp = 10; % plot every so many RFs
    figure; lum_ah = axes; hold on;
end
lum_t = zeros(size(stim_movie, 3), length(x_centers));
trans_times = NaN*zeros(nRFs, 1);
trans_starts = NaN*zeros(nRFs, 1);
tt_samps = NaN*zeros(nRFs, 1);
RFwidths = NaN*zeros(nRFs, 1);
for i=1:nRFs
    % set the RF sizes either using a set overlap with the closest RF or a constant acceptance angle
    if (variableSizeRF) %find the spacing in between the facets
        RFwidths(i) = findClosestPt([x_centers(i), y_centers(i)], cat(2, x_centers, y_centers))*2; % a little wider than the nearest RF center
        if RFwidths(i) > maxRFwidth
            RFwidths(i) = maxRFwidth;
        end
    else
        RFwidths(i) = acceptanceAngle/degperpx;
    end
    % computes the luminance accepted by a single photoreceptor over time
    [lum_t(:, i), RFweights] = photRF([x_centers(i) y_centers(i)], RFwidths(i), ...
                               1:size(stim_movie,1), 1:size(stim_movie,2), stim_movie);
    if pb
        if (~mod(i,rf_subsamp)) %plot only every rf_subsamp RF
            line('parent', lum_ah, 'xdata', mov_t(1:t_subsamp:end), 'ydata',  lum_t(1:t_subsamp:end, i));
        end
    end
    maxl = max(lum_t(:,i)); minl = min(lum_t(:,i));
    transi = find(lum_t(:,i) > minl & lum_t(:,i) < maxl);
    startedi = find(lum_t(:,i) < maxl); %specific to a light to dark stimulus - need to alter if want to work both ways
    if (isempty(transi) && isempty(startedi))
        trans_times(i) = 0;
        trans_starts(i) = NaN;
        tt_samps(i) = 0;
        %disp(['There was an RF with no lum change:' sprintf('Xcenter : %g  Ycenter: %g', x_centers(i)*degperpx, y_centers(i)*degperpx)]);
        %figure; plot(mov_t, lum_t(:,i));
    else
        if(isempty(transi)) % transition is too fast - no intermediate samples
            trans_starts(i) = startedi(1);
            trans_times(i) = 0;
            tt_samps(i) = 0;
        else
            trans_times(i) = mov_t(transi(end)) - mov_t(transi(1)-1); % time between extremes, in ms
            tt_samps(i) = length(transi)-1;
            trans_starts(i) = transi(1);
        end
    
        % This is to visualize the weighting of the stimulus across the visual field
        RFcorner = max([1, 1], floor([x_centers(i) y_centers(i)] - size(RFweights)/2)+1);
        RFedge_x = min(RFcorner(1) + size(RFweights,1) - 1, size(stim_movie,1));
        RFedge_y = min(RFcorner(2) + size(RFweights,2) - 1, size(stim_movie,2));
        RFsize = [RFedge_x, RFedge_y] - RFcorner + [1, 1];
        global_weights(RFcorner(1):RFedge_x, RFcorner(2):RFedge_y) = ...
            global_weights(RFcorner(1):RFedge_x, RFcorner(2):RFedge_y) + RFweights(1:RFsize(1), 1:RFsize(2));
    end
end
clear stim_movie;

if (pb & dbg) %plot the spatial distribution of RFs as a set of circles, before culling those without luminance changes
    pos = centers_deg;
    width_deg = RFwidths .* degperpx;
    figure; hold on;
    xlabel('xpos, deg'); ylabel('ypos, deg');
    for i=1:size(pos,1)
        ch = circle(pos(i,:), width_deg(i)/2, inf);
        set(ch, 'FaceColor', 'none', 'EdgeColor', 'r');
    end
    clear('pos', 'width_deg');
end

% delete any RFs that don't have a luminance change, saves memory/computation
nn = ~isnan(trans_starts);
nRFs = sum(nn);
trans_starts = trans_starts(nn);
trans_times = trans_times(nn);
tt_samps = tt_samps(nn);
lum_t = lum_t(:,nn);
RFwidths = RFwidths(nn) .* degperpx; %convert to deg from px
RFpos = centers_deg(nn,:);

if (pb & dbg) 
    %plot the spatial distribution of RFs as a set of circles
    xlabel('xpos, deg'); ylabel('ypos, deg');
    for i=1:size(RFpos,1)
        ch = circle(RFpos(i,:), RFwidths(i)/2, inf, 'k');
        set(ch, 'FaceColor', 'none');
    end
    if (variableSizeRF)
    % plot a histogram of the RF widths to see how they came out, if they're variable
    figure; hist(RFwidths, 20); xlabel('RF size, deg');
    end
end
%RFspeed = RFwidths./(trans_times/1000); % the simple way of getting speed which doesn't work so well.

% calculate RF speed using a fixed, standard RF size. If we want to factor in 
% different sized RFs to the model, than computing luminance change from a variable size then speed by a constant
% is a way to make the sizes matter. Calculating both by the same variable RF cancels the two out.
disp('Fitting Gaussian CDFs to luminance changes to determine instantaneous stimulus speed...');
RFspeed = NaN*ones(size(RFwidths)); 
for i = 1:nRFs
    %[fit_mu, fit_sigma, fit_speed, fit_trace, fit_t] = fitStimulusGaussian(t, lum_t(:,i), [], [], RFwidths(i)/3);
    [fit_mu, fit_sigma, fit_speed, fit_trace, fit_t] = fitStimulusGaussian(mov_t, lum_t(:,i), [], [], acceptanceAngle/4);
    RFspeed(i) = fit_speed;
    RFsigma(i) = fit_sigma;
    if (~mod(i,rf_subsamp)) %plot only every rf_subsamp RF
        line('parent', lum_ah, 'xdata', fit_t(1:t_subsamp:end), 'ydata',  fit_trace(1:t_subsamp:end), 'Color', 'r');
    end 
end
RFspeed_varsize = RFspeed .* RFwidths./acceptanceAngle; %this is the speed scaled using the actual RF width. Calculated for comparison. 

%output structure assignments - all the timings are in ms and the size parameters are in deg
mparams.transition_start_times = mov_t(trans_starts);
mparams.transition_durations = trans_times;
mparams.RFspeeds = RFspeed;
mparams.RFwidths = RFwidths;
mparams.RFspeeds_varsize = RFspeed_varsize;
mparams.RFpos = RFpos;
mparams.loverv = loverv;
mparams.mov_t = mov_t; % parameters for the stimulus itself
mparams.theta = theta;
mparams.vel = vel;


if (pb) % optional plotting of the RFs
    
    % 3D surface plot of the spatial sampling weight 
    figure;
    gx = linspace(0,180, size(global_weights,1));
    [X,Y] = meshgrid(gx);
    so = surf(X,Y,global_weights);
    set(so,'EdgeColor', 'none', 'LineStyle', 'none');
    clear('X','Y', 'global_weights');
    % histogram of speeds
    figure; hist(RFspeed, 40); xlabel('Speed in RF (deg/sec)');
    % double y-axis plot of transition duration and speed over time
    figure;
    [ah, h1, h2] = plotyy(mparams.transition_start_times, mparams.transition_durations,mparams.transition_start_times, RFspeed);
    set(ah(1), 'Ydir', 'reverse'); set(h1, 'Linestyle', 'none', 'marker','.'); set(h2, 'Linestyle', 'none', 'marker','.'); 
    line(mparams.transition_start_times, 3/2*mparams.RFwidths./mparams.transition_durations*1000, ...
        'Parent', ah(2), 'Marker','.', 'Color','r', 'Linestyle', 'none');
    % compare the fit SDs with the transition durations 
    figure; 
    plot(mparams.transition_start_times,mparams.transition_durations, '.b', mparams.transition_start_times, 6*RFsigma, '.r');
else
    clear('global_weights');
end


