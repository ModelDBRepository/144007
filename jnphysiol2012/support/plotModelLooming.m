function fh = plotExtracellularLooming(time, ifr, stim_time, stimulus, varargin)
% fh = plotExtracellularLooming(time, ifr, varargin)
%
% Varargin args - xlimits(2 element vector), spike counts (2D matrix, 1:trialtype 2:trialnumber)
% This function is to plot the looming responses for the 9 trial type looming 
% experiment with 3 types of looming stimulus and 3 looming speeds. 
% The varargin is for the spike rasters if they are to be included.  This function 
% could be used for multiple trials in a single experiment or for population data.
% The ifr that it takes is a 3D matrix with multiple ifr's.  Dimensions are trial type, time, and repetition.  
% Variability of this ifr is also plotted.

% This colorlist is by trial type, and already has the appropriate colors lightened
colorlist = {[0.3 0.6 1], [0.65 0.8 1],  [1 0 0], [1 .5 .5], [0 .8 0], [.5 .9 .5], [0 0 0], [0 0 0], [0 0 0]};
loverv_labels = {'l/|V| = 10ms', 'l/|V| = 40ms', 'l/|V| = 80ms'};
loverv = [10 40 80];
box_positions = [40 45 20 25 10 15]; 
box_positions = [37 40 17 20 7 10 43 23 13]; 
plot_type = 'bars'; % could equal 'bars', produces either depending on the value
% The optional argument is the x limits
if ~isempty(varargin) % give x limits
    xl = varargin{1};
else
    xl = [nanmin2(time) nanmax2(time)];
end

if (length(varargin) >= 2) % give spike count data
    spk_counts = varargin{2};
    plot_spike_counts = 1;
else
    plot_spike_counts = 0;
end

if (length(varargin) >= 3) % designate a gray area to indicate spike count area
    highlightArea = varargin{3};
else
    highlightArea = [];
end
if (length(varargin) >= 4) % giving spike vectors to plot
    spks = varargin{4};
    plot_spk_rasters = 1;
else
    spks = [];
    plot_spk_rasters = 0;
end


ti = find(time >= xl(1) & time <= xl(2));
mean_traces = nanmean2(ifr(:,ti,:), 3);
std_traces = nanstd2(ifr(:,ti,:), 0, 3);
sem_traces = std_traces/sqrt(size(ifr,3));
time = time(ti);
if(plot_spk_rasters)
    spks = spks(:,ti,:);
end
all_peaks = [];

yl = [nanmin2(nanmin2(mean_traces))+nanmin2(nanmin2(sem_traces)) nanmax2(nanmax2(mean_traces))+nanmax2(nanmax2(sem_traces))]; %y limits

% Preallocate some handle vectors
mean_h = zeros(9,1); mean_err_h = zeros(9,1); stim_lh = zeros(9,1);
stim_h = zeros(3,1); ifr_h = zeros(3,1); label_h = zeros(3,1);
% The summary boxplot figure!
boxfig = figure('Position',[1 1 900 800]);
boxpeaks_h = subplot(2,2,1);
boxcounts_h = subplot(2,2,2);
%Create the figure and axes
fig_h = figure('Position',[1 1 900 800]);
label_h = axes('Parent', fig_h,  'Position', [0 0 1 1], 'Visible', 'off');
for i = 1:3 %loop through the 3 l/v values
    os = (3-i)*.3;
    if (plot_spike_counts) 
        pwidth=.4;
        spkcount_ax_h(i) = axes('Parent', fig_h, 'Position', [.75 .1/3+os .15 .60/3], 'Color','w');
        text('parent', label_h, 'Position', [.75 .85], 'String', 'Spike Counts');
    else pwidth = .6;
    end
    if (plot_spk_rasters)
        pheight = .3;
        raster_ax_h(i) = axes('Parent', fig_h, 'Position', [.1 .4/3+os pwidth .25/3], 'Visible','off', 'Ydir', 'reverse');
    else
        pheight = .5;
    end
        
    stim_h(i) = axes('Parent', fig_h, 'Position', [.1 .65/3+os pwidth .15/3], 'Color','w');
    ifr_h(i) = axes('Parent', fig_h, 'Position', [.1 .1/3+os pwidth pheight/3], 'Color','w', 'TickDir', 'out');
    peakrate_ax_h(i) = axes('Parent', fig_h, 'Position', [pwidth+.15 .1/3+os .15 .60/3], 'Color','w', 'TickDir', 'out');
    text('parent', label_h, 'Position', [.55 .85], 'String', 'Peak Firing Rates');
    text('parent', label_h, 'Position', [.1 .7/3+os], 'String', loverv_labels{i});
    
    if (~isempty(highlightArea))
        ah = area(ifr_h(i), highlightArea(i,1:2), [200 200]);
        set(ah, 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    end
    tt = i;  
     %plot the normal loom
    [mean_h(i+1) mean_err_h(i+1)]= plot_err_poly(ifr_h(i), time, mean_traces(i,:)', sem_traces(i,:)', ...
        colorlist{i+6}, (colorlist{i+6}+[1 1 1])/2, 1, 2);
   
   
    % Now, plot the stimulus
    stim_lh(i) = line('Parent', stim_h(i), 'Xdata', stim_time(i,:), 'Ydata', stimulus(i,:), ...
        'Color', colorlist{i}, 'LineWidth', 2); 

    % Now plot the peak firing rates
    labels = {'chunky'; 'constant rate'; 'normal'};
    trials = [i];
    peak_ifr = squeeze(nanmax2(ifr(trials,:,:),[],2));
    zi = find(peak_ifr == 0); peak_ifr(zi) = NaN; %zero peaks are not valid peak rates
    mean_peak_ifr = nanmean2(peak_ifr, 2);
    std_peak_ifr = nanstd2(peak_ifr, 0, 2);
    sem_peak_ifr = std_peak_ifr./sqrt(size(peak_ifr,2));
    bcolors = {colorlist{tt}, (colorlist{tt+1}+[1 1 1])/2, colorlist{i+6}};
    %plot the l/v summary axis now, in box or bar form
    if (strcmp(plot_type, 'bars'))
        axes(peakrate_ax_h(i)); hold on;
        for j=1:3
            bh = bar(j, mean_peak_ifr(j));
            set(bh, 'FaceColor', bcolors{j}, 'EdgeColor', 'k');
        end
        addErrorBars(gca, 1:3, mean_peak_ifr, sem_peak_ifr, [.2 .2 .2], .2);
        set(gca, 'XTickLabel', labels);
    else % 'boxes'
        axes(peakrate_ax_h(i)); hold on;
        X = peak_ifr';
        cm = cat(1, bcolors{:});
        boxplot(X, 'colors', cm, 'Notch', 'on', 'labels', labels, 'labelorientation', 'horizontal');
    end
    all_peaks = cat(2, all_peaks, peak_ifr');
    % align things!
    set(ifr_h(i), 'Xlim', xl, 'Ylim', yl, 'TickDir', 'out')
    set(stim_h(i), 'Xlim', xl, 'Ylim', [-inf inf], 'visible', 'off');    
    
    if plot_spk_rasters
        spk_offset = 0;
        prots = [i];
        for j=1
            sh = plotSpikeRasters(raster_ax_h(i), time, squeeze(spks(prots(j),:,:))', spk_offset);
            set(sh, 'color', colorlist{prots(j)});
            spk_offset = spk_offset + size(spks,3);
        end
        set(raster_ax_h(i), 'Xlim', xl, 'Ylim', [-inf inf], 'visible', 'off');
    end
    if plot_spike_counts
        % This is the version that makes bar plots for the 
        if (strcmp(plot_type, 'bars'))
            mean_spk_counts = nanmean2(spk_counts(trials,:), 2);
            std_spk_counts = nanstd2(spk_counts(trials,:), 0, 2);
            sem_spk_counts = std_spk_counts./sqrt(size(spk_counts,2));
            %plot now
            axes(spkcount_ax_h(i)); hold on;
            for j=1:3
                bh = bar(j, mean_spk_counts(j));
                set(bh, 'FaceColor', bcolors{j}, 'EdgeColor', 'k');
            end
            addErrorBars(gca, 1:3, mean_spk_counts, sem_spk_counts, [.2 .2 .2], .2);
            set(gca, 'XTickLabel', labels);
        else %boxes
            axes(spkcount_ax_h(i));
            X = spk_counts(trials,:)';
            cm = cat(1, bcolors{:}); 
            boxplot(X, 'colors', cm, 'Notch', 'on', 'labels', labels, 'labelorientation', 'horizontal');
        end

    end
end

 % do the box plot now
axes(boxpeaks_h); 
box_positions = [37, 40, 43, 17, 20, 23, 7, 10 13];
[box_positions, sorti] = sort(box_positions);
%X = all_peaks(:,[1:2, 4:5, 7:8]);
X = all_peaks(:,[1:9]);
X = X(:,sorti);
%cm = cat(1, bcolors{:});
boxplot(X, box_positions,'notch', 'on', 'positions', box_positions, 'colors', 'k');
set(boxpeaks_h, 'xtick', [10 20 40]);
ylabel('Peak Spike Rate (Hz)', 'FontSize', 18);
if plot_spike_counts
    axes(boxcounts_h);
    sorti = [5, 6, 9, 3, 4, 8, 1, 2, 7];
    %X = spk_counts(1:6,:)';
    X = spk_counts';
    X = X(:,sorti);
    %cm = cat(1, bcolors{:}); 
    boxplot(X, box_positions,'Notch', 'on','positions', box_positions, 'colors', 'k');
    ylabel('Spike Count', 'FontSize', 18);
    set(boxcounts_h, 'xtick', [10 20 40]);
end

setPresentationDefaults(fig_h, 0);    

    
