function lgmd_resps = speed2lgmdResponse(speeds, type)
% function lgmd_resps = transitionDuration2lgmdResponse(transitionDuration, type)
%
% 7/8/2010, Pete Jones
% This function is used to translate a transition duration of luminance change 
% on the monitor (from an actual looming stimulus), to a response strength resulting
% an from an LGMD input (either current clamp or voltage clamp). TYPE is a string specifying
% either 'cc' or 'vc' depending on the relationship the caller wants.

% This is the saturating function that I fit to the lower part of the relation between 
% photoreceptor slope during a projector based luminance change and the normalized response
% height (peak depolarization, OFF response) of the LGMD recordings.  Both are fit using the population
% mean values, and are made to saturate at the response level predicted by the linear fit at 250 deg/s.  The
% slope at this speed was obtained from the linear fit to photoreceptor responses to translating edges.
mmfun = @(mmparam, xdata)-1.*((mmparam(1).*xdata)./(mmparam(2) + xdata)); %michealis-menten type saturating formula
% Weibull function, 3 parameters - the advantage is that it is restricted to go to zero at zero, whereas the logistic is not
wblfun = @(p, xd) weibull_func(p, -xd); %since we expect negative inputs to the functions, but weibull_func expects positive.

 % Parameter 1 is the saturation level, and 2 is the x value of max/2
speed_slope_mmparams = [641.1321 542.7785]; % These are the values from the photoreceptor recordings - light intertrial intervals

% Weibull function fit parameters
lmc.wblp = [214.3636  0.8896];
vc.wblp = [166.9140    1.0877];
cc.wblp = [188.2623  1.3035];

ts = eval(type);

slopes = mmfun(speed_slope_mmparams, speeds);
lgmd_resps = wblfun(ts.wblp, slopes);


