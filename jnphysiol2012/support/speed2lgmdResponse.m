function resps = speed2lgmdResponse(speeds, type)
% function lgmd_resps = transitionDuration2lgmdResponse(transitionDuration, type)
%
% 2010, Peter W. Jones
% This function is used to translate a transition duration of luminance change 
% on the monitor (from an actual looming stimulus), to a response strength resulting
% an from an LGMD input (either current clamp or voltage clamp). TYPE is a string specifying
% 'lmc', 'cc' or 'vc' depending on the relationship the caller wants.

% This is the saturating function that I fit to the lower part of the relation between 
% photoreceptor slope during a projector based luminance change and the normalized response
% height (peak depolarization, OFF response) of the LGMD recordings.  Both are fit using the population
% mean values, and are made to saturate at the response level predicted by the linear fit at 250 deg/s.  The
% slope at this speed was obtained from the linear fit to photoreceptor responses to translating edges.
mmfun = @(mmparam, xdata)-1.*((mmparam(1).*xdata)./(mmparam(2) + xdata)); %michealis-menten type saturating formula
% logistic function, 2 parameters: 1) slope, 2) position.  Will vary from 0-1, so height is not a parameter.
sigfun = @(sigparam, xdata) 1./(1+exp((xdata-sigparam(2)).*sigparam(1))); %because of the lack of a negative in the exp, expects negative x

% It turns out that I'm gonna try to use the linear relationships between photorecteptor slope and the normalized responses, and put the 
% saturation at the stage of the speed to slope relationship.
speed_slope_mmparams = [211.9025 193.3486]; %[162.0280  203.8908]; % Parameter 1 is the saturation level, and 2 is the x value of max/2
speed_slope_mmparams = [552.5687 354.7390]; % These are the values from the new photoreceptor recordings - dark intertrial intervals
speed_slope_mmparams = [641.1321 542.7785]; % These are the values from the new photoreceptor recordings - light intertrial intervals

% sigmoidal 
cc.sigp = [0.0273 -174.8450]; 
vc.sigp = [0.0137 -142.3484];
lmc.sigp= [0.0102 -177.9359];

ts = eval(type);
slopes = mmfun(speed_slope_mmparams, speeds);
resps = sigfun(ts.sigp, slopes);
