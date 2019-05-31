function new_sig = resampleSignal(time_vect, sig_vect, new_time_vect)
% function new_sig = resampleSignal(time_vect, sig_vect, new_time_vect)
%
% Resample within reason without an even sample ratio.  Uses the built
% in resampling for the timeseries object, by constructing a timeseries
% object then returning only the signal, not the time vector.

% Create a timeseries object:
ts1 = timeseries(sig_vect, time_vect);

% Resample ts1 using its default interpolation method:
new_ts=resample(ts1, new_time_vect);
new_sig = new_ts.data;

if(1==0) %plot the signals if you want
    figure;hold on;
    plot(time_vect, sig_vect, 'b');
    plot(new_time_vect, new_sig, 'g');
end