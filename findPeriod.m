function [T, autocorr, lags] = findPeriod(x,Ts)
% Gives the time period of a function
[autocorr,lags] = xcorr(x,x);

% find the period of the signal
[~,locs] = findpeaks(autocorr);


difference = diff(locs);
T = mean(difference(3:end-2))*Ts;

end