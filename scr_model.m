function y = scr_model(x,onsets,param)
%y = scr_model(x,onsets,param)
%
%x: time 
%onsets: stimulus onsets
%params: [amplitudes for each onsets tau1 dtau]; where tau2 = tau1 - dtau;

x      = x(:);
amp    = param(1:end-3);
latency= param(end-2);
tau1   = param(end-1);
tau2   = param(end);
%y = scr_model(x,onsets,amp,tau1,dtau)
% tau2   = tau1-dtau;
maxx   = tau1 * tau2 * log(tau1/tau2) / (tau1 - tau2);  %b' = 0
maxamp = abs(exp(-maxx/tau2) - exp(-maxx/tau1));
c      =  amp/maxamp;
%
y      = zeros(length(x),1);
for e = 1:length(onsets);
    xx         = x - onsets(e)-latency;
    c          = amp(e)/maxamp;
% c = amp(e);
    y(xx>=0)   = y(xx>=0) + c * (exp(-xx(xx>=0)/tau1) - exp(-xx(xx>=0)/(tau2)));
end