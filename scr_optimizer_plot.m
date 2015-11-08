function [stop] = plotiter(x, optimValues, state)


stop  = false;
global onsets 
global time
global original
y = scr_model(time,onsets,x);
plot(y);
hold on
plot(original,'r')
plot((y-original).^2,'k');
hold off;
drawnow
% pause;