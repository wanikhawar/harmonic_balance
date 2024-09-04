function ydot = duffing(t,y)
% % duffing oscillator

m = 1;
c = 0.2;
k = 1;
alpha = 1;

% m = 1.0734;
% c = 2.7417;
% k = -1.9223;
% alpha = 2.4490;

F = 2;
omega = 0.4;

ydot = zeros(2,1);
ydot(1) = y(2);
ydot(2) = (1/m) * (F*cos(omega.*t)-c*y(2)-k*y(1)-alpha*y(1)^3);

end