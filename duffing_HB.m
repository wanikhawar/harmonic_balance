clc
clear
close all

% sampling time
Ts = 0.01;

% number of terms considered in the fourier series
n_terms = 31;

% known parameters
F = 2;

% ode generation tspan
tspan = 0:Ts:500;

% initial condition
y0 = zeros(2,1);

% use ODE45 to generate the data
[tval, yval] = ode45(@duffing, tspan, y0, odeset("RelTol",1e-6));

displacement = yval(:,1);
velocity = yval(:,2);

% length of the signal
L = length(tval);

% start of steady state
% truncate the first 2/5 of the signal
% determine the lag using cross correlation of the reconstructed signal and
% the original signal
% [crosscorr, lags] = xcorr(displacement, p3)
start = ceil((2/5)*length(displacement));
start = start+420;

% steady state displacement
displacement_ss = displacement(start:end);

% find the period of the signal if unknown
T = findPeriod(displacement_ss,Ts);
omega = 2*pi/T;

% length of the period in time steps
periodLength = T/Ts;

% One period from start index
OPFS = floor(start+periodLength);

% plots

figure(1)
plot(tval, displacement)
xlabel("$t$",'interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca, 'fontsize',16)


figure(2)
plot(tval, velocity)
xlabel("$t$",'interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
set(gca, 'fontsize',16)

% Phase Portrait
figure(3)
plot(displacement(start:end),velocity(start:end))
xlabel("$x$",'interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
set(gca, 'fontsize',16)

% one cycle of displacement
figure(4)
plot(tval(start:OPFS),displacement(start:OPFS))
xlabel("$t$",'interpreter','latex')
ylabel('$x$','Interpreter','latex')
set(gca, 'fontsize',16)

% get the fourier series coefficients of the data
[A,B] = fourierCoeff(displacement(start:OPFS));

% one cycle 
singlePeriod = displacement(start:OPFS);

% find the fourier series reconstruction of the cycle
fourierReconstructionDisp = zeros(1,length(singlePeriod));
t_fourier = (0:floor(periodLength))*Ts;

for i = 1:n_terms
    fourierReconstructionDisp = fourierReconstructionDisp + A(i)*cos(2*pi*(i-1)*t_fourier/T) ...
        + B(i)*sin(2*pi*(i-1)*t_fourier/T);
end


figure(5)
plot(tval(start:OPFS), displacement(start:OPFS),'b-','linewidth',2)
hold on
plot(tval(start:OPFS), fourierReconstructionDisp,'r--','linewidth',3)
xlabel("$t$",'interpreter','latex')
ylabel('$x$','Interpreter','latex')
legend('original signal', 'fourier reconstruction')
set(gca, 'fontsize',16)
 
% HARMONIC BALANCE

p1 = zeros(length(tval),1);
for i = 1:n_terms
    p1 = p1 + ((i-1)^2)*(A(i)*cos((i-1)*omega*tval)+B(i)*sin((i-1)*omega*tval));
end
p1 = -(omega^2)*p1;


p2 = zeros(length(tval),1);
for i = 1:n_terms
    p2 = p2 + (i-1)*(-A(i)*sin((i-1)*omega*tval)+B(i)*cos((i-1)*omega*tval));
end
p2 = omega*p2;

p3 = zeros(length(tval),1);
for i = 1:n_terms
    p3 = p3 + A(i)*cos((i-1)*omega*tval)+B(i)*sin((i-1)*omega*tval);
end

p4 = p3.^3;

p5 = F*cos(omega*tval);

G = [p1 p2 p3 p4];

D = G'*G;

r = (D\(G'))*p5;

figure(6)
plot(displacement(start:end),velocity(start:end),'b-','linewidth',3)
hold on 
plot(p3,p2,'r--','linewidth',2)
xlabel("$x$",'interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
set(gca, 'fontsize',16)
legend('original signal','fourier reconstruction')

