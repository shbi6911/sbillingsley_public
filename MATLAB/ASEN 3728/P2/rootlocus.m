clear all; close all;
%define constants
g = 9.81;
Iy = 3e-4;
m = 2;

%% vertical motion system

%define target poles
zeta = 0.7;
period = 0.5;

%define gains
[k5,k6,~] = kvalues(m,zeta,period);

disp("K5 = " + string(k5));
disp("K6 = " + string(k6));

%% longitudinal motion system
%define target poles
zeta = 0.87;
period = 0.155;

%define gains
[k1,k2,k4] = kvalues(Iy,zeta,period);
k3vec = linspace(0.001,0.02,10);

disp("K1 = " + string(k1));
disp("K2 = " + string(k2));
disp("K4 = " + string(k4));

figure(); grid on; hold on;
for i = 1:length(k3vec)
    k3 = k3vec(i);
    %define A matrix
    A = [0 1 0 0; 0 0 -g 0; 0 0 0 1; (k3*k4)/Iy k3/Iy -k2/Iy -k1/Iy];
    vals = eig(A);
    scatter(real(vals),imag(vals))
end
legend(string(k3vec),'Location','eastoutside');
axis equal;
hold off;

%% functions
function [k1,k2,k4] = kvalues(inertia, zeta, period)

omega_n = pi/period;

k1 = (zeta*omega_n)*2*inertia;

x = omega_n*sqrt(1-zeta^2);
y = k1^2/(4*inertia^2);

k2 = (x^2 + y)*inertia;

z = ((2*pi)/(zeta*omega_n))*10;

k4 = (1/z) - 0.05*(1/z);

end
