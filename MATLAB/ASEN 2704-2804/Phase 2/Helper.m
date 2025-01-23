clc; clear; close all;

W = [0.163 0.166 0.166 0.182 0.235 0.235 0.289 0.257 0.257]; % Weight [kg]
h = 17.5; % Launch height [m]
R = [15.4 17.1 18.1 18.5 15 16 15 31.8 29.7]; % Flight range [m]
tf = [4.6 4.75 4.53 2.57 3.28 3.16 3.94 3.83 3.21]; % Flight time [s]
rho = 1.17723; % Atmospheric density [kg/m^3]
S = (18490.5*2)*10^-6; % Planform area [m^2]
a = 0.2132; 

for i = 1:length(W)
V(i) = sqrt((h^2)+(R(i)^2))/tf(i);
LD(i) = R(i)/h;
CL(i) = (W(i))/(0.5*rho*S*(V(i)^2));
CD(i) = CL(i)/LD(i);
alpha(i) = (CL(i)/a)+a;
end
