%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

clear; clc;

const.rho = (6.967e-13+1.454e-13)/2;
const.cd = 2.2;                       %coefficient of drag
const.A = 2;                         %total cross sectional area
               
const.d = sqrt(2)/4;                  %distance to center, each side
const.m = 8;                          %mass of 8 kg for 6U cubesat
    
const.r = 6371;
    %const.a = const.r + [150;250;350;450;550]; %semimajor axis range in km
const.a = 6371 + 550;
const.mu = 3.986004418e5;             % grav parameter in km^3/s^2
const.T = 2.*pi.*sqrt(const.a.^3./const.mu);          %orbital period range in s
const.v = (2.*pi.*const.a)./const.T.*1000;        %velocity range in m/s

const.gamma = deg2rad(linspace(0,45,1000));         %cant angles
const.theta = deg2rad(2);
%const.R1 = const.L./cos(const.gamma);
%const.AR = 
    
const.F = 0.5.*const.rho.*const.v.^2.*const.cd.*const.A;           %drag force

gamma_1 = atan(0.05);
gamma_2 = atan(0.1);

tau = -const.F.*const.d.*sin(const.gamma).*sin(const.theta); %torque
accel = tau./(const.m.*(const.d.*cos(const.gamma)).^2);   %accel

figure();
plot(rad2deg(const.gamma),tau,'b','LineWidth',2);
line1 = xline(rad2deg(gamma_1),'r','LineWidth',2);
line2 = xline(rad2deg(gamma_2),'r','LineWidth',2);
legend([line1,line2],"5 cm height","10 cm height");
title("Torque vs. Angle");

figure();
plot(rad2deg(const.gamma),accel,'b','LineWidth',2);
%line1 = xline(rad2deg(gamma_1),'r','LineWidth',2);
%line2 = xline(rad2deg(gamma_2),'r','LineWidth',2);
%legend([line1,line2],"5 cm height","10 cm height");
yline(-((2*pi)/const.T^2));
title("Acceleration vs. Angle");