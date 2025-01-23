%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

%experiment with angled sails
clear; clc;
const = getConst();
%test a vector of angular velocities
ang_vel = 1.00E-4;
%one day time span
tspan = [0 (365*24*3600)];

%set initial state
IC = [0,ang_vel];
var = 50;

% call ode45, generate trajectory
[t,y] = ode45(@(t,state) PanelEOM(t,state,var),tspan,IC);

%plot trajectory results
figure();
plot(t,y(:,1),'LineWidth',2);
%ylim([-deg2rad(10)-0.01,deg2rad(10)+0.01]);
line1 = yline(deg2rad(5),'r','LineWidth',2);
line2 = yline(deg2rad(-5),'r','LineWidth',2);

function [state_prime] = PanelEOM(t,state,var)
    
    const = getConst();
    tau = -const.F.*const.d.*sin(const.gamma(var)).*sin(state(1));
    accel = tau/(const.m*(const.d*cos(const.gamma(var)))^2);
    state_prime = [state(2);accel];

end

function const = getConst()
    %const.rho = [2.070e-9;7.248e-11;9.518e-12;...
    %    1.585e-12;(6.967e-13+1.454e-13)/2];    %density range from Vallado
    const.rho = (6.967e-13+1.454e-13)/2;
    const.cd = 2.2;                       %coefficient of drag
    const.A = 2;                          %total cross sectional area
               
    const.d = sqrt(2)/4;                  %distance to center, each side
    const.m = 8;                          %mass of 8 kg for 6U cubesat
    
    const.r = 6371;
    %const.a = const.r + [150;250;350;450;550]; %semimajor axis range in km
    const.a = 6371 + 550;
    const.mu = 3.986004418e5;             % grav parameter in km^3/s^2
    const.T = 2.*pi.*sqrt(const.a.^3./const.mu);          %orbital period range in s
    const.v = (2.*pi.*const.a)./const.T.*1000;        %velocity range in m/s
    
    const.F = 0.5.*const.rho.*const.v.^2.*const.cd.*const.A;           %drag force
    const.gamma = deg2rad(linspace(0,10,100));        %cant angles
    %const.gamma = deg2rad(5);
end