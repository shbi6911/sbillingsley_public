%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

%experiment with angled sails
clear; clc;
const = getConst();
%test a vector of angular velocities
ang_vel = linspace(1e-14,2.5e-4,100);
%one day time span
tspan = [0 (24*3600)];
results = zeros(1,length(ang_vel));
final_results = zeros(1,length(const.gamma));
%loop over angles
for i = 1:length(const.gamma)
    %loop over angular velocities
    for j = 1:length(ang_vel)
        IC = [0,ang_vel(j)];
        var = i;
        [t,y] = ode45(@(t,state) PanelEOM(t,state,var),tspan,IC);
        results(j) = max(y(:,1));
    end
    index = find(results > deg2rad(5), 1);
    final_results(i) = ang_vel(index);
end

figure();
plot(rad2deg(const.gamma),final_results,'LineWidth',2);
%ylim([-deg2rad(10)-0.01,deg2rad(10)+0.01]);
%gamma_1 = atan(0.05);
%gamma_2 = atan(0.1);
%line1 = xline(rad2deg(gamma_1),'r','LineWidth',2);
%line2 = xline(rad2deg(gamma_2),'r','LineWidth',2);
title("Maximum Angular Velocity before Violation");
xlabel("Gamma Angle in Degrees");
ylabel("Instantaneous Angular Velocity (rad/s)");
yline(0.4e-4);
yline(0.8e-4);
yline(1.2e-4);
yline(1.6e-4);
yline(2.0e-4);

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
    const.gamma = deg2rad(linspace(0,20,100));        %cant angles
    %const.gamma = deg2rad(5);
end