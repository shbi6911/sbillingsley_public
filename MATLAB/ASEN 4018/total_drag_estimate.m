%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

%order of magnitude estimates for total drag force
clear; clc;
%density values from Vallado 2022
rho = [2.070e-9;7.248e-11;9.518e-12;...
    1.585e-12;(6.967e-13+1.454e-13)/2];    %density range
cd = 2.2;                       %coefficient of drag
A = 2;                          %total cross sectional area

r = 6371;
a = r + [150;250;350;450;550]; %semimajor axis range
mu = 3.986004418e5;             % grav parameter in km^3/s^2
v = sqrt(mu./a).*1000;          %orbital velocity range in m/s

F = (v.^2).*((0.5*cd*A).*rho);           %drag force


