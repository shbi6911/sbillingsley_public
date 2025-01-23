%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

function const = getConst()
%getConst defines a struct of constants for calculating drag force on a
%CubeSat
%
    %const.rho = [2.070e-9;7.248e-11;9.518e-12;...
        %1.585e-12;(6.967e-13+1.454e-13)/2];    %density range from Vallado
    const.rho = (6.967e-13+1.454e-13)/2;
    const.cd = 2.2;                       %coefficient of drag
    const.A = 2;                          %total cross sectional area
               
    const.d = sqrt(2)/4;                  %distance to center, each side
    const.m = 8;                          %mass of 8 kg for 6U cubesat
    
    const.r = 6371;
    const.a = const.r + 550; %semimajor axis range in km
    const.mu = 3.986004418e5;             % grav parameter in km^3/s^2
    const.T = 2.*pi.*sqrt(const.a.^3./const.mu);          %orbital period range in s
    const.v = (2.*pi.*const.a)./const.T.*1000;        %velocity range in m/s
    const.gamma = deg2rad(linspace(0,10,100));        %cant angle

    const.F = 0.5.*const.rho.*const.v.^2.*const.cd.*const.A;
end