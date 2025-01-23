%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/25/2024
% Last Modified: 3/25/2024

% The getParams function outputs the minidrone parameters given in the lab
% document "ASEN 3801 Quadrotor Simulation and Control". The function takes
% no input, and is meant to simplify defining constants over various
% scripts.

% INPUTS: none
%
% OUTPUTS
% m: minidrone mass (kilograms)
% g: acceleration of gravity ( meters/(second^2) )
% d: distance between motors and center of mass (meters)
% km: control moment coefficient ( Newton*meters/Newtons )
% I: 3x3 matrix containing diagonal values of the mass moments of inertia
%       about the x, y, and z body axes ( kilograms*(meter^2) ). The
%       minidrone is assumed to be symmetric about its body axes, so the
%       off-diagonal values are zero.
% nu: aerodynamic force coefficient ( Newtons/(meters/(second^2)) ). Defined 
%       as 1/2 of the product of air density, drag coefficient, and 
%       cross-sectional area normal to the air-relative velocity vector.
% mu: aerodynamic moment coefficient ( Newton*meters/(radians/(second^2)) )

%% Begin Function Definition
function [m, g, d, km, I, nu, mu]= getParams()
    m= 0.068; % [kg] minidrone mass
    g= 9.81; % [m/s^2] gravitational acceleration
    d= 0.060; % [m] Radial distance from CG to propeller
    km= 0.0024; % [N*m/(N)] Control moment coefficient 
    Ix= 5.8*(10^-5); % [kg*m^2] Body x-axis Moment of Inertia
    Iy= 7.2*(10^-5); % [kg*m^2] Body y-axis Moment of Inertia
    Iz= 1.0*(10^-4); % [kg*m^2] Body z-axis Moment of Inertia
        I= diag([Ix, Iy, Iz]); % [kg*m^2] Moment of inertia matrix
    nu= 1*(10^-3); % [N/(m/s)^2] Aerodynamic force coefficient
    mu= 2*(10^-6); % [N*m/(rad/s)^2] Aerodynamic moment coefficient
end