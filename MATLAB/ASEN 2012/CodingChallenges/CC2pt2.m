%% ASEN 2012: Coding Challenge #2
% by: [Shane Billingsley]
% SID: [110231742]
% last modified: [9/16/22]
%

%% Problem 2: 
% for this problem, you'll be examining the error propagated in the
% Tsiolkovsky rocket equation with a set of given uncertainties. It should
% be noted that you MUST use the general method for this part of the coding
% challenge. If you're struggling with these partial derivatives, go ahead
% and ask a TA or a classmate for help. 
%
% Tsiolkovsky Eqn: Δv = Isp*g0*ln(m0/mf)
% note: in MATLAB, ln(x) is written as log(x)

%% Housekeeping
clc; close all

%% Declare given parameters
% given values
g0 = 9.81; % acceleration due to gravity [m/s2]
Isp = 459; % specific impulse [s]
ms = 13050; % mass of structural elements [kg]
mp = 71800; % mass of the propellant [kg]

% uncertainties
d_Isp = 11; 
d_ms = 60;
d_mp = 300; 

% dependant quantities
m0 = ms+mp;
mf = ms;

% dependant uncertainties
d_m0 = norm([d_ms d_mp]); % the norm() function is useful for adding in quadrature
d_mf = d_ms;

%% Compute best guess for Δv
v = Isp*g0*log(m0/mf);

%% Define error equations 

% define each error equation with best guess values
dvdIsp = d_Isp/Isp; % uncertainty in Δv due to specfici impulse
dvdm0 = d_m0/m0; % uncertainty in Δv due to initial mass
dvdmf = d_mf/mf; % uncertainty in Δv due to final mass

% combine these relative errors in quadrature using the general method
GM_d_m0 = ((Isp*g0)/m0)*d_m0;
GM_d_mf = ((-Isp*g0)/mf)*d_mf;
GM_d_Isp = g0*log(m0/mf)*d_Isp;
d_v = norm([GM_d_Isp,GM_d_m0,GM_d_mf]);
