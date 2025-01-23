%% ASEN 2012: Coding Challenge #2
% by: [Shane Billingsley]
% SID: [110231742]
% last modified: [9/16/22]
%

%% Problem 3: 
% for the final portion of the coding challenge, you are asked to run a
% Monte Carlo simulation on a set of data to see the sensitivities of each
% variable on the overall outcome, converging to a true value for the final
% value for Δv

%% Housekeeping
clc; close all

%% Carry-over from problem 2
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

%% Create random distributions
% NOTE: the randn() function is essential for this part

N = 1e7; % number of random Monte Carlo samples
% re-define arrays to Monte Carlo vectors: x = dx*randn(1,N) + x
Isp = d_Isp*randn(1,N)+Isp;
m0 = d_m0*randn(1,N)+m0;
mf = d_mf*randn(1,N)+mf;

% run simulation with N samples (compute Δv)
v = (Isp*g0).*(log(m0./mf));

%% Make a histogram to visualize the distribution
% histogram
figure(); hold on
    histogram(v,100)
    
    title('Gaussian Curve of Monte Carlo Simulation')
    xlabel('Change in Velocity (\Deltav) [m/s]')
    ylabel('Frequency of Occurances')
hold off

% compute the mean: is it close to your value from part 2?
v_bar = mean(v);

% bonus: try changing the value of N: how does this affect your results?