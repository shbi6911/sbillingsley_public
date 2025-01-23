%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

clear all ; clc
tic
% Subject to Change
h_initial = 550e3; % Altitude [m]
m = 47 ; % Mass [kg]
Cd = 2.2; %C oefficient of Drag
A = 2 ; % Cross-Sectional Area m^2

% Constants
G = 6.674e-11;
Me = 5.98e24; % [kg]
Re = 6378.1e3; % [m]

% IMPORTANT:
% Densities from MSIS are in g/cm^3
% 1 g/cm^3 = 1000 kg/m^3
deorbit = 120e3;
timestep = 60 * 60 * 3 * 1;
timestep_days = 60*60*24;

% Steps
steps = 14610;

% Preallocate
velocity = zeros(1,steps);
alt = h_initial * ones(1, steps);
time = zeros(1,steps);

period(1) = sqrt(4*pi^2 * (Re+alt(1))^3 / G / Me);
period(2:steps) = 0;

for i = 1:steps
    if (alt(i) < deorbit)
        break
    end

    msis = MSISatmosphere1000(alt(i)/1000);
    density = msis.mass * 1000;

    % change in period
    dP = 3*pi*A / m * (Re + alt(i)) * density * Cd * timestep;
    period(i+1) = period(i) - dP;
    alt(i+1) = (period(i+1)^2 * G * Me / (4*pi^2))^(1/3) - Re;
    time(i+1) = i*timestep;

    velocity(i+1) = sqrt(G*Me/(Re + alt(i+1)));
end

seconds_to_years = 60*60*24*365.25;

ind = find(alt ~= h_initial);
plotY = alt(ind) / 1000;
plotX = time(ind)/seconds_to_years;

figure()
plot(plotX, plotY)
title("Deorbit Altitude of " + string(m) + " kg CubeSat from " + string(h_initial/1000) + " km altitude")
subtitle("Atmospheric Model: NRL MSISv00, Date Input: 2016/001")
xlabel("Time [Years]")
ylabel("Altitude [km]")
yline(deorbit/1000, 'r--') 
toc
