%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

clear all; clc;
tic

% Subject to Change
altitudes = [400e3, 450e3, 500e3, 550e3]; % Array of altitudes [m]
masses = [50, 75, 100, 150]; % Array of masses [kg]
Cd = 2.2; % Coefficient of Drag
A = 2; % Cross-Sectional Area [m^2]

% Constants
G = 6.674e-11; % Gravitational constant
Me = 5.98e24; % Mass of Earth [kg]
Re = 6378.1e3; % Radius of Earth [m]
deorbit = 120e3; % Deorbit altitude [m]
timestep = 60 * 60 * 3; % Timestep in seconds
seconds_to_years = 60 * 60 * 24 * 365.25; % Seconds to years conversion factor

% Loop through each altitude
for alt_idx = 1:length(altitudes)
    h_initial = altitudes(alt_idx);
    
    % Create a new figure for each altitude
    figure();
    hold on;
    
    % Loop through each mass
    for m_idx = 1:length(masses)
        m = masses(m_idx);
        
        % Preallocate variables
        steps = 14610;
        velocity = zeros(1, steps);
        alt = h_initial * ones(1, steps);
        time = zeros(1, steps);
        period = zeros(1, steps);
        
        % Initial period
        period(1) = sqrt(4 * pi^2 * (Re + alt(1))^3 / G / Me);
        
        % Simulation loop
        for i = 1:steps
            if alt(i) < deorbit
                break
            end
            
            % Atmospheric density from MSIS model
            msis = MSISatmosphere1000(alt(i) / 1000); % Altitude in km
            density = msis.mass * 1000; % Convert g/cm^3 to kg/m^3
            
            % Change in period
            dP = 3 * pi * A / m * (Re + alt(i)) * density * Cd * timestep;
            period(i + 1) = period(i) - dP;
            
            % Update altitude
            alt(i + 1) = (period(i + 1)^2 * G * Me / (4 * pi^2))^(1/3) - Re;
            time(i + 1) = i * timestep;
            
            % Update velocity
            velocity(i + 1) = sqrt(G * Me / (Re + alt(i + 1)));
        end
        
        % Plot for current mass
        ind = find(alt ~= h_initial); % Only plot relevant data
        plotY = alt(ind) / 1000; % Altitude in km
        plotX = time(ind) / seconds_to_years; % Time in years
        
        plot(plotX, plotY, 'DisplayName', "Mass = " + string(m) + " kg");
    end
    
    % Customize plot for the current altitude
    title("Deorbit Altitude from " + string(h_initial / 1000) + " km")
    subtitle("Atmospheric Model: NRL MSISv00, Date Input: 2016/001")
    xlabel("Time [Years]")
    ylabel("Altitude [km]")
    yline(deorbit / 1000, 'r--', 'DisplayName', 'Deorbit Altitude');
    legend('Location','southwest');
    hold off;
end

toc
