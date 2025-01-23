%% ASEN 3801, Quadrotor Simulation and Control: Task 3
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 4/14/2024
% Last Modified: 4/14/2024

clear; close all; clc % housekeeping

% Deviations to be plotted in Task 3, Problems 3 & 4
%   a. Deviation by +5 deg in roll
%   b. Deviation by +5 deg in pitch
%   c. Deviation by +0.1 rad/sec in roll rate
%   d. Deviation by +0.1 rad/sec in pitch rate

%% Begin Task 3 Main Script
% Getting minidrone parameters
    [m, g, d, km, I, nu, mu]= getParams(); % all in SI units
% Time to integrate over (constant throughout all of Task 1)
    tspan= [0, 10]; % [s]
% Tolerance settings for ode45
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% Reference height
    zE_ref= -10; % (m)

% Figure variables/settings
    % Struct containing vectors of figure numbers
        figs.T3_Q3a=88:93; % problem 3a, 5 deg. roll deviation
        figs.T3_Q3b=94:99; % problem 3b, 5 deg. pitch deviation
        figs.T3_Q3c=100:105; % problem 3c, 0.1 rad/sec roll rate deviation
        figs.T3_Q3d=106:111; % problem 3d, 0.1 rad/sec pitch rate deviation
        figs.T3_Q5_lat=112; % problem 5, locus plot for lateral gains
        figs.T3_Q5_lon=113; % problem 5, locus plot for longitudinal gains
        figs.T3_Q7_lat=114:119; % problem 7, lateral translation
        figs.T3_Q7_lon=120:125; % problem 7, longitudinal translation
        figs.T3_Q8=126:130; % question 3.8, lateral response
    % Line colors
    linFeedback_col= "#D95319"; % linear EOM w/feedback
    nonLinFeedback_col= "#0072BD"; % nonlinear EOM w/feedback
% Cell arrays of legend entries
    legendLabels= {'Nonlinear EOM', 'Linear EOM' };
    legendLabels_3D= {'Nonlinear EOM', '', '', 'Linear EOM', ...
        'Initial Position', 'Final Position'};

% Design Constraints
    tauMax=1.25; % [s] maximum time constant

%% Calculating Lateral Velocity Gain
%       K1 = roll rate gain (p)
%       K2 = roll angle gain (phi)
%       K3 = lateral velocity gain (v, in body coordinates)

% Time constants used to obtain eigenvalues
    tau1= 0.5; % [s] dominating time constant
    tau2= 0.01; % [s]
% K3 gains to plot
    n= 1000;
    K3max= 1*(10^-3); % initial (larger) K3
    K3min= -5*(10^-4); % final (smaller) K3
% Calling subfunctions
    % All eigenvalues/gains
        [eigVals, gains] = latGains(tau1, tau2, K3max, K3min, n);
    % Plotting all eigenvalues
        rootLocus(eigVals, gains(1,1:2),[K3max, K3min], 'Lateral', ...
            figs.T3_Q5_lat)
    % Negative, real eigenvalues/gains
        [realEigs, realGains]= getReals(eigVals, gains);
    % Time constants/gains that satisfy design objective
        [tau_lat, lat_gains] = goodTaus(realEigs, realGains, tauMax);
% Creating table to simplify viewing
    labels= ["Kp", "Kphi", "Kv", "Time Constants"];
latResults= table(lat_gains(:,1), lat_gains(:,2), lat_gains(:,3), tau_lat, ...
    'VariableNames', labels)

%% Calculating Longitudinal Velocity Gain
%       K1 = pitch rate gain (q)
%       K2 = pitch angle gain (theta)
%       K3 = longitudinal velocity gain (u, in body coordinates)

% Time constants used to obtain eigenvalues
    tau1= 0.5; % [s] dominating time constant
    tau2= 0.02; % [s]
% K3 gains to plot
    n= 1000;
    K3max= 1*(10^-3); % initial (larger) K3
    K3min= -1*(10^-3); % final (smaller) K3
% Calling subfunctions
    % All eigenvalues/gains
        [eigVals, gains] = lonGains(tau1, tau2, K3max, K3min, n);
    % Plotting all eigenvalues
        rootLocus(eigVals, gains(1,1:2),[K3max, K3min], 'Longitudinal', ...
            figs.T3_Q5_lon)
    % Negative, real eigenvalues/gains
        [realEigs, realGains]= getReals(eigVals, gains);
    % Time constants/gains that satisfy design objective
        [tau_lon, lon_gains] = goodTaus(realEigs, realGains, tauMax);
% Creating table to simplify viewing
    labels= ["Kq", "Ktheta", "Ku", "Time Constants"];
lonResults= table(lon_gains(:,1), lon_gains(:,2), lon_gains(:,3), tau_lon, ...
    'VariableNames', labels)

%% Task 3, Problem 3a: Deviation by +5 deg in Roll
% Given disturbance
    deltaPhi= 5*pi/180; % [rad]
% Initial conditions for disturbances
    IC= getICs(z=zE_ref, phi=deltaPhi);

% NONLINEAR CLOSED-LOOP ode45 call for 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_CL(t, var, g, m, I, nu, mu), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3a, nonLinFeedback_col)

% LINEAR CLOSED-LOOP ode45 call for 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linear_CL(t, var, g, m, I), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3a, linFeedback_col, LineStyle=...
    '-.', LegendEntries=legendLabels, LegendEntries_3D= legendLabels_3D, ...
    figLabel= 'Task 3, Problem 3a & 4a: ')

%% Task 3, Problem 3b: Deviation by +5 deg in Pitch
% Given disturbance
    deltaTheta= 5*pi/180; % [rad]
% Initial conditions for disturbances
    IC= getICs(z=zE_ref, theta=deltaTheta);

% NONLINEAR CLOSED-LOOP ode45 call for 5 deg. pitch disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_CL(t, var, g, m, I, nu, mu), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3b, nonLinFeedback_col)

% LINEAR CLOSED-LOOP ode45 call for 5 deg. pitch disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linear_CL(t, var, g, m, I), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3b, linFeedback_col, LineStyle=...
    '-.', LegendEntries=legendLabels, LegendEntries_3D= legendLabels_3D, ...
    figLabel= 'Task 3, Problem 3b & 4b: ')

%% Task 3, Problem 3c: Deviation by 0.1 rad/sec in Roll Rate
% Given disturbance
    delta_p= 0.1; % [rad/s]
% Initial conditions for disturbances
    IC= getICs(z=zE_ref, p=delta_p);

% NONLINEAR CLOSED-LOOP ode45 call for 0.1 rad/sec roll rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_CL(t, var, g, m, I, nu, mu), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3c, nonLinFeedback_col)

% LINEAR CLOSED-LOOP ode45 call for 0.1 rad/sec roll rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linear_CL(t, var, g, m, I), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3c, linFeedback_col, LineStyle=...
    '-.', LegendEntries=legendLabels, LegendEntries_3D= legendLabels_3D, ...
    figLabel= 'Task 3, Problem 3c & 4c: ')

%% Task 3, Problem 3d: Deviation by 0.1 rad/sec in Pitch Rate
% Given disturbance
    delta_q= 0.1; % [rad/s]
% Initial conditions for disturbances
    IC= getICs(z=zE_ref, q=delta_q);

% NONLINEAR CLOSED-LOOP ode45 call for 0.1 rad/sec pitch rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_CL(t, var, g, m, I, nu, mu), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3d, nonLinFeedback_col)

% LINEAR CLOSED-LOOP ode45 call for 0.1 rad/sec pitch rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linear_CL(t, var, g, m, I), ...
    tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = InnerLoopFeedback(x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T3_Q3d, linFeedback_col, LineStyle=...
    '-.', LegendEntries=legendLabels, LegendEntries_3D= legendLabels_3D, ...
    figLabel= 'Task 3, Problem 3d & 4d: ')

%% Task 3, Problem 7: Testing Simulation for Lateral Feedforward Command
% Lateral reference velocity for design objective
    v_ref= 0.5; % [m/s] 
% Initial conditions
    IC= getICs(v=v_ref, z=zE_ref);

% NONLINEAR CLOSED-LOOP ode45 call for 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_velocityRef(t, var, g, m, I, nu, ...
    mu), tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = VelocityReferenceFeedback(t, x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q7_lat, nonLinFeedback_col, ...
        figLabel='Problem 3.7, Lateral Translation: ')

%% Task 3, Problem 7: Testing Simulation for Longitudinal Feedforward Command
% Lateral reference velocity for design objective
    u_ref= 0.5; % [m/s] 
% Initial conditions
    IC= getICs(u=u_ref, z=zE_ref);

% NONLINEAR CLOSED-LOOP ode45 call for 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_velocityRef(t, var, g, m, I, nu, ...
    mu, 'longitudinal'), tspan, IC, options);
% Calculating non-constant control forces
    N= length(x);
    [Fc, Gc] = VelocityReferenceFeedback(t, x');
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T3_Q7_lon, nonLinFeedback_col, ...
        figLabel='Problem 3.7, Longitudinal Translation: ')

%% Task 3, Problem 8: Comparing Simulated Response with both Data Sets
% Control input start time (for recorded data sets)
    t_trim= 6; % [s]
% Line colors
    sim_col= "#0072BD"; % simulation line color
    dat11_49_col= "#D95319"; % data file 11_49 line color
    dat11_20_col= "#77AC30"; % data file 11_20 line color
% 2D legend entries
    legend_2D= {'Simulated Flight Path', 'Data from Submitted Gains', ...
        'Default Data'};
% 3D legend entries
    legend_3D= {'Simulated Flight Path', '','', ...
        'Data from Submitted Gains', '','', 'Default Data'};

% Lateral reference velocity for design objective
    v_ref= 0.5; % [m/s] 
% Initial conditions
    zE_ref= -1.1; % [m]
    IC= getICs(v=v_ref, z=zE_ref);

% NONLINEAR CLOSED-LOOP ode45 call for 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_velocityRef(t, var, g, m, I, nu, ...
    mu), tspan, IC, options);
% Plotting with PlotAircraftSim function
    plot_4pt8(t, x', figs.T3_Q8, sim_col)

% Plotting data from calculated control gains (file: "RSdata_11_49.mat")
    droneDataFile= 'RSdata_11_49.mat';
    [x, t] = getDataState(droneDataFile, t_trim);
    plot_4pt8(t, x', figs.T3_Q8, dat11_49_col, LineStyle='-.')

% Plotting data from calculated control gains (file: "RSdata_11_20.mat")
    droneDataFile= 'RSdata_11_20.mat';
    [x, t] = getDataState(droneDataFile, t_trim);
    plot_4pt8(t, x', figs.T3_Q8, dat11_20_col, LineStyle=':', LegendEntries ...
        =legend_2D, LegendEntries_3D=legend_3D, ...
        figLabel='Problem 3.8, Lateral Translation:')

    