%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/23/2024
% Last Modified: 4/13/2024

close all; clear; clc % housekeeping

%% Begin Task 1
% Getting minidrone parameters
    [m, g, d, km, I, nu, mu]= getParams(); % all in SI units
% Time to integrate over (constant throughout all of Task 1)
    tspan= [0, 10]; % [s]
% Tolerance settings for ode45
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% Figure variables/settings
    % Struct containing vectors of figure numbers
        figs.Q1pt2=1:6; % question 1.2
        figs.Q1pt3=7:12; % question 1.3
        figs.Q1pt4_yaw0=13:18; % question 1.4 w/ psi=0 deg
        figs.Q1pt4_yaw90=19:24; % question 1.4 w/ psi=90 deg
        figs.Q1pt5_data=25:27; % question 1.5, Minidrone data plots
        figs.Q1pt5=28:33; % question 1.5, stability simulation plots

    nonlin_col= "#0072BD"; % line color for nonlinear EOM
    drone_col= "#0072BD"; % line color for Minidrone data plots

%% Task 1, Question 2: Steady, Hovering Flight (No Aerodynamic Forces)
% Control forces for Steady, Level, Hovering flight
    [fc_4x1, motor_forces] = openLoopForces(m, g, d, km);

% Initial conditions
    IC = getICs(z=-1.1); % (m) using same initial height as minidrone

% To remove aerodynamic forces without changing the EOM function, the
% aerodynamic and moment coefficients (nu and mu) are set to zero in the
% function call
% Calling ode45 without aerodynamic forces
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, 0, 0, ...
    motor_forces), tspan, IC, options);

% Verification Plots for Hovering Trim State
    N= length(t); % time steps from ode45 integration
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);

% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x', fc_4xN, figs.Q1pt2, nonlin_col, figLabel= ...
        'Task 1, Problem 2: ')

%% Task 1, Question 3: Steady, Hovering Flight (With Aerodynamic Forces)
% Aerodynamic forces and moments are now to be accounted for, so mu and nu
% are returned to the valued given in the lab document. The trim condition
% is the same as in Task 1, Q1.2, so the control forces and initial
% conditions remain the same.

% Calling ode45 with aerodynamic forces
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu,...
    motor_forces), tspan, IC, options);

% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x', fc_4xN, figs.Q1pt3, nonlin_col, figLabel= ...
        'Task 1, Problem 3: ')

%% Task 1, Question 4: Trimmed Flight at 5m/s East & Yaw = 0 deg
    psi= 0; % [rad] zero degree/radian yaw angle
% Solving derived equations for new trim thrust and euler angles
    V_E= [0; 5; 0]; % [m/s] given inertial trim velocity 
    D= nu*(norm(V_E)^2); % [N] magnitude of drag force
    phi= atan(D/(m*g)); % [rad] roll angle req'd to oppose drag force

% Control forces required for trim at calculated roll angle
[yaw0_fc_4x1, yaw0_motor_forces] = openLoopForces(m, g, d, km, angle=phi);

% Initial conditions
    P0_E= [0; 0; -1.1]; % [m] inertial position
    euler0= [phi; 0; psi]; % [rad] euler angles (roll,pitch, yaw)
    V0_B= R_E2B(euler0)*V_E; % [m/s] velocity rotated into body coordinates
    ICs_1pt4 = vec2ICs(P_E=P0_E, O=euler0, V_B=V0_B);

% Calling ode45 with aerodynamic forces, V_E=5m/s East, psi=0deg
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu,...
    yaw0_motor_forces), tspan, ICs_1pt4, options);

% Verification Plots for Trim State at 5m/s East, psi=0deg 
    N= length(t); % time steps from ode45 integration
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= yaw0_fc_4x1*ones(1,N);

% Defining plot limits manually based on observation
    zE_lims= [1, 1.2]; % [m]
    v_lims= [4.4, 5.5]; % m/s
    w_lims= [-0.2, -.15]; % [m/s]

% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x', fc_4xN, figs.Q1pt4_yaw0, nonlin_col, ...
        zE_ylim=zE_lims, v_ylim=v_lims, w_ylim=w_lims, zE_3D_zlim=zE_lims, ...
        figLabel= 'Task 1, Problem 4, \psi=0^{\circ}: ')

%% Task 1, Question 4: Trimmed Flight at 5m/s East & Yaw = 90 deg
    psi= 90*(pi/180); % [rad] 90 degree yaw angle
% Solving derived equations for new trim thrust and euler angles
    V_E= [0; 5; 0]; % [m/s] given inertial trim velocity 
    D= nu*(norm(V_E)^2); % [N] magnitude of drag force
    theta= -atan(D/(m*g)); % [rad] pitch angle req'd to oppose drag force

% Control forces required for trim at calculated roll angle
[yaw90_fc_4x1, yaw90_motor_forces] = openLoopForces(m, g, d, km, angle=theta);

% Initial conditions
    P0_E= [0; 0; -1.1]; % [m] inertial position
    euler0= [0; theta; psi]; % [rad] euler angles (roll,pitch, yaw)
    V0_B= R_E2B(euler0)*V_E; % [m/s] velocity rotated into body coordinates

  ICs_1pt4 = vec2ICs(P_E=P0_E, O=euler0, V_B=V0_B);

% Calling ode45 with aerodynamic forces, V_E=5m/s East, psi=0deg
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu,...
    yaw90_motor_forces), tspan, ICs_1pt4, options);

% Verification Plots for Trim State at 5m/s East, psi=90deg 
    N= length(t); % time steps from ode45 integration
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= yaw90_fc_4x1*ones(1,N);

% Defining plot limits manually based on observation
    xE_lims= [0, 0.1]; % [m]
    zE_lims= [1, 1.2]; % [m]
    v_lims= [0, 0.1]; % m/s
    w_lims= [-0.19, -.18]; % [m/s]

% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x', fc_4xN, figs.Q1pt4_yaw90, nonlin_col, ...
        xE_ylim=xE_lims, zE_ylim=zE_lims, v_ylim=v_lims, w_ylim=w_lims, ...
        xE_3D_xlim=xE_lims, zE_3D_zlim=zE_lims, ...
        figLabel= 'Task 1, Problem 4, \psi=90^{\circ}: ')

%% Task 1, Question 5: Stability of Steady, Hovering Flight
% Analyzing stability of simulation by testing various non-zero ICs. The
% following calculations assumes that the same initial conditions (and 
% motor forces) are used as in Q1.2-1.3, except for a disturbance to 
% analyze the natural response of% the minidrone.

% Control forces for Steady, Level, Hovering flight
    [fc_4x1, motor_forces] = openLoopForces(m, g, d, km);

% Initial conditions with small angular disturbance
    deltaTheta_deg= 2; % [deg]
    deltaTheta= deltaTheta_deg*(pi/180); % [rad]

% Initial conditions
    IC= getICs(z=-1.1, theta=deltaTheta);

[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu,...
    motor_forces), tspan, IC, options);

% Verification Plots for Hovering Trim State
    N= length(t); % time steps from ode45 integration
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);

% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x', fc_4xN, figs.Q1pt5, nonlin_col, figLabel= ...
    ['Task 1, Problem 5, \Delta\theta=' num2str(deltaTheta_deg) '^{\circ}: '])

%% Plotting Minidrone data
% Analyzing stability using Minidrone data
    droneDataFile= 'RSdata_nocontrol.mat';
% Extracting state matrix and ref. commands from data
    [x, t, refMat] = getDataState(droneDataFile) ;
% Plotting translational, rotational, & 3D plots with Minidrone data
    plotMiniDrone(x, t, refMat, figs.Q1pt5_data, figLabel= 'Task 1, Problem 5: ')





