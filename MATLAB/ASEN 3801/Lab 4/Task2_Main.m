%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/26/2024
% Last Modified: 4/14/2024

clear; close all; clc % housekeeping

% Deviations to be plotted
%   a) Deviation by +5 deg in roll
%   b) Deviation by +5 deg in pitch
%   c) Deviation by +5 deg in yaw
%   d) Deviation by +0.1 rad/sec in roll rate
%   e) Deviation by +0.1 rad/sec in pitch rate
%   f) Deviation by +0.1 rad/sec in yaw rate

%% Begin Task 2
% Getting minidrone parameters
    [m, g, d, km, I, nu, mu]= getParams(); % all in SI units
% Open Loop control forces (constant for all of Problem 1)
    [fc_4x1, motor_forces] = openLoopForces(m, g, d, km);
% Time to integrate over (constant throughout all of Task 1)
    tspan= [0, 10]; % [s]
% Tolerance settings for ode45
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% Reference height
    zE_ref= -10; % (m)
% For ideal system, control disturbances are all zero (thrust output is
% exactly equal to desired thrust input)
    delta_Fc= zeros(3,1); % [N]
    delta_Gc= delta_Fc; % [N*m]

% Figure variables/settings
    % Struct containing vectors of figure numbers
        figs.T2_Q1a=34:39; % problem 1a
        figs.T2_Q1b=40:45; % problem 1b
        figs.T2_Q1c=46:51; % problem 1c
        figs.T2_Q1d=52:57; % problem 1d
        figs.T2_Q1e=58:63; % problem 1e
        figs.T2_Q1f=64:69; % problem 1f
        figs.T2_Q5d=70:75; % problem 5, repeating 2.1d deviation
        figs.T2_Q5e=76:81; % problem 5, repeating 2.1e deviation
        figs.T2_Q5f=82:87; % problem 5, repeating 2.1f deviation

    % Line colors
        nonlin_col= "#0072BD"; % nonlinear EOM
        lin_col= "#D95319"; % linearized EOM
        rateFeedback_col= "#77AC30"; % nonlinear EOM w/ rate feedback

% Cell arrays containing legend entries
% Parts 1 & 2
    legendLabels= {'Nonlinear EOM', 'Linear EOM' };
    legendLabels_3D= {'Nonlinear EOM', '', '', 'Linear EOM', ...
        'Initial Position', 'Final Position'};
% Part 5
    legendLabels_Part5= {'Nonlinear EOM', 'Nonlinear EOM w/Rate Feedback' };
    legendLabels_3D_Part5= {'Nonlinear EOM', '', '', 'Nonlinear EOM w/Rate Feedback', ...
        'Initial Position', 'Final Position'};

%% Task 2, Problem 1a: 5 Degree Deviation in Roll
% Given disturbance
    deltaPhi= 5*pi/180; % [rad]
% Initial conditions
    IC = getICs(z=zE_ref, phi=deltaPhi);

% NONLINEAR ode45 call with aerodynamic forces & 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1a, nonlin_col)

% LINEAR ode45 call with aerodynamic forces & 5 deg. roll disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1a, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D= legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1a & 2a: ')

%% Task 2, Problem 1b: 5 Degree Deviation in Pitch
% Given disturbance
    deltaTheta= 5*pi/180; % [rad]
% Initial conditions
    IC = getICs(z=zE_ref, theta=deltaTheta);

% NONLINEAR ode45 call with aerodynamic forces & 5 deg. pitch disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1b, nonlin_col)

% LINEAR ode45 call with aerodynamic forces & 5 deg. pitch disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1b, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D=legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1b & 2b: ')

%% Task 2, Problem 1c: 5 Degree Deviation in Yaw
% Given disturbance
    deltaPsi= 5*pi/180; % [rad]
% Initial conditions
    IC = getICs(z=zE_ref, psi=deltaPsi);

% NONLINEAR ode45 call with aerodynamic forces & 5 deg. yaw disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1c, nonlin_col)

% LINEAR ode45 call with aerodynamic forces & 5 deg. yaw disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1c, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D= legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1c & 2c: ')

%% Task 2, Problem 1d: 0.1 rad/sec Deviation in Roll Rate
% Given disturbance
    delta_p= 0.1; % [rad/s]
% Initial conditions
    IC = getICs(z=zE_ref, p=delta_p);

% NONLINEAR ode45 call with aerodynamic forces & 0.1 rad/sec roll rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1d, nonlin_col) % Problem 1d
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5d, nonlin_col) % Problem 5d

% LINEAR ode45 call with aerodynamic forces & 0.1 rad/sec roll rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1d, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D=legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1d & 2d: ')

% NONLINEAR ode45 call with rate feedback, aerodynamic forces, and 
% 0.1 rad/sec roll rate disturbance
[t, x]= ode45( @(t, var) QuadrotorEOMwithRateFeedback(t, var, g, m, I, nu, ...
    mu), tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = RotationDerivativeFeedback(x', m, g);
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5d, rateFeedback_col, ...
    LegendEntries=legendLabels_Part5, LegendEntries_3D=legendLabels_3D_Part5, ...
    LineStyle= '-.', figLabel= 'Task 2, Problem 5d: ')

%% Task 2, Problem 1e: 0.1 rad/sec Deviation in Pitch Rate
% Given disturbance
    delta_q= 0.1; % [rad/s]
% Initial conditions
    IC = getICs(z=zE_ref, q=delta_q);

% NONLINEAR ode45 call with aerodynamic forces & 0.1 rad/sec pitch rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1e, nonlin_col) % Problem 1e
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5e, nonlin_col) % Problem 5e

% LINEAR ode45 call with aerodynamic forces & 0.1 rad/sec pitch rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1e, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D=legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1e & 2e: ')

% NONLINEAR ode45 call with rate feedback, aerodynamic forces, and
% 0.1 rad/sec pitch rate disturbance
[t, x]= ode45( @(t, var) QuadrotorEOMwithRateFeedback(t, var, g, m, I, nu, ...
    mu), tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = RotationDerivativeFeedback(x', m, g);
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5e, rateFeedback_col, ...
    LegendEntries=legendLabels_Part5, LegendEntries_3D=legendLabels_3D_Part5, ...
    LineStyle= '-.', figLabel= 'Task 2, Problem 5e: ')

%% Task 2, Problem 1f: 0.1 rad/sec Deviation in Yaw Rate
% Given disturbance
    delta_r= 0.1; % [rad/s]
% Initial conditions
    IC = getICs(z=zE_ref, r=delta_r);

% NONLINEAR ode45 call with aerodynamic forces & 0.1 rad/sec yaw rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, ...
    motor_forces), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1f, nonlin_col)
    PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5f, nonlin_col)

% LINEAR ode45 call with aerodynamic forces & 0.1 rad/sec yaw rate disturbance
[t, x]= ode45(@(t, var) QuadrotorEOM_Linearized(t, var, g, m, I, delta_Fc, ...
    delta_Gc), tspan, IC, options);
% Number of time steps from ode45 integration
    N= length(t);
% Converting to 4xN via matrix multiplication (control forces are constant
% since there is no feedback)
    fc_4xN= fc_4x1*ones(1,N);
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q1f, lin_col, LegendEntries= ...
    legendLabels, LegendEntries_3D=legendLabels_3D, LineStyle= '-.', ...
    figLabel= 'Task 2, Problem 1f & 2f: ')

% NONLINEAR ode45 call with rate feedback, aerodynamic forces, and 
% 0.1 rad/sec yaw rate disturbance
[t, x]= ode45( @(t, var) QuadrotorEOMwithRateFeedback(t, var, g, m, I, nu, ...
    mu), tspan, IC, options);
% Calculating non-constant control forces
    [Fc, Gc] = RotationDerivativeFeedback(x', m, g);
    fc_4xN= [Fc(3,:); Gc];
% Plotting with PlotAircraftSim function
PlotAircraftSim(t, x',fc_4xN, figs.T2_Q5f, rateFeedback_col, ...
    LegendEntries=legendLabels_Part5, LegendEntries_3D=legendLabels_3D_Part5, ...
    LineStyle= '-.', figLabel= 'Task 2, Problem 5f: ')

%%printing plots
fig_start = figs.T2_Q1a(1);
fig_end = figs.T2_Q5f(end);
counter = 0;
for ii = fig_start:fig_end
    names = string(fieldnames(figs));
    probnum = floor(counter/6)+1;
    fignum = mod(counter,6)+1;
    filename = names(probnum) + "_" + string(fignum);
    disp(filename);
    saveas(figure(ii),'C:\Users\calda\MATLAB Drive\ASEN 3801\Lab 4\figures\' + filename,'png');
    counter = counter+1;
end