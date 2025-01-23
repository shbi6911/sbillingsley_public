% ASEN 3801 LAB 5
% Group 11

clear all
close all
clc

%% Problem 2
t_span_2 = [0 300];
aircraft_parameters = ttwistor;
% % 2.1 
% aircraft_state_init_2_1 = [0; 0; -1609.34; 0; 0; 0; 21; 0; 0; 0; 0; 0];
% aircraft_surfaces_2_1= [0; 0; 0; 0];
% wind_inertial_2_1 = [0; 0; 0];
% [tout_2_1, Xout_2_1] = ode45(@(time,aircraft_state) AircraftEOM(time, aircraft_state, aircraft_surfaces_2_1, wind_inertial_2_1, aircraft_parameters), t_span_2, aircraft_state_init_2_1) ;
% n_2_1 = length(tout_2_1);
% mat_2_1 = ones(4,n_2_1);
% control_input_array_2_1 = mat_2_1 .* [aircraft_surfaces_2_1(1); aircraft_surfaces_2_1(2); aircraft_surfaces_2_1(3); aircraft_surfaces_2_1(4)];
% fig = [figure(),figure(),figure(),figure(),figure(),figure()];
% PlotAircraftSim(tout_2_1, Xout_2_1, control_input_array_2_1, fig, "-b")

% %2.2 
% aircraft_state_init_2_2 = [ 0; 0; -1800; 0; 0.02780; 0; 20.99; 0; 0.5837; 0; 0; 0];
% aircraft_surfaces_2_2= [0.1079; 0; 0; 0.3182];
% wind_inertial_2_2 = [0; 0; 0];
% [tout_2_2, Xout_2_2] = ode45(@(time,aircraft_state) AircraftEOM(time, aircraft_state, aircraft_surfaces_2_2, wind_inertial_2_2, aircraft_parameters), t_span_2, aircraft_state_init_2_2) ;
% n_2_2 = length(tout_2_2);
% mat_2_2 = ones(4,n_2_2);
% control_input_array_2_2 = mat_2_2 .* [aircraft_surfaces_2_2(1); aircraft_surfaces_2_2(2); aircraft_surfaces_2_2(3); aircraft_surfaces_2_2(4)];
% fig = [figure(),figure(),figure(),figure(),figure(),figure()];
% PlotAircraftSim(tout_2_2, Xout_2_2, control_input_array_2_2, fig, "-b")

%2.3 
aircraft_state_init_2_3 =  [ 0; 0; -1800; deg2rad(15); deg2rad(-12); deg2rad(270); 19; 3; -2; deg2rad(0.08); deg2rad(-0.2); 0];
aircraft_surfaces_2_3= [deg2rad(5); deg2rad(2); deg2rad(-13); deg2rad(0.3)];
wind_inertial_2_3 = [0; 0; 0];
[tout_2_3, Xout_2_3] = ode45(@(time,aircraft_state) AircraftEOM(time, aircraft_state, aircraft_surfaces_2_3, wind_inertial_2_3, aircraft_parameters), t_span_2, aircraft_state_init_2_3) ;
n_2_3 = length(tout_2_3);
mat_2_3 = ones(4,n_2_3);
control_input_array_2_3 = mat_2_3 .* [aircraft_surfaces_2_3(1); aircraft_surfaces_2_3(2); aircraft_surfaces_2_3(3); aircraft_surfaces_2_3(4)];
fig = [figure(),figure(),figure(),figure(),figure(),figure()];
PlotAircraftSim(tout_2_3, Xout_2_3, control_input_array_2_3, fig, "-b")

%% Problem 3
t_span_3_1 = [0 3];
t_span_3_2 = [0 100];
aircraft_parameters = ttwistor;

% %3.1
% aircraft_state_init_3_1 = [ 0; 0; -1800; 0; 0.02780; 0; 20.99; 0; 0.5837; 0; 0; 0];
% aircraft_surfaces_3_1= [0.1079; 0; 0; 0.3182];
% wind_inertial_3_1 = [0; 0; 0];
% doublet_size = 15;
% doublet_time = 0.25;
% [tout_3_1, Xout_3_1] = ode45(@(time,aircraft_state) AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces_3_1, doublet_size,doublet_time,wind_inertial_3_1, aircraft_parameters), t_span_3_1, aircraft_state_init_3_1) ;
% n_3_1 = length(tout_3_1);
% mat_3_1 = ones(4,n_3_1);
% control_input_array_3_1 = mat_3_1 .* [aircraft_surfaces_3_1(1); aircraft_surfaces_3_1(2); aircraft_surfaces_3_1(3); aircraft_surfaces_3_1(4)];
% fig = [figure(),figure(),figure(),figure(),figure(),figure()];
% PlotAircraftSim(tout_3_1, Xout_3_1, control_input_array_3_1, fig, "-b")

%3.2
aircraft_state_init_3_2 = [ 0; 0; -1800; 0; 0.02780; 0; 20.99; 0; 0.5837; 0; 0; 0];
aircraft_surfaces_3_2= [0.1079; 0; 0; 0.3182];
wind_inertial_3_2 = [0; 0; 0];
doublet_size = 15;
doublet_time = 0.25;
[tout_3_2, Xout_3_2] = ode45(@(time,aircraft_state) AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces_3_2, doublet_size,doublet_time,wind_inertial_3_2, aircraft_parameters), t_span_3_2, aircraft_state_init_3_2) ;
n_3_2 = length(tout_3_2);
mat_3_2 = ones(4,n_3_2);
control_input_array_3_2 = mat_3_2 .* [aircraft_surfaces_3_2(1); aircraft_surfaces_3_2(2); aircraft_surfaces_3_2(3); aircraft_surfaces_3_2(4)];
fig = [figure(),figure(),figure(),figure(),figure(),figure()];
PlotAircraftSim(tout_3_2, Xout_3_2, control_input_array_3_2, fig, "-b")