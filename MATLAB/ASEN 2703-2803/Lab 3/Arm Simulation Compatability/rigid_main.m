%Bobby Hodgkinson
%3/31/2020
%Load parameters for simulink model to matlab workspace.

%bobby.Figures = [];
if exist('bobby')
close(bobby.Figures)
end
clc; clear all; close all;


beep off %turns off the annoying beep from simulink completion

%% User inputs - control gains
Kpt = 10.500;   %Proportional gain theta
Kdt = 0.5;    %Derivative gain theta

Kpd = 0;    %Proportional gain displacement (tip)
Kdd = 0;    %Derivative gain displacement (tip)

%% User inputs - simulation parameters
Sim_time = 10;       %[Seconds] Simulation run time
saturation = 10;    %[Volts] Motor saturation voltage. Default = 10
Amplitude = 0.5;    %[Rads] Amplitude of hub reference step command
Frequency = 0.2;    %[Hz] Frequency of hub reference step command
dead_zone = .7;   %[Volts] Motor dead zone voltage. Default = 0.25

%% run the parameters.m script to load the parameters to the workspace for simulink
evalin('base','parameters') 

%% Load the simulink and Run it

handle = load_system('rigid_comp');     %Load the rigid simulink simulation. The handle argument prevents the model from opening
open_system('rigid_comp/Hub Angle');    %Opens the Hub Angle scope
open_system('rigid_comp/Motor V');      %Opens the Motor Voltage scope

myworld = vrworld('rigid_2020.wrl');
view(myworld);
%melvin = view(myworld);
%set(melvin,'azimuth',30,'elevation',90)
bobby = get(myworld);

data = sim('rigid_comp');               %This will save the data from the workspace to 'data'
%% Comment the following line to disable autoclose of the VR Animation
%close(bobby.Figures);

%Time (ms)	Hubangle(Theta in rad)	Tip Deflection(m)	Hub Angular Velocity (rads/s)	Tip Velocity (m/s)	Position Reference (rad)	Output Voltage (V)	K1 (Hub Prop)	K2 (Tip Prop)	K3 (Hub Deriv)	K4 (Tip Deriv)
