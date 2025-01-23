%Bobby Hodgkinson
%3/31/2020
%Load parameters for simulink model to matlab workspace.


clc; clear all;

beep off %turns off the annoying beep from simulink completion

%% User inputs - control gains
Kpt = 15;           %Proportional gain theta
Kdt = 0;            %Derivative gain theta
Kpd = 0;            %Proportioanl gain displacement (tip)
Kdd = 0;            %Derivative gain displacement (tip)

%% User inputs - simulation parameters
Sim_Time = 5;       %[Seconds] Simulation run time
saturation = 10;    %[Volts] Motor saturation voltage. Default = 10
Amplitude = 0.2;    %[Rads] Amplitude of hub reference step command
Frequency = 0.1;    %[Hz] Frequency of hub reference step command
dead_zone = 0.25;   %[Volts] Motor dead zone voltage. Default = 0.25


%% run the parameters.m script to load the parameters to the workspace for simulink
evalin('base','parameters') 

%% Load the simulink and Run it

handle = load_system('flexible_comp'); %the handle argument prevents the model from opening
open_system('flexible_comp/Hub Angle');
open_system('flexible_comp/Tip');
open_system('flexible_comp/Motor V');

myworld = vrworld('flexible_2020.wrl');
view(myworld);
bobby = get(myworld);

data = sim('flexible_comp'); %This will save the data from the workspace to 'data'
%% Comment the following line to disable autoclose of the VR Animation
%close(bobby.Figures);

%Time (ms)	Hubangle(Theta in rad)	Tip Deflection(m)	Hub Angular Velocity (rads/s)	Tip Velocity (m/s)	Position Reference (rad)	Output Voltage (V)	K1 (Hub Prop)	K2 (Tip Prop)	K3 (Hub Deriv)	K4 (Tip Deriv)
