%% ASEN 3801, Lab 4: Function for Calculating Level Controls Forces
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/26/2024
% Last Modified: 4/14/2024

% The openLoopCtrls function is used to shorten preliminary calculations
% required before calling ode45 and the PlotAircraftSim functions.

% INPUTS
%       m: minidrone mass (kilograms)
%       g: acceleration of gravity ( meters/(second^2) )
%       d: distance between motors and center of mass (meters)
%       km: control moment coefficient ( Newton*meters/Newtons )
%       angle: struct containing options for non-level trim conditions
%       angle.tiltAngle: angle (radians) at which the drone has deviated 
%           from a level condition (assumes yaw=0deg and pitch or roll are 
%           not non-zero simultaneously).
%
% OUTPUTS
%       fc: 4x1 vector of the thrust forces in body coordinates (Newtons) 
%           provided by each motor (assumes symmetric thrust due to zero 
%           net moment required for trimmed, non-manuevering flight)
%       ctrl_in= 4x1 vector of the net thrust force Zc (Newtons) and the 
%           three control moments Lc, Mc, Nc (Netwon*meters), all in body 
%           coordinates

%% Begin Function Definition
function [fc_4x1, motor_forces] = openLoopForces(m, g, d, km, trim)
    arguments
        m (1,1) double
        g (1,1) double
        d (1,1) double
        km (1,1) double
        trim.angle (1,1) double = 0 % [rad] default orientation (level)
    end
% Control forces for Steady, Hovering flight
    fc_mag= m*g/cos(trim.angle); % [N] thrust opposes weight
    motor_forces = (fc_mag/4)*ones(4,1); % [N] symmetric thrust

% Calculating Zc, Lc, Mc, Nc with matrix equation
coeff_mat= [     -1              -1            -1            -1      ; ...
            -(d/sqrt(2))    -(d/sqrt(2))   (d/sqrt(2))   (d/sqrt(2)) ; ...
             (d/sqrt(2))    -(d/sqrt(2))  -(d/sqrt(2))   (d/sqrt(2)) ; ...
                  km             -km            km          -km     ];

% Control force & moments vector in body coordinates [Zc; Lc; Mc; Nc]
    fc_4x1= coeff_mat*motor_forces ; % [N and N*m]

end



