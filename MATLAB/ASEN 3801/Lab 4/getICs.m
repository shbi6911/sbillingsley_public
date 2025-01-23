%% ASEN 3801, Creating an Initial Condition State Vector
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/25/2024
% Last Modified: 3/25/2024

% The getICs function can take up to 12 inputs. The inputs correspond
% to the 12 scalar entries within the 12x1 state vector. When no inputs are
% given, the initial state vector is a 12x1 containing all zeros.

% SYNTAX
%       Unlike typical function inputs, the use of name-value arguments
%       requires values to be given to the structure field within the
%       function call. Field values should be assigned without the use of
%       dot notation (x=value instead of initialState.x=value).
%
%   Example: IC = getICs(x=3, z=-5); assigns values for the initial x and z
%       positions, while all other state vector entries default to zero.
%   Alternate syntax: val1=3; val2=-5; IC=getICs(x=val1, z=val2);

% INPUTS (optional)
%       x: initial x-position in inertial coordinates (meters)
%       y: initial y-position in inertial coordinates (meters)
%       z: initial z-position in inertial coordinates (meters)
%       phi: initial roll angle (radians)
%       theta: initial pitch angle (radians)
%       psi: initial yaw angle (radians)
%       u: initial x-velocity in body coordinates (meters/second)
%       v: initial y-velocity in body coordinates (meters/second)
%       w: initial z-velocity in body coordinates (meters/second)
%       p: initial roll rate in body coordinates (radians/second)
%       q: initial pitch rate in body coordinates (radians/second)
%       r: initial yaw rate in body coordinates (radians/second)
%
% OUTPUTS
%       IC: 12x1 vector containing the initial aircraft state

%% Begin Function Definition
function IC = getICs(initialState)
% Defining default initial conditions for state vector
    arguments
    % Initial positions in inertial coordinates, seen by inertial frame
        initialState.x (1,1) double = 0 % [m] x-position
        initialState.y (1,1) double = 0 % [m] y-position
        initialState.z (1,1) double = 0 % [m] z-position
  
    % Initial 3-2-1 Euler angles phi, theta, psi
        initialState.phi (1,1) double = 0 % [rad] roll
        initialState.theta (1,1) double = 0 % [rad] pitch 
        initialState.psi (1,1) double = 0 % [rad] yaw

    % Initial velocities in body coordinates, seen by inertial frame
        initialState.u (1,1) double = 0 % [m/s] x-velocity
        initialState.v (1,1) double = 0 % [m/s] y-velocity
        initialState.w (1,1) double = 0 % [m/s] z-velocity

    % Initial rotation rates in body coordinates, seen by inertial frame
        initialState.p (1,1) double = 0 % [rad/s] about body x-axis
        initialState.q (1,1) double = 0 % [rad/s] about body y-axis
        initialState.r (1,1) double = 0 % [rad/s] about body z-axis
    end

% Subvariables for 3x1 vectors
    P_E= [initialState.x ; initialState.y; initialState.z]; % [m]
    O= [initialState.phi ; initialState.theta; initialState.psi]; % [rad]
    V_B= [initialState.u ; initialState.v; initialState.w]; % [m/s]
    omega_B= [initialState.p ; initialState.q; initialState.r]; % [rad/s]
% Concatenating into state vector form
    IC = [P_E; O; V_B; omega_B];
end


