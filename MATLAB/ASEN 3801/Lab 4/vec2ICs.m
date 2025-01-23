%% ASEN 3801, Creating an Initial Condition State Vector
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/25/2024
% Last Modified: 3/25/2024

% The vec2ICs (vector-to-ICs) function can take up to 4 inputs. The 4 
% inputs correspond to the four 3x1 vectors that are contained within the
% state vector (inertial position, 3-2-1 euler angles, body frame 
% velocities, body frame rotations). When no inputs are given, the initial 
% state vector is a 12x1 containing all zeros.

% SYNTAX
%       Unlike typical function inputs, the use of name-value arguments
%       requires values to be given to the structure field within the
%       function call. Field values should be assigned without the use of
%       dot notation (x=value instead of initialState.x=value).
%
%   Example: IC = getICs(P_E=[0;0;-5]); assigns values for the position
%       vector, while all other state vector entries default to zero.

% INPUTS (optional)
%       P_E: 3x1 vector of the form [x; y; z] (meters)
%       O: 3x1 vector of the form [phi; theta; psi] (radians)
%       V_B: 3x1 vector of the form [u; v; w] (meters/second)
%       omega_B: 3x1 vector of the form [p; q; r] (radians/second)
%
% OUTPUTS
%       IC: 12x1 vector containing the initial aircraft state

%% Begin Function Definition
function IC = vec2ICs(initialState)
% Defining default initial conditions for state vector
    arguments
    % Initial positions in inertial coordinates, seen by inertial frame
        initialState.P_E (3,1) double = [0; 0; 0] % [m] vector form
    
    % Initial 3-2-1 Euler angles phi, theta, psi
        initialState.O (3,1) double = [0; 0; 0] % [rad] vector form

    % Initial velocities in body coordinates, seen by inertial frame
        initialState.V_B (3,1) double = [0; 0; 0] % [m/s] vector form

    % Initial rotation rates in body coordinates, seen by inertial frame
        initialState.omega_B (3,1) double = [0; 0; 0] % [rad/s] vector form
    end

% Concatenating into state vector form
    IC= [initialState.P_E; initialState.O; initialState.V_B; ...
        initialState.omega_B];
end


