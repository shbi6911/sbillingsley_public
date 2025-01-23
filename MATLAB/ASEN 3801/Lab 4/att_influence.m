%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/24/2024
% Last Modified: 3/24/2024

% The att_influence function calculates and outputs the 3x3 attitude 
% influence matrix when given a 3x1 vector of Euler angles in the order 
% roll (phi), pitch (theta), yaw (psi).

% INPUTS: 3x1 vector [roll (phi); pitch (theta); yaw (psi)] in radians

% OUTPUTS: 3x3 attitude influence matrix

%% Begin Function Definition
function T = att_influence(eulerVec)
    % Extracting angles from input vector phiThetaPsi
        phi= eulerVec(1); % [rad] roll angle phi
        theta= eulerVec(2); % [rad] pitch angle theta
    
    % Attitude influence matrix T
        T= [ 1 sin(phi)*tan(theta) cos(phi)*tan(theta);
            0 cos(phi) -sin(phi);
            0 sin(phi)*sec(theta) cos(phi)*sec(theta)] ;
end