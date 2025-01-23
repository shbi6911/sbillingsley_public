% Create a function to calculate the control vectors Fc and Gc. The function takes as input
% the 12x1 aircraft state var, aircraft mass m, and gravitational acceleration g. The control
% force in the body z-direction should still equal the weight of the quadrotor. Set the control
% moments about each body axis proportional to the rotational rates about their respective
% axes, but in the opposite sign of the angular velocity with a gain of 0.004 Nm/(rad/sec): 

function [Fc, Gc] = RotationDerivativeFeedback(var, m, g)
    % Extracting rotation rates p, q, & r from state vector
    omega_B= var(10:12,:); % [rad/s]
    N= width(var); % number of columns in var
    Fc= [0; 0; -m*g]*ones(1,N); % [N] Zc equal and opposite to weight
    Gc= -0.004*omega_B; % [N*m]
end