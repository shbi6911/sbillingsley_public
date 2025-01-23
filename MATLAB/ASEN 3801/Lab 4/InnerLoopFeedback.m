% Create a function to calculate the control vectors Fc and Gc. The function takes as input
% the 12x1 aircraft state var. The control force in the body z-direction should still equal the
% weight of the quadrotor (hard code the values in the function along with the control
% gains). Set the control moment using the control laws from Problem 3.1. 

function [Fc, Gc] = InnerLoopFeedback(var)
    % Dimensional variable
        N= width(var); % number of columns in var
    % Disturbance angles & rates from state vector
        phi= var(4,:); % [rad] roll
        theta= var(5,:); % [rad] pitch
        p= var(10,:); % [rad/s] roll rate
        q= var(11,:); % [rad/s] pitch rate
        r= var(12,:); % [rad/s] spin rate

    % Getting gravitational acceleration, minidronw inertia matrix & mass
        [m, g, ~, ~, ~, ~, ~]= getParams(); % [kg*m^2]

    % Control force (thrust) vector in body coordinates
        Fc= [0; 0; -m*g]*ones(1,N); % [N] Zc equal and opposite to weight
    
    % Control gains for real eigenvalues & time constant tau=0.5s
        K_phi= 0.0116; % [s^-2] proportional lateral gain
        K_p= 0.005916; % [s^-1] derivative lateral gain
        K_theta= 0.0072; % [s^-2] proportional longitudinal gain
        K_q= 0.003744; % [s^-1] derivative longitudinal gain

    % Roll control moment
        Lc= -K_phi*phi - K_p*p; % [N*m]
    % Pitch control moment
        Mc= -K_theta*theta - K_q*q; % [N*m]
    % Yaw control moment
        Nc= -0.004*r; % [N*m]
    % Vertically concatenating for final 3x1 control moment vector
        Gc= [Lc; Mc; Nc]; % [N*m]
end


