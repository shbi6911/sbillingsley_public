% Using the function from Problem 3.2, create a new EOM function to simulate the
% response of the closed loop linearized system.

function var_dot = QuadrotorEOM_Linear_CL(t, var, g, m, I)
    % Extracting disturbance roll and pitch Euler angles
        deltaPhi= var(4); % [rad] roll
        deltaTheta = var(5); % [rad] pitch
    % Extracting disturbances in body frame velocities, [u_B, v_B, w_B]
        deltaV_B = var(7:9); % [m/s]
    % Extracting disturbances in angular rates, [p; q; r]
        deltaOmega = var(10:12);

    % Steady hover controls at trim
        Fc0= [0; 0; -m*g]; % [N] trim control forces
        Gc0= [0; 0; 0]; % [N*m] trim control moments
    % Calculating control vectors for current iteration
        [Fc, Gc] = InnerLoopFeedback(var);
        deltaFc= Fc - Fc0; % [N] disturbance control thrust
        deltaGc= Gc - Gc0; % [N*m] disturbance control moments
    % Derivative of inertial position disturbances
        deltaP_E = deltaV_B; % [m/s]
    % Derivative of 3-2-1 Euler angle disturbances
        deltaEulerDot = deltaOmega; % [rad/s]
    % Derivative of body frame velocity disturbances
        deltaV_B= g*[-deltaTheta; deltaPhi; 0] + (1/m)*deltaFc; % [m/s^2]
    % Derivative of angular rate disturbances
        deltaOmegaDot = I*(deltaGc); % [rad/s]

    % Vertical concatenation to form derivative of state vector
        var_dot = [deltaP_E; deltaEulerDot; deltaV_B; deltaOmegaDot];
end