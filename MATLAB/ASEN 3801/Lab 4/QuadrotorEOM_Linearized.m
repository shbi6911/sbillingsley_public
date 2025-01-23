% Function: QuadrotorEOM_Linearized()
% Use: Linearized Steady State Hover Dynamics
% Changelog:
%   3/20/2024 -- Adrian Bryant -- Created/Commented
%   4/14/2024 -- Kyle Goodall
% Inputs:
%   t --> time [s]
%   var --> 12x1 state matrix:
%       [delta_x_e ; delta_y_e; delta_z_e] 3x1 of inertial XYZ coordinate deviations [m]
%       [delta_phi; delta_theta; delta_psi] 3x1 of deviations of Euler Angles from inertial->body frame
%       [rad/s]
%           (MUST BE RADIANS) 
%       [delta_u_e ; delta_v_e ; delta_w_e] 3x1 Body velocity deviations (inertially measured) [m/s]
%       [delta_p ; delta_q ; delta_r] 3x1 Body Rotation rate deviations [rad/s]
%   g --> gravity [m/s^2]
%   m --> mass [kg]
%   I --> 3x3 Moment of Inertia matrix [kg*m^2]
%   deltaFc -->  Vector of Deviation of Control Forces of the form:
%       [0; 0; delta_Zc]
%       [N]
%   deltaGc -->  Vector of Deviation of Control Moments of the form:
%       [delta_Lc ; delta_Mc ; delta_Nc]
%       [N-m]
% Outputs:
%   var_dot --> 12x1 state matrix derivative:
%       [delta_x_e_dot ; delta_y_e_dot ; delta_z_e_dot ] 3x1 of derivative of intertial XYZ coordinate deviations [m/s]
%       [delta_phi_dot; delta_theta_dot; delta_psi_dot] 3x1 of Euler Angle Derivatives from inertial->body frame
%       [rad/s^2]
%       [delta_u_e_dot ; delta_v_e_dot ; delta_w_e_dot] 3x1 Derivatives of Body velocities (inertially measured) [m/s^2]
%       [delta_p_dot ; delta_q_dot ; delta_r_dot] 3x1 Derivatives of Body Rotation rates [rad/s^2]

function var_dot = QuadrotorEOM_Linearized(t, var, g, m, I, deltaFc, deltaGc)    
    % Extracting disturbance roll and pitch Euler angles
        deltaPhi= var(4); % [rad] roll
        deltaTheta = var(5); % [rad] pitch
    % Extracting disturbances in body frame velocities, [u_B, v_B, w_B]
        deltaV_B = var(7:9); % [m/s]
    % Extracting disturbances in angular rates, [p; q; r]
        deltaOmega = var(10:12);

    % Derivative of inertial position disturbances
        deltaP_E = deltaV_B; % [m/s]
    % Derivative of 3-2-1 Euler angle disturbances
        deltaEulerDot = deltaOmega; % [rad/s]
    % Derivative of body frame velocity disturbances
        deltaV_B= g*[-deltaTheta; deltaPhi; 0] + (1/m)*deltaFc; % [m/s^2]
    % Derivative of angular rate disturbances
        deltaOmegaDot = I*(deltaGc); % [rad/s]

    % Return all derivatives in the output
        var_dot = [deltaP_E; deltaEulerDot; deltaV_B; deltaOmegaDot];
end

