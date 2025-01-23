% Function: QuardotorEOM()
% Changelog:
%   3/13/2024 -- Adrian Bryant -- Created/Commented
%   4/13/2024 -- Kyle Goodall
% Inputs:
%   t --> time [s]
%   var --> 12x1 state matrix:
%       [x_e ; y_e; z_e] 3x1 of intertial XYZ coordinates [m]
%       [phi; theta; psi] 3x1 of Euler Angles from inertial->body frame
%       [rad/s]
%           (MUST BE RADIANS) 
%       [u_e ; v_e ; w_e] 3x1 Body velocities (inertially measured) [m/s]
%       [p ; q ; r] 3x1 Body Rotation rates [rad/s]
%   g --> gravity [m/s^2]
%   m --> mass [kg]
%   I --> 3x3 Moment of Inertia matrix [kg*m^2]
%   d --> distance from rotors to CM [m]
%   km --> Control Moment Coefficient [N*m/N]
%   nu --> Aerodynamic Force Coefficient [N/(m/s)^2]
%   mu --> Aerodynamic Moment Coefficient [Nm/(rad/s)^2]
%   motor_forces --> 4x1 force matrix
% Outputs:
%   var_dot --> 12x1 state matrix derivative:
%       [x_e_dot ; y_e_dot ; z_e_dot ] 3x1 of derivative of intertial XYZ coordinates [m/s]
%       [phi_dot; theta_dot; psi_dot] 3x1 of Euler Angle Derivatives from inertial->body frame
%       [rad/s^2]
%       [u_e_dot ; v_e_dot ; w_e_dot] 3x1 Derivatives of Body velocities (inertially measured) [m/s^2]
%       [p_dot ; q_dot ; r_dot] 3x1 Derivatives of Body Rotation rates [rad/s^2]
function var_dot= QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, motor_forces)

% Extracting 3x1 vectors from current state
    % % Inertial position seen by inertial frame, P_E
    %     P_E= var(1:3); % [m]
    % 321 Euler Angles phi (roll), theta (pitch), and psi (yaw)
        euler = var(4:6); % [rad]
    % Body coordinate velocities seen by inertial frame, V_B
        V_B= var(7:9); % [m/s]
    % Angular rates in body coordinates seen by inertial frame, w_B
        w_B= var(10:12); % [rad/s]

% Attitude influence matrix
    T= att_influence(euler);

% Calculating Zc, Lc, Mc, Nc with matrix equation
coeff_mat= [     -1              -1            -1            -1      ; ...
            -(d/sqrt(2))    -(d/sqrt(2))   (d/sqrt(2))   (d/sqrt(2)) ; ...
             (d/sqrt(2))    -(d/sqrt(2))  -(d/sqrt(2))   (d/sqrt(2)) ; ...
                  km             -km            km          -km     ];
% Control force & moments vector in body coordinates [Zc; Lc; Mc; Nc]
    fc_4x1= coeff_mat*motor_forces ; % [N and N*m]

% Extracting 3x1 force and moment control vectors
    Fc= [0; 0; fc_4x1(1)]; % [N] thrust force
    Gc= fc_4x1(2:4); % [N*m] control moments

% Position derivatives
    Pdot_E= (R_E2B(euler)')*V_B; % [m/s^2]
% Euler angle derivatives
    eulerDot= T*w_B; % [rad/s^2]

% Components of velocity derivative vector
    % First term: accelerations due to rotation
        a_rot= cross(w_B, V_B); % [m/s^2]
    % Second term: acceleration due to gravity
        a_grav= R_E2B(euler)*[0; 0; g]; % [m/s^2]
    % Third term: accelerations due to external (aerodynamic) forces
        % Air-relative velocity from current state (ASSUMES NO WIND)
            Va= norm(V_B); % [m/s]
        % External (aerodynamic) force & acceleration vector
            f_ext= -nu*Va*V_B; % [N]
            a_ext= (1/m)*f_ext; % [m/s^2]
    % Fourth term: acceleration due to control force (thrust)
        a_c= (1/m)*Fc; % [m/s^2]
% Velocity derivative vector
    Vdot_B= a_rot + a_grav + a_ext + a_c ; % [m/s^2]

% Components of angular rate derivative
    % External (aerodynamic) moments
        G_ext= -mu*norm(w_B)*w_B; % [N*m]
    % Sum of moments due to external forces (control+aerodynamic)
        G= Gc + G_ext; % [N*m]
    % "Rotational" term from kinetic transport theorem
        G_rot= cross(w_B,I*w_B); % [N*m]
% Angular rate derivative
    wDot_B= (I^-1)*(G - G_rot); % [rad/s^2]

% Vertically concatenating for final 12x1 state vector derivative
    var_dot= [Pdot_E; eulerDot; Vdot_B; wDot_B];

end