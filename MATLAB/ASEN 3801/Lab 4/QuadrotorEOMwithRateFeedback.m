% Add the feedback controller to your nonlinear quadrotor equations of motion. Create a
% new equation of motion function with the feedback controller from Problem 2.4 in it:
% function var_dot = QuadrotorEOMwithRateFeedback(t, var, g, m, I, nu, mu)
% Repeat simulations for the setups in Problems 2.1.d through 2.1.f. What is the effect of
% this control law? Support your argument with new plots for each setup, overlaying the
% controlled and uncontrolled nonlinear systems. Also calculate and plot the motor thrust
% forces for the controlled and uncontrolled systems.

function var_dot = QuadrotorEOMwithRateFeedback(t, var, g, m, I, nu, mu)
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
% Control force and control moment vectors (Zc, Lc, Mc, & Nc calculations)
    [Fc, Gc] = RotationDerivativeFeedback(var, m, g);

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
        % External (aerodynamic) force vector
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


