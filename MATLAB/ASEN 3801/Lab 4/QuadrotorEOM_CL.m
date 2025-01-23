% Repeat Problem 3.3 (i.e. create a new EOM function that adds the control to the nonlinear
% quadrotor dynamics) using the nonlinear dynamics model together with the feedback
% control design for the linearized system, and compare the closed loop linearized and
% nonlinear behaviors. 


function var_dot= QuadrotorEOM_CL(t, var, g, m, I, nu, mu)

% Extracting 3x1 vectors from current state
    % 321 Euler Angles phi (roll), theta (pitch), and psi (yaw)
        euler = var(4:6); % [rad]
    % Body coordinate velocities seen by inertial frame, V_B
        V_B= var(7:9); % [m/s]
    % Angular rates in body coordinates seen by inertial frame, w_B
        w_B= var(10:12); % [rad/s]

% Attitude influence matrix
    T= att_influence(euler);

% Control force & moment vectors with roll & pitch feedback (Zc,Lc,Mc,& Nc)
    [Fc, Gc] = InnerLoopFeedback(var); % [N, N*m]
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