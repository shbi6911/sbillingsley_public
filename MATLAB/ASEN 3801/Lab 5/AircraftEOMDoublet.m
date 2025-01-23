function xdot = AircraftEOMDoublet(time, aircraft_state, aircraft_surfaces, doublet_size,doublet_time,wind_inertial, aircraft_parameters)
% INPUTS:
    % time [s]
    % aircraft_state = 12 x 1 aircraft state vector
        % [xi; yi; zi; roll; pitch; yaw; u; v; w; p; q; r]
    % aircraft_surfaces = 4 x 1 control surface vector [de; da; dr; dt]
    % wind_inertial = 3 x 1 inertial wind velocity in inertial coordinates
    % aircraft_parameters struct
    % doublet_size = size of applied doublet in degrees
    % doublet_time = time to apply each half of doublet in seconds
% OUTPUT:
    % xdot = 12 x 1 derivative of state vector

    % Setup
    ap = aircraft_parameters;
    g = ap.g;

    %Apply doublet to elevator control
    delta_e_trim = aircraft_surfaces(1);
    if time <= doublet_time
        aircraft_surfaces(1) = delta_e_trim + deg2rad(doublet_size);
    elseif time <= 2*doublet_time
        aircraft_surfaces(1) = delta_e_trim - deg2rad(doublet_size);
    else
        aircraft_surfaces(1) = delta_e_trim;
    end

    % Get aero forces and moments
    density = stdatmo(-aircraft_state(3)); %pull z component of state vector 
    [aero_forces, aero_moments] = AeroForcesAndMoments(aircraft_state, ...
        aircraft_surfaces, wind_inertial, density, aircraft_parameters);
    X = aero_forces(1);
    Y = aero_forces(2);
    Z = aero_forces(3);
    L = aero_moments(1);
    M = aero_moments(2);
    N  = aero_moments(3);
    
    % Setup state vector
        % x_pos = x(1);
        % y_pos = x(2);
        % z_pos = x(3);
    pos = [aircraft_state(1); aircraft_state(2); aircraft_state(3)];
        phi = aircraft_state(4);
        theta = aircraft_state(5);
        Psi = aircraft_state(6);
    orientation = [aircraft_state(4); aircraft_state(5); aircraft_state(6)];
        u = aircraft_state(7);
        v = aircraft_state(8);
        w = aircraft_state(9);
    vel = [aircraft_state(7); aircraft_state(8); aircraft_state(9)];
        p = aircraft_state(10);
        q = aircraft_state(11);
        r = aircraft_state(12);
    rates = [aircraft_state(10); aircraft_state(11); aircraft_state(12)];
    
    % Position derivatives
        % pos_dot = [ cos(theta)*cos(Psi), sin(phi)*sin(theta)*cos(Psi) - cos(phi)*sin(Psi), cos(phi)*sin(theta)*cos(Psi) + sin(phi)*sin(Psi);
        %                  cos(theta)*sin(Psi), sin(phi)*sin(theta)*sin(Psi) + cos(phi)*cos(Psi), cos(phi)*sin(theta)*sin(Psi) - sin(phi)*cos(Psi);
        %                  -sin(theta),              sin(phi)*cos(theta),                                    cos(phi)*cos(theta) ] ...
        %            * pos;
    R_supB_subE = rotation321(orientation);
    R_supE_subB = R_supB_subE';
    pos_dot = R_supE_subB * vel;
    
    % Orientation derivatives
        phi_dot = p + (sin(phi) * tan(theta))*q + (cos(phi) * tan(theta))*r ;
        theta_dot = cos(phi)*q - sin(phi)*r ;
        psi_dot = (sin(phi) * sec(theta))*q + (cos(phi) * sec(theta))*r;
    orientation_dot = [phi_dot; theta_dot; psi_dot];
    
    % Velocity derivatives
        u_dot = r*v - q*w - g*sin(theta) + 1/ap.m * X;
        v_dot = p*w - r*u + g*cos(theta)*sin(phi) + 1/ap.m * Y;
        w_dot = q*u - p*v  + g*cos(theta)*cos(phi) +1/ap.m * Z;
    vel_dot = [u_dot; v_dot; w_dot];
    
    % Rate derivatives
    gamma = ap.Ix*ap.Iz - ap.Ixz^2;
    gamma1 = (ap.Ixz * (ap.Ix - ap.Iy + ap.Iz))/gamma; 
    gamma2 = (ap.Iz * (ap.Iz - ap.Iy) + ap.Ixz^2)/gamma;
    gamma3 = ap.Iz / gamma;
    gamma4 = ap.Ixz / gamma;
    gamma5 = (ap.Iz - ap.Ix) / ap.Iy;
    gamma6 = ap.Ixz / ap.Iy;
    gamma7 = (ap.Ix * (ap.Ix - ap.Iy) + ap.Ixz^2)/gamma;
    gamma8 = ap.Ix / gamma;
    
    mat1 = [gamma1*p*q - gamma2*q*r; %intermediate matrices
            gamma5*p*r - gamma6*(p^2 - r^2);
            gamma7*p*q - gamma1*q*r];
    mat2 = [gamma3*L + gamma4*N;
            M / ap.Iy;
            gamma4*L + gamma8*N];
    
    rates_dot = mat1 + mat2;
    
    % Return
    xdot = [pos_dot; orientation_dot; vel_dot; rates_dot];

end

