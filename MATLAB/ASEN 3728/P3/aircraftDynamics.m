function xdot = aircraftDynamics(x, u, aircraft_parameters)
%aircraftDynamics - computes the state derivative of a conventional aircraft
%
% Inputs:
%    x - state vector [x, y, z, phi, theta, psi, u, v, w, p, q, r] in SI units and radians
%    u - control vector [de, da, dr, dt] in radians
%    aircraft_parameters - struct containing nondimensional derivatives and other aircraft parameters
%
% Outputs:
%    xdot - state derivative vector [xdot, ydot, zdot, phidot, thetadot, psidot, udot, vdot, wdot, pdot, qdot, rdot] in SI units and radians

%x= [1;1;1;2;2;2;3;3;3;4;4;4];
% u = [0;0;0;0];
% aircraft_parameters = ttwistor;

ap = aircraft_parameters;
%break out state vector for convenient reference
p = x(1:3); % Position (m)
o = x(4:6); % Orientation (rad)
v_E_B = x(7:9); % Velocity (m/s)
omega_B = x(10:12); % Angular velocity (rad/s)

%break out input angles for convenience
phi = o(1);   theta = o(2); psi = o(3);
%get rotation matrix from body to inertial
R = rotation321(o);
%construct attitude influence matrix using given formula for 3-2-1
T = attitudeInfluence321(o);
%construct inertia matrix
I_B = [ap.Ix 0 ap.Ixz;0 ap.Iy 0; ap.Ixz 0 ap.Iz];

%call functions to get aero forces and moments
[X, Z, M] = lonAeroForcesAndMoments(x, u, ap);
[Y, L, N] = latAeroForcesAndMoments(x, u, ap);

% Construct derivative state vector
xdot = zeros(12, 1);    %preallocate
xdot(1:3) = R'*v_E_B;   %derivative of position is velocity, rotated
xdot(4:6) = T*omega_B;  %attitude influence times angular rate
%acceleration is omega cross v + gravity + aero forces over mass
xdot(7:9) = +...-cross(omega_B,v_E_B)
    (ap.g*[(-sin(theta));(cos(theta)*sin(phi));(cos(theta)*cos(phi))])+...
    (1/ap.m)*[X;Y;Z];
%angular acceleration is inverse inertia matrix times cross product and
%moment
xdot(10:12) = I_B\([L;M;N]-cross(omega_B,(I_B*omega_B)));

end
