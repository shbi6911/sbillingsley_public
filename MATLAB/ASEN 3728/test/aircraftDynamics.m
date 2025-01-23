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

ap = aircraft_parameters;

[X, Z, M] = lonAeroForcesAndMoments(x, u, ap);
[Y, L, N] = latAeroForcesAndMoments(x, u, ap);

% Replace with your code
o = x(4:6);
phi = o(1);
theta = o(2);
psi = o(3);

v_EB = x(7:9);
uE = v_EB(1);
vE = v_EB(2);
wE = v_EB(3);
% 
% p1 = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
%       cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
%       -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
% p_dot = p1*v_EB;

omega_B = x(10:12);
p = omega_B(1);
q = omega_B(2);
r = omega_B(3);

% o1 = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
%       0 cos(phi) -sin(phi);
%       0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
% o_dot = o1*omega_B;

R = rotation321(o);
T = attitudeInfluence321(o);

p_dot = transpose(R)*v_EB;
o_dot = T*omega_B;

% fg = R*[0; 0; ap.m*ap.g];
% 
% f = [X; Y; Z]+fg;

% V_dot = (f/ap.m)-cross(omega_B,v_EB);

% I_B = [ap.Ix 0 -ap.Ixz; 0 ap.Iy 0; -ap.Ixz 0 ap.Iz];
% 
% G_B = [L; M; N];

% omega_dot = I_B\(G_B-cross(omega_B,(I_B*omega_B)));

v1 = [r*vE-q*wE; p*wE-r*uE; q*uE-p*vE];
v2 = ap.g*[-sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi)];
v3 = [X; Y; Z]/ap.m;

V_dot = v1+v2+v3;
% 
Ix = ap.Ix; Iy = ap.Iy; Iz = ap.Iz; Ixz = ap.Ixz;

omega1 = [((Iy-Iz)/Ix)*q*r; ((Iz-Ix)/Iy)*p*r; ((Ix-Iy)/Iz)*p*q];
omega2 = [L/Ix; M/Iy; N/Iz];
omega_dot = omega1+omega2;

xdot = [p_dot; o_dot; V_dot; omega_dot];

end
