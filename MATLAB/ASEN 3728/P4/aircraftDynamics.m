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

p_E = x(1:3);
o = x(4:6);
v_E_B = x(7:9);
omega = x(10:12);

ap = aircraft_parameters;

I_B = [ap.Ix 0 -ap.Ixz;
       0 ap.Iy 0;
       -ap.Ixz 0 ap.Iz];

[X, Z, M] = lonAeroForcesAndMoments(x, u, ap);
[Y, L, N] = latAeroForcesAndMoments(x, u, ap);

f_g_E = [0; 0; 9.81]*ap.m;
f_aero_B = [X; Y; Z];
f_B = f_aero_B + rotation321(o)*f_g_E;
G_B = [L; M; N];

pdot_E = rotation321(o)'*v_E_B;
odot = attitudeInfluence321(o)*omega;
vdot_E_B = -cross(omega, v_E_B) + f_B/ap.m;
omegadot_B = I_B\(-cross(omega, I_B*omega) + G_B);

xdot = [pdot_E; odot; vdot_E_B; omegadot_B];
end
