function [X, Z, M] = lonAeroForcesAndMoments(x, u, aircraft_parameters)
%lonAeroForcesAndMoments - Calculate the longitudinal aerodynamic forces and moments
%
% Inputs:
%    x - state vector [x, y, z, phi, theta, psi, u, v, w, p, q, r] in SI units and radians
%    u - control vector [de, da, dr, dt] in radians
%    aircraft_parameters - struct containing nondimensional derivatives and other aircraft parameters
%
% Outputs:
%    X - force in x direction in N
%    Z - force in z direction in N
%    M - moment in y direction in Nm


ap = aircraft_parameters;

rho = stdatmo(-x(3));

de = u(1);
dt = u(4);

V = norm(x(7:9));
alpha = atan2(x(9), x(7));
Q = 0.5*rho*V^2;

q = x(11);

sa = sin(alpha);
ca = cos(alpha);

CL = ap.CL0 + ap.CLalpha*alpha + ap.CLq*q*ap.c/(2*V) + ap.CLde*de;
CD = ap.CDmin + ap.K*(CL - ap.CLmin)^2;

Thrust = rho*ap.Sprop*ap.Cprop*(V + dt*(ap.kmotor - V))*dt*(ap.kmotor-V); %% Thrust model described in http://uavbook.byu.edu/lib/exe/fetch.php?media=shared:propeller_model.pdf

CX = -CD*ca + CL*sa;
CZ = -CD*sa - CL*ca;

X = Q*ap.S*CX + Thrust;
Z = Q*ap.S*CZ;

Cm = ap.Cm0 + ap.Cmalpha*alpha + ap.Cmq*q*ap.c/(2*V) + ap.Cmde*de;
M = Q*ap.S*ap.c*Cm;

end
