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

rho = stdatmo(-x(3));   %air density at altitude
ap = aircraft_parameters;   %for brevity

dt = u(4);          %control thrust, incorporated in Thrust below

V = norm(x(7:9));   %airspeed
Q = 0.5*rho*V^2;    %dynamic pressure
alpha = atan2(x(9),x(7));   %angle of attack
q_hat = (x(11)*ap.c)/(2*V);    %normalized q vector

Thrust = rho*ap.Sprop*ap.Cprop*(V + dt*(ap.kmotor - V))*dt*(ap.kmotor-V); %% Thrust model described in http://uavbook.byu.edu/lib/exe/fetch.php?media=shared:propeller_model.pdf

%find nondimensional coefficients for lift, drag and moment
CL = ap.CL0 + ap.CLalpha*alpha + ap.CLq*q_hat + ap.CLde*u(1);
CD = ap.CDmin + ap.K*(CL - ap.CLmin)^2;
Cm = ap.Cm0 + ap.Cmalpha*alpha + ap.Cmq*q_hat + ap.Cmde*u(1);

%find lift and drag forces and moment
lift = Q*ap.S*CL;    drag = Q*ap.S*CD;  moment = Q*ap.S*ap.c*Cm;

%convert lift and drag (in stability frame) to X and Z (in body frame)
aero_forces = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]*[-drag;-lift];

% Find final 
X = Thrust + aero_forces(1);
Z = aero_forces(2) ;
M = moment;

end
