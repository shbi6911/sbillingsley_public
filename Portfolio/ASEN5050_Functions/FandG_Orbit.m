%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function [orbit_state] = FandG_Orbit(Rvec,Vvec,mu,n)
%FandG_Orbit generates a 6xn matrix of state vectors encompassing a single
%2BP orbit, using f and g functions
%INPUTS     Rvec        3x1 radius vector of Cartesian elements, assumed
%                       to be in and inertial equatorial frame
%           Vvec        3x1 velocity vector in same frame as Rvec
%           mu          gravitational parameter of central body in 2BP
%                       orbit
%           n           number of points along the orbit to be calculated
%                       and output in the orbit_state matrix
%
%OUTPUTS    orbit_state 6xn matrix of state values where each 6x1 column is
%                       a state vector of the form [x;y;z;vx;vy;vz]

r_0 = norm(Rvec);     v_0 = norm(Vvec);     %magnitudes
Hvec = cross(Rvec,Vvec);    h = norm(Hvec); %specific angular momentum
p = h^2/mu;                                 %semi-latus rectum
evec = cross(Vvec,Hvec)./mu - Rvec./r_0;      %eccentricity
e = norm(evec);     
if e == 0
    e_hat = [1;0;0];
else
    e_hat = evec./e;
end
a = p/(1-e^2);                              %semimajor axis
theta_star = acos(dot(Rvec./r_0, e_hat));     %true anomaly
if dot(Rvec,Vvec) < 0                       %true anomaly quadrant check
    theta_star = -theta_star;
end % if quad check
theta = [linspace(0,pi,n/2),linspace(-pi,0,n/2)];
delta_theta = theta - theta_star;

r1 = p./(1+e.*cos(theta));

%find radius and velocity vectors at impact using f and g functions
f = 1-(r1./p).*(1-cos(delta_theta));
g = (r1.*r_0.*sin(delta_theta))./sqrt(mu*p);
coeff = ((1-cos(delta_theta))./p) - (1./r1) - (1./r_0); 
f_dot = sqrt(mu/p).*tan(delta_theta./2).*coeff;
g_dot = 1 - (r_0./p).*(1-cos(delta_theta));

Rvec1 = Rvec*f + Vvec*g;                    %calculate arrays of vectors
Vvec1 = Rvec*f_dot + Vvec*g_dot;
orbit_state = [Rvec1;Vvec1];
end %function