%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function [a,e,i,RAAN,omega,theta_star] = CartesianToKepler(Rvec,Vvec,mu)
%CartesianToKepler takes in a Cartesian state vector and generates
% a set of Keplerian orbital elements, assuming an inertial frame centered on
%the primary body, using an equatorial plane of same
%OUTPUTS    R = three element XYZ Cartesian position vector     km
%           V = three element XYZ Cartesian velocity vector     km/s
%           mu = gravitational parameter of the central body    km^3/s^2
%
%INPUTS:    a = semimajor axis  (km)
%           e = eccentricity
%           i = inclination     (rad)
%           RAAN = right ascension of the ascending node (rad)
%           omega = argument of periapsis (rad)
%           theta_star = true anomaly   (rad)
%
%norm input vectors
R = norm(Rvec);     V = norm(Vvec);
%define h vector
Hvec = cross(Rvec,Vvec);  h = norm(Hvec);
%define orbital energy
energy = V^2/2 - mu/R;
%define semi-major axis
a = -mu/(2*energy);
%find eccentricity vector, unit vector, and scalar eccentricity
evec = cross(Vvec,Hvec)./mu - Rvec./R;
e = norm(evec);     e_hat = evec./e;
%find line of nodes vector and unit vector
nvec = cross([0;0;1], Hvec);  n_hat = nvec./norm(nvec);
%inclination
i = acos(Hvec(3)/h);
%RAAN
RAAN = acos(n_hat(1));
if n_hat(2) < 0
    RAAN = -RAAN;
end
%argument of periapsis
omega = acos(dot(n_hat,e_hat));
if e_hat(3) < 0
    omega = -omega;
end
%find true anomaly at this point
theta_star = acos(dot(Rvec./R, e_hat));
if dot(Rvec,Vvec) < 0
    theta_star = -theta_star;
end

end