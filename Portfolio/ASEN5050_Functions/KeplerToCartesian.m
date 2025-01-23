%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function [R,V] = KeplerToCartesian(a,e,i,RAAN,omega,theta_star, mu)
%KeplertoCartesian takes in a set of Keplerian orbital elements and
%generates a Cartesian state vector, assuming an inertial frame centered on
%the primary body of the Keplerian orbit, using an equatorial plane of same
%
%INPUTS:    a = semimajor axis  (km)
%           e = eccentricity
%           i = inclination     (rad)
%           RAAN = right ascension of the ascending node (rad)
%           omega = argument of periapsis (rad)
%           theta_star = true anomaly   (rad)
%           mu = gravitational parameter of the central body    km^3/s^2
%
%OUTPUTS    R = three element XYZ Cartesian position vector     km
%           V = three element XYZ Cartesian velocity vector     km/s
r = (a*(1-e^2))/(1+(e*cos(theta_star)));    %find radius from true anomaly
R_rot = [r;0;0];    %radius vector in rotating frame
h = sqrt(mu*a*(1-e^2)); %define specific angular momentum
v_r_2 = (mu/h)*e*sin(theta_star);   %radial velocity
v_theta_2 = (mu/h)*(1+(e*cos(theta_star))); %transverse velocity
V_rot = [v_r_2;v_theta_2;0];   %velocity vector in rotating frame

%define quantities for rotation matrix
co = cos(RAAN);    ci = cos(i);    ct = cos(theta_star + omega);
so = sin(RAAN);    si = sin(i);    st = sin(theta_star + omega);
%define rotation matrix
DCM = [((co*ct)-(so*ci*st)),((-co*st)-(so*ci*ct)),(so*si);
       ((so*ct)+(co*ci*st)),((-so*st)+(co*ci*ct)),(-co*si);
       (si*st), (si*ct), (ci)];
R = DCM*R_rot;  %rotate into Cartesian
V = DCM*V_rot;
end