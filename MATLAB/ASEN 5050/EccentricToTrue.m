%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function theta_star = EccentricToTrue(E,e)
%EccentricToTrue takes in eccentric anomaly and eccentricity and calculates
%true anomaly
%INPUTS         E       eccentric anomaly of a satellite in radians
%               e       eccentricity of the orbit
%
%OUTPUTS        theta_star     true anomaly of the satellite in
%                              radians
coeff = sqrt((1+e)/(1-e));
right = coeff*tan(E/2);
theta_star = 2*atan(right);
end