%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function E = TrueToEccentric(theta_star,e)
%TrueToEccentric takes in true anomaly and eccentricity and calculates
%eccentric anomaly
%INPUTS         theta_star      true anomaly of a satellite in radians
%               e               eccentricity of the orbit
%
%OUTPUTS        E              eccentric anomaly of the satellite in
%                              radians
coeff = sqrt((1-e)/(1+e));
right = coeff*tan(theta_star/2);
E = 2*atan(right);
end