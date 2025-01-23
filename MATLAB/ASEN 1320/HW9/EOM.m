%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function [statechange] = EOM(T,V,P)
%EOM This function takes in a time step and projectile state, and computes
%the rate of change of the projectile's state over the specified time step
% Input: Time step T, initial projectile state V = [vx, vy, x, y]
%        parameters P = [Cd, radius, mass, rho, g]
%Output: rates of change at current time step [dvx/dt, dvy/dt, vx, vy]
D = 0.5*P(4)*P(1)*(pi*P(2)*P(2))*(V(1)^2+V(2)^2);
theta = atan2(V(2),V(1));
statechange = zeros(4,1);
statechange(1,1) = -((D*cos(theta))/P(3));
statechange(2,1) = (-((D*sin(theta))/P(3))-P(5));
statechange(3,1) = V(1);
statechange(4,1) = V(2);
end