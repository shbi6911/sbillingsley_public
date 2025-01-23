%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function E = kepler_E(e, M)
% kepler_E iteratively solves Kepler's equation for a given eccentricity and
% mean motion, using an intial guess defined by mean motion
%
%INPUTS:    e       orbital eccentricity
%           M       orbital mean motion in seconds
%
%OUTPUTS:   E       eccentric anomaly    
error = 1e-8;
if M < pi
    E = M+ e/2;
else
    E = M-e/2;
end
ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1-e*cos(E));
    E = E-ratio;
end
end