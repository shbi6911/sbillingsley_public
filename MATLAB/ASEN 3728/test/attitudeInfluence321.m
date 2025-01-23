function T = attitudeInfluence321(phiThetaPsi)
% ATTITUDEINFLUENCE321  Calculate the attitude influence matrix for a 3-2-1 rotation sequence.
% Inputs:
%   phiThetaPsi - 3x1 vector of Euler angles [phi; theta; psi] (rad)
% Outputs:
%   T - 3x3 attitude influence matrix

%break out input angles for convenience
phi = phiThetaPsi(1);   theta = phiThetaPsi(2); psi = phiThetaPsi(3);
%construct attitude influence matrix using given formula for 3-2-1
T = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
    0 cos(phi) -sin(phi);
    0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
end
