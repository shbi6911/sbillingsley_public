function R = rotation321(phiThetaPsi)
% ROTATION321 Creates a rotation matrix from Euler angles for a 3-2-1 rotation sequence.
% Inputs:
%   phiThetaPsi = 3x1 vector of Euler angles [phi; theta; psi] (rad)
% Outputs:
%   R = 3x3 direction cosine matrix

%break out input angles for convenience
phi = phiThetaPsi(1);   theta = phiThetaPsi(2); psi = phiThetaPsi(3);
%calculate individual rotation matrices
R_phi = [1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
R_theta = [cos(theta) 0 -sin(theta);0 1 0; sin(theta) 0 cos(theta)];
R_psi = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0;0 0 1];
%multiply rotation matrices to arrive at final DCM
R = R_phi*R_theta*R_psi;
end
