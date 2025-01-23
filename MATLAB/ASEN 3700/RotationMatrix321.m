%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function DCM = RotationMatrix321(angles)
%RotationMatrix321 generates a rotation matrix for a 3-2-1 rotation in a
%Euclidean x-y-z space
%
% INPUTS:   alpha:  Rotation about the x-axis [degrees]
%           beta:   Rotation about the y-axis [degrees]
%           gamma:  Rotation about the z-axis [degrees]
% 
% OUTPUTS:  DCM:    3x3 3-2-1 Rotation Matrix
alpha = angles(1);  beta = angles(2);   gamma = angles(3);
R1 = [1 0 0; 0 cosd(gamma) sind(gamma); 0 -sind(gamma) cosd(gamma)];
R2 = [cosd(beta) 0 -sind(beta); 0 1 0; sind(beta) 0 cosd(beta)];
R3 = [cosd(alpha) sind(alpha) 0; -sind(alpha) cosd(alpha) 0; 0 0 1];
DCM = R1 * R2 * R3;
end