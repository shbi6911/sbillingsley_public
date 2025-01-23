% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% RotationMatrix321
% Created: 1/31/24

function DCM = RotationMatrix321(attitude321)
% INPUTS:   attitude312: 3x1 array of angles [deg]
% 
% OUTPUTS:  DCM:    3x3 3-2-1 Rotation Matrix

alpha = attitude321(1);
beta = attitude321(2);
gamma = attitude321(3);
R1 = [1 0 0; 0 cosd(alpha) sind(alpha); 0 -sind(alpha) cosd(alpha)];
R2 = [cosd(beta) 0 -sind(beta); 0 1 0; sind(beta) 0 cosd(beta)];
R3 = [cosd(gamma) sind(gamma) 0; -sind(gamma) cosd(gamma) 0; 0 0 1];
DCM = R1 * R2 * R3;

end