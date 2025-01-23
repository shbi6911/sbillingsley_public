%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function DCM = RotationMatrix313(angles)
%RotationMatrix313 generates a rotation matrix for a 3-1-3 rotation in a
%Euclidean x-y-z space
%
%INPUTS     attitude313 a column vector of Euler angles in order 3-1-3
%
%OUTPUTS    DCM         a 3x3 Direction Cosine Matrix for the specified
%                       rotation
%
%METHODOLOGY    RotationMatrix313 takes in a rotation specified by a 3-1-3
%sequence of Euler angles and outputs the associated DCM.
alpha = angles(1);  beta = angles(2);   gamma = angles(3);
R3_1 = [cosd(gamma) sind(gamma) 0; -sind(gamma) cosd(gamma) 0; 0 0 1];
R2 = [1 0 0; 0 cosd(beta) sind(beta); 0 -sind(beta) cosd(beta)];
R3_2 = [cosd(alpha) sind(alpha) 0; -sind(alpha) cosd(alpha) 0; 0 0 1];
DCM = R3_1*R2*R3_2;
end