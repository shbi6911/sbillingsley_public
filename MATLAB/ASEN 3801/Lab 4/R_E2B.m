%% ASEN 3801 Lab 2, Function RotationMatrix321
% Notation, Orientation, Rotation Matrices, & Direction Cosine Matrices
% Group Members: Milo Casey, Kyle Goodall, Jacob
% Lab Section-Group Number: 001-12
% Date Created: 2/2/2024
% Last Modified: 2/10/2024

function DCM = R_E2B(euler)
% Inputs:       3 x 1 matrix of Euler angles in radians
%
% Outputs:      Rotation DCM
%
% Methodology:  Use individual DCM's to create one 321 DCM

% Define Variables from attitude vector (in radians)
    phi = euler(1);
    theta = euler(2);
    psi = euler(3);

% Principal Axis Rotation Matrices
    R1phi = [1     0           0;
              0     cos(phi)   sin(phi);
              0     -sin(phi)  cos(phi)];

    R2theta = [ cos(theta)    0     -sin(theta);
               0            1     0;
               sin(theta)    0     cos(theta)];

    R3psi = [ cos(psi)      sin(psi)     0;
                -sin(psi)     cos(psi)      0;
                0               0               1];

% DCM for inertial to body rotation (when 3-2-1 Euler angles are used)
    DCM = R1phi*R2theta*R3psi;
end

