%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function [axis, angle] = AxisAngle(DCM)
% AxisAngle outputs an axis and an angle for a rotation defined by an input
% Direction Cosine Matrix.  Unity eigenvector is identified programatically
% using a tolerance of 1 +- 10^-3.  BEWARE of problems which might have
% another eigenvalue within this tolerance of 1.
%
%INPUTS     DCM         a 3x3 direction cosine matrix
%
%OUTPUTS    axis        a 3x1 unit vector representing the axis of rotation
%           angle       a scalar value of the angle of rotation in radians

[V,D]=eig(DCM);     %find eigenvalues/eigenvectors
[~,y] = find(abs(real(D) -1) < 10^-3);  %find the column # of unity eigen
axis = V(:,y);  %set axis to the corresponding column of eigenvector matrix
if y ==1        %find a column for non-unity eigenvalue
    angle_ind = 2;
elseif y == 2
    angle_ind = 3;
elseif y == 3
    angle_ind = 1;
end
%find angle using arccos on the real part of a complex eigenvalue
angle = acos(real(D(angle_ind,angle_ind)));
end