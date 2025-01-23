%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%Homework 3

%Problem 11.21
%define givens
I_G = [20 -10 0;-10 30 0;0 0 40];
omega = [10;20;30];
H_G = I_G*omega;
%disp(H_G);
KE = 0.5*dot(omega,H_G);
%disp(KE);

%Problem 4
%define givens
angles_deg = [60;-30;120];  angles_rad = angles_deg'*(pi/180);
DCM_321 = RotationMatrix321(angles_deg);
%DCM_321_check = eul2rotm(angles_rad);
%disp(DCM_321);
%disp(DCM_321_check);

DCM_313 = RotationMatrix313(angles_deg);
%DCM_313_check = eul2rotm(angles_rad,"ZXZ");
%disp(DCM_313);
%disp(DCM_313_check);

[axis_321,angle_321] = AxisAngle(DCM_321); %find axis and angle
angle_321_deg = angle_321*(180/pi); %convert to degrees
%disp(axis_321);
%disp(angle_321_deg);

[axis_313,angle_313] = AxisAngle(DCM_313); %find axis and angle
angle_313_deg = angle_313*(180/pi); %convert to degrees
%disp(axis_313);
%disp(angle_313_deg);

quat_321 = Quat(axis_321,angle_321);
%disp(quat_321);
quat_313 = Quat(axis_313,angle_313);
%disp(quat_313);

Z = [0;0;1]; %set originial z-axis
Z_321 = DCM_321*Z;  %new Z-axis
Z_angle_321 = acosd(dot(Z_321,Z));  %find angle between old and new
%disp(Z_angle_321);

Z_313 = DCM_313*Z;  %new Z-axis
Z_angle_313 = acosd(dot(Z_313,Z));  %find angle between old and new
%disp(Z_angle_313);

%Problem 5
%define givens
omega = [1.25;0.5;0.25];
omega_321 = DCM_321'*omega;
%disp(omega_321);
omega_313 = DCM_313'*omega;
%disp(omega_313);

%Problem 6
%define givens
angles_P6_deg = [-35;45;135];
angles_P6_rad = angles_P6_deg.*(180/pi);
omega_P6_deg = [20;30;10];
omega_P6_rad = angles_P6_rad.*(180/pi);

omega_P6_body = RotationMatrix313(angles_P6_deg)*omega_P6_deg;
disp(omega_P6_body);

%% functions
function quaternion = Quat(axis,angle)
%INPUTS     axis    vector representing axis of rotation
%           angle   scalar value of rotation angle in radians
%OUTPUTS    quaternion      associated quaternion for the rotation,
%                           represented as a 4-element column vector
if isrow(axis)
    axis = axis';   %check and convert to column vector
end
quaternion = [sin(angle/2).*axis;cos(angle/2)];
end

function [a_tilde] = tilde(a)
%INPUTS     a           a 3-element vector in Euclidean space
%
%OUTPUTS    a_tilde     the corresponding cross product matrix

a_tilde = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end

function [axis, angle] = AxisAngle(DCM)
%INPUTS     DCM         a 3x3 direction cosine matrix
%OUTPUTS    axis        a 3x1 unit vector representing the axis of rotation
%           angle       a scalar value of the angle of rotation in radians

[V,D]=eig(DCM);     %find eigenvalues/eigenvectors
[~,y] = find(abs(real(D) -1) < 10^-3);  %find the column # of unity eigen
axis = V(:,y);  %set axis to the corresponding column of eignevector matrix
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

