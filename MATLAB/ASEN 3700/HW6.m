%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

a = 1/sqrt(2);

Q_BI_1 = SunAndMoon([0;-1;0],[0;a;a],[1;0;0],[-a;-a;0]);

[axis1, angle1] = AxisAngle(Q_BI_1);
angle1deg = rad2deg(angle1);

Q_BI_2 = SunAndMoon([0;-1;0],[0;a;a],[a;0;-a],[0;0;1]);

[axis2, angle2] = AxisAngle(Q_BI_2);
angle2deg = rad2deg(angle2);

Q_BI_3 = SunAndMoon([0;-1;0],[0;a;a],[0;a;-a],[0;0;1]);

[axis3, angle3] = AxisAngle(Q_BI_3);
angle3deg = rad2deg(angle3);

function [Q] = SunAndMoon(s_I,m_I,s_b,m_b)

t_I1 = s_I;
t_I2 = cross(s_I,m_I)/norm(cross(s_I,m_I));
t_I3 = cross(t_I1,t_I2)/norm(cross(t_I1,t_I2));
Q_TI = [t_I1 t_I2 t_I3];

t_B1 = s_b;
t_B2 = cross(s_b,m_b)/norm(cross(s_b,m_b));
t_B3 = cross(t_B1,t_B2)/norm(cross(t_B1,t_B2));
Q_TB = [t_B1,t_B2,t_B3];

Q = Q_TI*Q_TB';

end

function [axis, angle] = AxisAngle(DCM)
%INPUTS     DCM         a 3x3 direction cosine matrix
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
