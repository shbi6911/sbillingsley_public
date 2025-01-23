%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%Question 11.7

    %define given matrix
sys = [10 1 1 1; 10 -1 -1 -1; 8 4 -4 4; 8 -2 2 -2; 12 3 -3 -3; 12 -3 3 3];
m_sys = sum(sys(:,1));  %total mass
cm_vec = [(sum(sys(:,1).*(sys(:,2))))/m_sys;(sum(sys(:,1).*(sys(:,3))))/m_sys;...
    (sum(sys(:,1).*(sys(:,4))))/m_sys]; %position of center mass
    %move all position vectors to center mass
sys_adjust = sys - (ones(size(sys)).*[0,cm_vec']);

I_matrix = zeros(3,3);  %preallocate
for i = 1:length(sys_adjust)    %sum up point moments of inertia
    I_matrix = I_matrix + inertia_matrix(sys_adjust(i,1),sys_adjust(i,2:4));
end
%disp(I_matrix);

    %Question 11.8
u = [1;2;2];    u_hat = u/norm(u);  %unit axis vector
I_matrix_origin = zeros(3,3);  %preallocate
for i = 1:length(sys)    %sum up point moments of inertia about the origin
    I_matrix_origin = I_matrix_origin + inertia_matrix(sys(i,1),sys(i,2:4));
end
I_u = u_hat'*I_matrix_origin*u_hat; %calculate moment of inertia about u-hat
disp(I_u);

    %Question 3
[V,D] = eig(I_matrix);
%disp(diag(D));
I_PA = V'*I_matrix*V;
%disp(I_PA);

    %Question 4
    %define constants
m_cube = 4; m_sphere = 2;   d = 2;  s_cube = 0.5;   r_sphere = 0.4;
    %Part A
    %sphere inertia tensor about sphere CM
I_sphere = (2/5)*m_sphere*r_sphere^2;
I_G_sphere = [I_sphere 0 0;0 I_sphere 0; 0 0 I_sphere];
    %cube inertia tensor about cube CM
I_cube = (1/6)*m_cube*s_cube^2;
I_G_cube = [I_cube 0 0;0 I_cube 0; 0 0 I_cube];
    %vector from cube to sphere
R_21 = [0;(d+(s_cube/2)+r_sphere);0];
    %composite inertia tensor about rigid body CM
I_G_total = I_G_cube + I_G_sphere -...
    (((m_cube*m_sphere)/(m_cube+m_sphere))*tilde(R_21)*tilde(R_21));
%disp(I_G_total);
    %Part B
alpha = 10/3600;    %desired angular acceleration in rad/s^2
moment = -alpha*I_G_total(3,3);
%disp(moment);
    %Part C
vec = [0;sind(45);cosd(45)];
I_vec = vec'*I_G_total*vec;
%disp(I_vec);

    %Question 5
    %Part A
    %define constants
m_ell = 500;    a = 3; b = 1.25;  c = 1;
I_1 = (m_ell/5)*(b^2+c^2);
I_2 = (m_ell/5)*(c^2+a^2);
I_3 = (m_ell/5)*(a^2+b^2);
I_PA_ell = diag([I_1;I_2;I_3]);
%disp(I_PA_ell);
    %Part B
    %define constants
w_1 = 5;    w_2 = -3;   w_3 = 1;    omega = [w_1;w_2;w_3];
H_G = I_PA_ell*omega;
KE = 0.5*dot(omega,H_G);
H_G_mag = norm(H_G);
%disp(H_G_mag);
%disp(KE);
    %Part C
moment_ell = tilde(omega)*I_PA_ell*omega;
%disp(moment_ell);

function [a_tilde] = tilde(a)
%INPUTS     a           a 3-element vector in Euclidean space
%
%OUTPUTS    a_tilde     the corresponding cross product matrix

a_tilde = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end

function I_matrix = inertia_matrix(m,pos)
x= pos(1); y = pos(2); z = pos(3);
I_matrix = [m*(y^2 + z^2),-m*x*y,-m*x*z;
            -m*x*y, m*(x^2+z^2), -m*y*z;
            -m*x*z, -m*y*z, m*(x^2+y^2)];
end