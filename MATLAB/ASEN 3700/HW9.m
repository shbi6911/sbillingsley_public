%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%% Q1
clear;
const = getConstOrbitz();
R = [6500;-7500;-2500];
V = [4;3;-3];
r = norm(R);
v = norm(V);
%find v_r to check quadrant
v_r = dot(R,V)/r;
%positive result means going from periapsis to apoapsis
H = cross(R,V);
h = norm(H);
i = rad2deg(acos(H(3)/h));
%line of nodes
N = cross([0;0;1],H);
n = norm(N);
omega = rad2deg(acos(N(1)/n));
E1 = cross(V,H)./const.Earth.mu;
E2 = R./r;
E = E1 - E2;
e = norm(E);
a = (h^2/const.Earth.mu)*(1/(1-e^2));
omega = rad2deg(acos(dot(N,E)/(n*e)));
theta = rad2deg(acos(dot(E,R)/(r*e)));

%% Q2
clear;
const = getConstOrbitz;
R = [-6600;-1300;-5200];
E = [-0.4;-0.5;-0.6];
W = cross(R,E)./norm(cross(R,E));
i = acosd(W(3));

%% Q3
clear;  const = getConstOrbitz;
a = 7016;   e = 0.05;   i = deg2rad(45);    Omega = 0;  omega=deg2rad(20);
theta = deg2rad(10);    phi = omega + theta;
r = (a*(1-e^2))/(1+(e*cos(theta)));
R3_omega = [cos(-Omega) sin(-Omega) 0;
            -sin(-Omega) cos(-Omega) 0;
            0           0           1];
R3_phi = [cos(-phi) sin(-phi) 0;
            -sin(-phi) cos(-phi) 0;
            0           0           1];
R1_i = [1       0       0; 
        0 cos(-i) sin(-i);
        0 -sin(-i) cos(-i)];
Q = R3_omega*R1_i*R3_phi;
disp(Q);
R_inertial = Q*[r;0;0];

%% Q4
clear;  const = getConstOrbitz;
R = [6.342550e6;-6.544023e5;3.284575e4];    %km
V = [-3.828922;1.940090;-2.156209e-2];      %km/s
r = norm(R);
v = norm(V);
%find v_r to check quadrant
v_r = dot(R,V)/r;
%negative result means going from apoapsis to periapsis
H = cross(R,V);
h = norm(H);
i = rad2deg(acos(H(3)/h));
%line of nodes
N = cross([0;0;1],H);
n = norm(N);
Omega = -rad2deg(acos(N(1)/n));
E1 = cross(V,H)./const.Jupiter.mu;
E2 = R./r;
E = E1 - E2;
e = norm(E);
a = (h^2/const.Jupiter.mu)*(1/(1-e^2));
omega = -rad2deg(acos(dot(N,E)/(n*e)));
theta = -rad2deg(acos(dot(E,R)/(r*e)));
v_theta = sqrt(v^2 - v_r^2);

R3_theta = [cosd(theta) sind(theta) 0;
            -sind(theta) cosd(theta) 0;
            0           0           1];
V_per = R3_theta'*[v_r;v_theta;0];
R_per = R3_theta*[r;0;0];

thingy = sqrt((1-e)/(1+e));
T = ((2*pi)/sqrt(const.Jupiter.mu))*a^(3/2);
E_0 = 2*atan(thingy*tand(theta/2));
M_0 = E_0 - e*sin(E_0);
t_0 = (M_0*T)/(2*pi);
delta_t = 14*24*3600;
t_1 = t_0 + delta_t;
M_1 = (2*pi*t_1)/T;
E_1 = kepler_E(e,M_1);
theta_1 = 2*atan((1/thingy)*tan(E_1/2));
r_1 = (a*(1-e^2))/(1+(e*cos(theta_1)));
v_1_r = (const.Jupiter.mu/h)*e*sin(theta_1);
v_1_theta = (const.Jupiter.mu/h)*(1+e*cos(theta_1));

phi = omega + rad2deg(theta_1);
R3_omega = [cosd(-Omega) sind(-Omega) 0;
            -sind(-Omega) cosd(-Omega) 0;
            0           0           1];
R3_phi = [cosd(-phi) sind(-phi) 0;
            -sind(-phi) cosd(-phi) 0;
            0           0           1];
R1_i = [1       0       0; 
        0 cosd(-i) sind(-i);
        0 -sind(-i) cosd(-i)];
Q = R3_omega*R1_i*R3_phi;
R_inertial = Q*[r_1;0;0];
V_inertial = Q*[v_1_r;v_1_theta;0];