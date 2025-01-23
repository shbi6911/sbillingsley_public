%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%% Question 1
clear; close all; const = getConstOrbitz;
V = [-4;3;-5];
uhat_r = [0.26726;0.53452;0.80178];
vhat = V./norm(V);
gamma = 90-acosd(dot(vhat,uhat_r));

hhat = cross(uhat_r,vhat);
theta_hat = cross(uhat_r,hhat);

v_r = dot(V,uhat_r);
v_theta = norm(cross(uhat_r,V));
gamma_check = rad2deg(atan2(v_r,v_theta));

%% Question 2
const = getConstOrbitz;
h = 60000;
r = h^2/const.Earth.mu;
T_2 = ((2*pi)/sqrt(const.Earth.mu))*r^(3/2);
T_hrs_2 = T_2/3600;

%% Question 3
const = getConstOrbitz;
apogee = 2000;
perigee = 500;
r_a = apogee + const.Earth.radius;
r_p = perigee + const.Earth.radius;
e = (r_a - r_p)/(r_p + r_a);
h = sqrt(r_p*const.Earth.mu*(1+e));
v_p = h/r_p;
v_a = h/r_a;
a = (r_a+r_p)/2;
T_3 = 2*pi*sqrt((a^3)/const.Earth.mu);
T_hrs_3 = T_3/3600;
T_min_3 = T_3/60;

%% Question 4
const = getConstOrbitz;
v_p_4 = 10;
theta = deg2rad(120);
h_4 = r_p*v_p_4;
e_4 = (h_4^2/(r_p*const.Earth.mu))-1;
a_4 = r_p/(1-e_4);
r_a_4 = (2*a_4) - r_p;
r_4 = (a_4*(1-e_4^2)/(1+(e_4*cos(theta))));
alt_4 = r_4 - const.Earth.radius;
gamma_4 = rad2deg(atan2((e_4*sin(theta)),(1+(e_4*cos(theta)))));

%% Question 5
const = getConstOrbitz;
R_5_vec = [-5613.648;142.5547;-4995.530];
V_5_vec = [-0.9121618;-0.7208981;1.596800];
h_5_vec = cross(R_5_vec,V_5_vec);
e_5_vec = (cross(V_5_vec,h_5_vec)/const.Mars.mu) - (R_5_vec./norm(R_5_vec));
energy_5 = (norm(V_5_vec)^2/2) - (const.Mars.mu/norm(R_5_vec));

%% Question 6
const = getConstOrbitz;
r_6 = 4953.786;
v_6 = 3.130917;
gamma_6 = -20.230104;
v_6_r = v_6*sind(gamma_6);
v_6_perp = v_6*cosd(gamma_6);
v_6_vec = [v_6_r;v_6_perp;0];
h_6 = r_6*v_6_perp;
p_6 = h_6^2/const.Mars.mu;
energy_6 = (v_6^2/2) - (const.Mars.mu/r_6);
e_6 = sqrt(1+((2*h_6^2*energy_6)/const.Mars.mu^2));
r_p_6 = p_6/(1+e_6);
r_a_6 = p_6/(1-e_6);
a_6 = (r_a_6+r_p_6)/2;
theta_6 = -rad2deg(acos((p_6-r_6)/(r_6*e_6)));
T_6 = 2*pi*sqrt(a_6^3/const.Mars.mu);
T_6_hrs = T_6/3600;
b_6 = a_6*sqrt(1 - e_6^2);
