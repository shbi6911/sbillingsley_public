%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc; const = getConst();
%givens
r_p = 600000;           %km
r_a = 1800000;
a_titan = 1221830;
r_titan = 2575;

m_titan = 1.3455e23;    %kg

%characterize spacecraft orbit around Saturn
a_1 = (r_p+r_a)/2;
e_1 = 1 - (r_p/a_1);
p_1 = a_1*(1 - e_1^2);
h_1 = sqrt(p_1*const.mu.Saturn);
%find true anomaly at Titan orbit intersection
theta_star_1 = -acos((p_1 - a_titan)/(a_titan*e_1));

%find velocity and flight path angle
v_1 = sqrt(const.mu.Saturn*((2/a_titan)-(1/a_1)));
gamma = -acos(h_1/(a_titan*v_1));
%find components of velocity and express as a vector
v_r1 = v_1*sin(gamma);   v_theta1 = v_1*cos(gamma);
V1_roh = [v_r1;v_theta1;0];

%velocity of Titan relative to Saturn
v_titan = sqrt(const.mu.Saturn/a_titan);
V_titan = [0;v_titan;0];

V_infty_in = V1_roh - V_titan;  %velocity of s/c relative to Titan

%characterize hyperbolic trajectory relative to Titan
r_p_h = 3000;       %km
mu_titan =  m_titan*const.G;
a_h = -mu_titan/norm(V_infty_in)^2;
e_h = 1 - (r_p_h/a_h);
delta = 2*asin(1/e_h);

%find outgoing velocity relative to Saturn
beta_1 = pi - acos(dot(V_infty_in,V_titan)/(norm(V_infty_in)*norm(v_titan)));
beta_2 = beta_1 - delta;
V_infty_out = [-norm(V_infty_in)*sin(beta_2);-norm(V_infty_in)*cos(beta_2);0];
V_out = V_titan + V_infty_out;

%characterize new orbit
h_2 = norm(cross([a_titan;0;0],V_out));
a_2 = -const.mu.Saturn/(dot(V_out,V_out) - (2*const.mu.Saturn)/a_titan);
e_2 = sqrt(1 - (h_2^2/(const.mu.Saturn*a_2)));
p_2 = a_2*(1-e_2^2);
theta_star_2 = -acos((p_2 - a_titan)/(a_titan*e_2));
%find equivalent impulsive maneuver
delta_V = V_out - V1_roh;