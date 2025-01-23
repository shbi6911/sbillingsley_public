%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% 9-10-24
clear; clc;
X1 = [4.7904e7; -7.6993e7; -3.8208e6];
V1 = [3.2840e1; -8.6630; -2.0141];

X2 = [1.0064e8; -6.8707e7; -6.7414e6];
V2 = [9.5730; 1.1668e1; -3.8757e-1];

h_1 = cross(X1, V1);
disp(norm(h_1));
h_2 = cross(X2, V2);
disp(norm(h_2));

const.mu.Sun = 1.32712428e11;

energy1 = dot(V1,V1)/2 - const.mu.Sun/(norm(X1));
energy2 = dot(V2,V2)/2 - const.mu.Sun/(norm(X2));
disp(energy1);
disp(energy2);

%% 9-17-24
clear; clc;
R_1 = [9420.5; 9420.5; 0];
n_hat_1 = -R_1./norm(R_1);
e_hat_1 = [-0.8624;-0.3624;-0.3535];
omega = acos(dot(n_hat_1,e_hat_1));
omega = -omega;
omega_deg = rad2deg(omega);
disp(omega_deg);

%% 10-8-24
clear; clc;
%givens
m = 1224;       %mass in kg
m_p = 157;      %propellant mass in kg
I_sp = 212;     %specific impulse of hydrazine propellant system
g0 = 9.81;
const.mu.Moon = 4902.799;
%initial orbit
a1 = 1938;      e1 = 0;
%final orbit
v2 = 1.5905;    gamma_2 = deg2rad(-1.4784);

%velocity magnitude of initial circular orbit
v1 = sqrt(const.mu.Moon/a1);
%find delta-v
delta_v = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(gamma_2))*1000;
%find propellant used
prop = m*(1-exp(-delta_v/(I_sp*g0)));

%% 10-24-24
clear; clc; const = getConst();
%givens
a_t = 5.75*const.AU;    e_t = 0.8104;   r1 = const.a.Mars*const.AU;
%velocity of Mars in the heliocentric frame
v_Mars = sqrt(const.mu.Sun/r1);

theta_star_t = acos((a_t*(1-e_t^2) - r1)/(r1*e_t));
v1 = sqrt(((2*const.mu.Sun)/r1) - (const.mu.Sun/(a_t)));
h_t = sqrt((a_t*(1-e_t^2))*const.mu.Sun);
gamma_1 = acos(h_t/(r1*v1));
vr_1 = v1*sin(gamma_1);     vt_1 = v1*cos(gamma_1);
V_in = [vr_1;vt_1;0];
V_Mars = [0;v_Mars;0];

V_inf_in = V_in - V_Mars;

r_p = 200 + const.radius.Mars;

a_hyper = -const.mu.Mars/dot(V_inf_in,V_inf_in);

e_hyper = 1-(r_p/a_hyper);
delta = 2*asin(1/e_hyper);

V_inf_theta = -norm(V_inf_in)*sin(delta/2);
V_inf_r = norm(V_inf_in)*cos(delta/2);
V_out = V_Mars + [V_inf_r;V_inf_theta;0];

%% 11-12-24
clear; clc; const = getConst;
%givens
J2_mars = 0.001964;

P_mars = 2*pi*sqrt((const.a.Mars*const.AU)^3/const.mu.Sun);
P_mars_days = P_mars/(3600*24);
theta_dot = (2*pi)/P_mars;

%find semimajor axis of satellite orbit
P_sat = 118*60;
a_sat = (((P_sat/(2*pi))^2)*const.mu.Mars)^(1/3);
%find inclination using solved J2 equation
i_sat = acos((theta_dot*a_sat^(7/2))/(sqrt(const.mu.Mars)*J2_mars*const.radius.Mars^2)*(-2/3));
i_sat_deg = rad2deg(i_sat);

J2_const_mars = (sqrt(const.mu.Mars)*J2_mars*const.radius.Mars^2);
i2_sat = asin(16/25);
thing = (3/2)*(J2_const_mars/(theta_dot*a_sat^(7/2)))*cos(i2_sat)+1;
e2_sat = sqrt(thing);
thing_2 = (3/2)*(J2_const_mars/(theta_dot*a_sat^(7/2)))*cos(pi-i2_sat)+1;
e2_sat_2 = sqrt(thing_2);

%% 11-21-24
clear; clc;
mu_Mars = 4.3e4;
r_Mars = 3400;

a_target = r_Mars + 200;
n = sqrt(mu_Mars/a_target^3);





