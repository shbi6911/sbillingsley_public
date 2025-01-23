%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1

clear; clc;
const.mu.Saturn = 3.794e7; %grav parameter in km^3/s^2
const.radius.Saturn = 60268; %radius of Saturn in km
Rvec_1 = [-720000; 670000; 310000];  %orbital radius in km, inertial
Vvec_1 = [2.160; -3.360; 0.620]; %velocity in km/s
R_1 = norm(Rvec_1); V_1 = norm(Vvec_1);     %magnitudes

%find Keplerian orbital elements and eccentric anomaly
[a_1,e_1,i_1,RAAN_1,omega_1,theta_star_1] = CartesianToKepler(Rvec_1,Vvec_1,const.mu.Saturn);
E_1 = TrueToEccentric(theta_star_1,e_1);

%find true anomaly and eccentric anomaly at impact
%impact point approximated as where orbital radius = planet radius
R_1_2 = const.radius.Saturn;
%find true anomaly at this point (negative b/c moving toward periapsis)
theta_star_1_2 = -acos((a_1*(1-e_1^2) - R_1_2)/(R_1_2*e_1));
%find eccentric anomaly at this true anomaly
E_1_2 = TrueToEccentric(theta_star_1_2,e_1);

%find mean motion
n_1 = sqrt(const.mu.Saturn/a_1^3);
%find elapsed time
E_1_pos = (2*pi)+E_1;
E_1_2_pos = (2*pi) + E_1_2;
t1 = (1/n_1)*((E_1_2_pos - e_1*sin(E_1_2_pos)) - (E_1_pos - e_1*sin(E_1_pos)));
t1_days = t1/(3600*24);

%find radius and velocity vectors at impact using f and g functions
delta_theta = theta_star_1_2 - theta_star_1;
p_1 = a_1*(1-e_1^2);
f = 1-(R_1_2/p_1)*(1-cos(delta_theta));
g = (R_1_2*R_1*sin(delta_theta))/sqrt(const.mu.Saturn*p_1);
coeff = ((1-cos(delta_theta))/p_1) - (1/R_1_2) - (1/R_1); 
f_dot = sqrt(const.mu.Saturn/p_1)*tan(delta_theta/2)*coeff;
g_dot = 1 - (R_1/p_1)*(1-cos(delta_theta));

Rvec_1_2 = f*Rvec_1 + g*Vvec_1;
Vvec_1_2 = f_dot*Rvec_1 + g_dot*Vvec_1;

%% Problem 2

clear; clc;
const.mu.Jupiter = 1.268e8; %grav parameter in km^3/s^2
const.radius.Jupiter = 71492;   %radius of Jupiter in km
Rvec_2 = [5.352950e6;7.053778e5;-4.059700e5];   %orbital radius in km, inertial
Vvec_2 = [-4.164248;1.963690;3.191257e-1];  %inertial velocity in km/s
R_2 = norm(Rvec_2); V_2 = norm(Vvec_2);     %magnitudes

%find Keplerian elements and eccentric anomaly
[a_2,e_2,i_2,RAAN_2,omega_2,theta_star_2] = CartesianToKepler(Rvec_2,Vvec_2,const.mu.Jupiter);
E_2 = TrueToEccentric(theta_star_2,e_2);

%find true anomaly at descending node as 180 deg - arg of periapsis
theta_star_2_2 = pi-omega_2;
%find eccentric anomaly at descending node
E_2_2 = TrueToEccentric(theta_star_2_2,e_2);

%find mean motion
n_2 = sqrt(const.mu.Jupiter/a_2^3);
%find elapsed time
t2 = (1/n_2)*((E_2_2 - e_2*sin(E_2_2)) - (E_2 - e_2*sin(E_2)));
t2_days = t2/(3600*24);

%find radius and velocity vectors at t2 using f and g functions
delta_theta = theta_star_2_2 - theta_star_2;
p_2 = a_2*(1-e_2^2);
R_2_2 = p_2/(1+e_2*cos(theta_star_2_2));
f = 1-(R_2_2/p_2)*(1-cos(delta_theta));
g = (R_2_2*R_2*sin(delta_theta))/sqrt(const.mu.Jupiter*p_2);
coeff = ((1-cos(delta_theta))/p_2) - (1/R_2_2) - (1/R_2); 
f_dot = sqrt(const.mu.Jupiter/p_2)*tan(delta_theta/2)*coeff;
g_dot = 1 - (R_2/p_2)*(1-cos(delta_theta));



Rvec_2_2 = f*Rvec_2 + g*Vvec_2;
Vvec_2_2 = f_dot*Rvec_2 + g_dot*Vvec_2;

%find eccentric and true anomalies at a time 20 days after t1
t1 = (1/n_2)*(E_2 - e_2*sin(E_2));
delta_t = (20*(3600*24)) + t1; %time after t1 in seconds
E_2_3 = EccentricSolver(a_2,const.mu.Jupiter,e_2,delta_t,10^-9);
theta_star_2_3 = EccentricToTrue(E_2_3,e_2);