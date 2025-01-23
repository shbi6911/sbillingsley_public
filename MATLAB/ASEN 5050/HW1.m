%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

% ASEN 5050 Homework 1 Due 9/10/24\
% Shane Billingsley

const.mu.Sun = 1.32712428e11;  %grav parameters in km^3/s^2
const.mu.Jupiter = 1.268e8;

%heliocentric s/c orbital parameters
r_1 = 1.6599e8;   %orbital radius in km
a_1 = 2.3132e8;   %semimajor axis in km
v_r_1 = -11.6485;   %radial component of velocity in km/s

energy_1 = -const.mu.Sun/(2*a_1);
v_1 = sqrt(2*(energy_1+(const.mu.Sun/r_1)));
gamma_1 = asin(v_r_1/v_1);
v_theta_1 = v_1*cos(gamma_1);
h_1 = r_1*v_1*cos(gamma_1);
p_1 = h_1^2/const.mu.Sun;
e_1 = norm(cross([v_r_1;v_theta_1;0],[0;0;h_1])/const.mu.Sun - [1;0;0]);
theta_star_1 = acos(((h_1^2/const.mu.Sun)-r_1)/(r_1*e_1));
P_1 = 2*pi*sqrt(a_1^3/const.mu.Sun);
P_1_hrs = P_1/3600;
P_1_days = P_1_hrs/24;

%inertial orbital vectors
R_X_1 = 1.0751e8;       %inertial orbital radius in km
R_Y_1 = -1.2647e8;       
R_Z_1 = 1.3644e5;
R_1_vec = [R_X_1;R_Y_1;R_Z_1];

V_X_1 = 1.5180e1;         %inertial velocity in km/s
V_Y_1 = 2.8193e1;
V_Z_1 = 1.0504e-2;
V_1_vec = [V_X_1;V_Y_1;V_Z_1];

h_inert_1_vec = cross(R_1_vec,V_1_vec);
e_inert_1 = (cross(V_1_vec,h_inert_1_vec)/const.mu.Sun)-R_1_vec./norm(R_1_vec);
energy_inert_1 = (0.5*dot(V_1_vec,V_1_vec))-(const.mu.Sun/norm(R_1_vec));

%% plotting
% b_1 = a_1*sqrt(1-e_1^2);
% t_vec = linspace(0,(2*pi),500);
% X = a_1*cos(t_vec);
% Y = b_1*sin(t_vec);
% plot(X,Y,'k','LineWidth',2); hold on; grid on;
% % Y2 = a_1*sin(t_vec);
% % plot(X,Y2,'r');
% plot(0,0,'.','MarkerSize',15,'Color','k');
% plot(a_1*e_1,0,'.','MarkerSize',15,'Color','k');
% % X3 = linspace(-a_1,a_1,1000);
% xlim([-(1.1*a_1),(1.1*a_1)]);
% ylim([-(1.1*a_1),(1.1*a_1)]);
% axis equal



%% Problem 2

v_infty_2 = 10.7527;    %hyperbolic velocity at infinity in km/s
theta_star_infty_2_deg = 139.3724;  %true anomaly at infinity in degrees
theta_star_infty_2 = deg2rad(theta_star_infty_2_deg);%convert to radians

a_2 = -const.mu.Jupiter/(v_infty_2^2);
e_2 = -1/(cos(theta_star_infty_2));
delta_2 = 2*asin(1/e_2);

r_p_2 = a_2*(1-e_2);
v_p_2 = sqrt(2*(const.mu.Jupiter/r_p_2)-(const.mu.Jupiter/(a_2)));

%energy_2 = (v_infty_2^2)/2;

%energy_2_p = (v_p_2^2)/2 - const.mu.Jupiter/r_p_2;
