%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1

%givens
const.mu.Mars = 4.305e4;    %grav parameter in km^3/s^2
const.radius.Mars = 3397.2; %radius in km
h = 2.4799e4;               %momentum in km^2/s
r1 = 15000;                 %radius in km
v1 = 1.8111;                %velocity in km/s
r2 = 19000;                 %radius in km

%find orbital energy
energy1 = ((v1^2)/2) - (const.mu.Mars/r1);
%find semimajor axis
a = -const.mu.Mars/(2*energy1);

%find eccentricity
e = sqrt(1-((h^2)/(const.mu.Mars*a)));

%find true anomaly at time t1
theta_star1 = -acos((a*(1-e^2) - r1)/(r1*e));
disp(rad2deg(theta_star1));
%find eccentric anomaly at time t1
E1 = TrueToEccentric(theta_star1,e);
disp(rad2deg(E1));

%find true anomaly at time t2
theta_star2 = -acos((a*(1-e^2) - r2)/(r2*e));
disp(rad2deg(theta_star2));
%find eccentric anomaly at time t1
E2 = TrueToEccentric(theta_star2,e);
disp(rad2deg(E2));

%find mean motion
n = sqrt(const.mu.Mars/a^3);
%find elapsed time
E2_pos = (2*pi)+E2;
t2 = (1/n)*((E2_pos - e*sin(E2_pos)) - (E1 - e*sin(E1)));
t2_days = t2/(3600*24);

%orbital period
P = 2*pi*sqrt(a^3/const.mu.Mars);

%% Problem 1 plotting
% b = a*sqrt(1-e^2);
% t_vec = linspace(0,(2*pi),500);
% X = a*cos(t_vec);
% Y = b*sin(t_vec);
% plot(X,Y,'k','LineWidth',2); hold on; grid on;
% % Y2 = a_1*sin(t_vec);
% % plot(X,Y2,'r');
% plot(0,0,'.','MarkerSize',15,'Color','k');
% plot(a*e,0,'.','MarkerSize',15,'Color','k');
% % X3 = linspace(-a_1,a_1,1000);
% xlim([-(1.1*a),(1.1*a)]);
% ylim([-(1.1*a),(1.1*a)]);
% axis equal

%% Problem 2
clear; clc;
const.mu.Mars = 4.305e4;    %grav parameter in km^3/s^2
const.radius.Mars = 3397.2; %radius in km
%given Cartesian vectors
R3 = [-7.665e3;6.5468e3;-4.574e2];
V3 = [1.6334;0.1226;-1.9455];

%find Keplerian elements
[a3,e3,i3,RAAN3,omega3,theta_star3] = CartesianToKepler(R3,V3,const.mu.Mars);

%find eccentric anomaly at time t3
E3 = TrueToEccentric(theta_star3,e3);
disp(rad2deg(E3));

%orbital period
P = 2*pi*sqrt(a3^3/const.mu.Mars);

%find time t3
t3 = sqrt(a3^3/const.mu.Mars)*(E3 - (e3*sin(E3)));
%add two hours to get t4
t4 = t3 + 7200;
%find eccentric anomaly at t4
E4 = EccentricSolver(a3,const.mu.Mars,e3,t4,10^-8);
%find true anomaly at t4
theta_star4 = EccentricToTrue(E4,e3);

%find radius and velocity vectors at t4 using f and g functions
delta_theta = theta_star4 - theta_star3;
p3 = a3*(1-e3^2);
r4 = p3/(1+e3*cos(theta_star4));
r3 = norm(R3);
f = 1-(r4/p3)*(1-cos(delta_theta));
g = (r4*norm(R3)*sin(delta_theta))/sqrt(const.mu.Mars*p3);
coeff = ((1-cos(delta_theta))/p3) - (1/r4) - (1/norm(R3)); 
f_dot = sqrt(const.mu.Mars/p3)*tan(delta_theta/2)*coeff;
g_dot = 1 - (norm(R3)/p3)*(1-cos(delta_theta));

R4 = f*R3 + g*V3;
V4 = f_dot*R3 + g_dot*V3;

%% Problem 3
%% Problem 2
clear; clc;
const.mu.Mars = 4.305e4;    %grav parameter in km^3/s^2
const.radius.Mars = 3397.2; %radius in km
%given unit vectors
nd_hat = [-0.64279;-0.76604;0];
ra_hat = [-0.02970;-0.97508;-0.21985];

%get line of nodes and eccentricity unit vectors
n_hat5 = -nd_hat;
e_hat5 = -ra_hat;

%RAAN
RAAN = acos(n_hat5(1));
if n_hat5(2) < 0
    RAAN = -RAAN;
end
%argument of periapsis
omega = acos(dot(n_hat5,e_hat5));
if e_hat5(3) < 0
    omega = -omega;
end