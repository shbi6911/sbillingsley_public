%By:        Shane Billingsley
%Class:     ASEN 3728 Aircraft Dynamics
%Date:      Spring 2024

%% Question 1
clear; close all;
%define constants
m = 6.77;   %mass in kg
c = 0.2;    %chord in meters
b = 3.3;    %span in meters
S = 0.67;   %planform area in m^2
Cd_min = 0.02;  %minimum coefficient of drag
K = 0.0224; %k constant (1/pi*e*AR)
Cl_alpha = 5.75;    %change in lift coefficient with angle of attack
Cl_q = 10.14;       %change in lift coefficient with pitch rate
Cl_delta_e = 0.0079;    %change in lift coefficient with change in elevator
Cm_q = -24.4;       %change in moment coefficient with pitch rate
Cm_delta_e = -0.02; %change in moment coefficient with change in elevator
Cm_0 = 0.12;    %initial moment coefficient
Cl_min = 0; %minimum lift coefficient

%find drag coefficient under given conditions
%define environmental variables
rho = 1.10; %air density in kg/m^3
V_a = 21;   %airspeed in m/s
g = 9.81;   %acceleration due to gravity in m/s^2

C_l = (2*m*g)/(rho*V_a^2);
C_D = Cd_min + (K*C_l^2);
disp("Coefficient of drag = " + string(C_D));

%under same conditions, and given alpha and elevator, move CG to trim
alpha = 4.07;   alpha = deg2rad(alpha); %angle of attack in degrees
delta_e = 1.65; delta_e = deg2rad(delta_e); %angle of attack in radians
h_n = 0.75; %location of neutral point

moment = -Cm_delta_e*delta_e - Cm_0;
lift = Cl_alpha*alpha + Cl_delta_e*delta_e;
h = moment/lift +h_n;
disp("Center of gravity is h = "+string(h));

%% Question 2
clear; close all;
%define constants
m = 22000;  %mass in kg
S = 87;     %planform area in m^2
c=3.9;      %chord in meters
b=22.3;     %span in meters
a_wb=4.4;   %lift slope in 1/rad
K_n=0.25;     %static margin
S_t = 21.4; %planform area of tail in m^2
de_da = 0.3;    %change in epsilon with change in alpha
a_t = 3.67; %lift slope of tail in 1/rad
e = 0.95;   %efficiency factor
rho = 1.225;    %air density in kg/m^3
alpha_0 = 1.2;  alpha_0 = deg2rad(alpha_0); %trim angle of attack in rads
u_0 = 125;  %trim forward velocity in m/s
gamma_a_0 = 0;  %yaw angle at trim
T_0 = 14000;    %thrust at trim
g = 9.81;   %acceleration due to gravity in m/s^2
AR = b^2/S; K = 1/(pi*e*AR);    %K constant using givens
lift = (m*g) - T_0*sin(alpha_0);    %find lift using force balance
drag = T_0*cos(alpha_0);            %find drag using force balance
V_a = u_0/cos(alpha_0);
q_infty = 0.5*rho*V_a^2;        %find dynamic pressure
C_l = lift/(q_infty*S);         %find coefficient of lift
C_d = drag/(q_infty*S);         %find coefficient of drag
disp("Coefficient of lift = "+string(C_l));
disp("Coefficient of drag = "+string(C_d));

C_l_alpha = a_wb*(1+((a_t*S_t)/(a_wb*S))*(1-de_da));
%C_d_0 = C_d - K*C_l^2;
C_z_alpha = -(C_l_alpha + C_d);
Z_w = 0.5*rho*u_0*S*C_z_alpha;
disp("Z_w = "+string(Z_w));

C_m_alpha = -C_l_alpha*(K_n);
M_w = 0.5*rho*u_0*c*S*C_m_alpha;
disp("M_w = "+string(M_w));

C_x_alpha = C_l*(1-2*K*C_l_alpha);
X_w = 0.5*rho*u_0*S*C_x_alpha;
disp("X_w = "+string(X_w));




