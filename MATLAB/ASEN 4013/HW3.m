%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 3
clear; clc;
%preallocate for final values
final = zeros(5,6);

h = 50000;      %altitude in feet
A_1 = 4.235;    %area of inlet in ft^2
M_1 = 3;        final(2,1) = M_1;   %freestream Mach number
psia = 0.006944;    %conversion factor from lb/ft^2 to psi
gamma = 1.4;    %ratio of specific heats

%Mach numbers at subsequent stations
M_2 = 1;        final(2,2) = M_2;
M_3 = 0.15;     final(2,3) = M_3;
M_4 = 0.15;     final(2,4) = M_4;
M_5 = 1;        final(2,5) = M_5;
M_6 = 3;        final(2,6) = M_6;

%standard sea-level atmospheric values
P_std = 2116;       %lbf/ft^2
T_std = 518.69;     %deg Rankine
rho_std = 0.07647;  %lbm/ft^3
a_std = 1116;       %ft/s

%table values at 50000 ft
delta_1 = 0.1151;   theta_1 = 0.7519;
%atmospheric values at altitude
T_1 = theta_1*T_std;    P_1 = delta_1*P_std;    
rho_1 = rho_std*(delta_1/theta_1);
a_1 = a_std*sqrt(theta_1);

%find mass flow rate
V_1 = M_1*a_1;      m_dot_1 = rho_1*A_1*V_1;
disp("mass flow rate is " + string(m_dot_1) + " lbm/s");

%fill in table station 1
final(1,1) = A_1;   final(3,1) = P_1*psia;
final(4,1) = T_1;   final(5,1) = V_1;

%find total pressure and temperature of the freestream (stations 1-3)
thing_1 = 1 + ((gamma-1)/2)*M_1^2;
thing_2 = 1 + ((gamma-1)/2)*M_2^2;
thing_3 = 1 + ((gamma-1)/2)*M_3^2;
P_exp = -gamma/(gamma-1);
T_exp = -1;
rho_exp = -1/(gamma-1);
T_tot_123 = T_1/(thing_1^T_exp);
P_tot_123 = P_1/(thing_1^P_exp);
rho_tot_123 = rho_1/(thing_1^rho_exp);

%find temp, pressure, density at stations 2 and 3
T_2 = T_tot_123*thing_2^T_exp;
P_2 = P_tot_123*thing_2^P_exp;
rho_2 = rho_tot_123*thing_2^rho_exp;
T_3 = T_tot_123*thing_3^T_exp;
P_3 = P_tot_123*thing_3^P_exp;
rho_3 = rho_tot_123*thing_3^rho_exp;

%find areas using continuity
V_2 = M_2*sqrt(1.4*32.174*53.34*T_2);
V_3 = M_3*sqrt(1.4*32.174*53.34*T_3);
A_2 = m_dot_1/(rho_2*V_2);
A_3 = m_dot_1/(rho_3*V_3);

%find areas using A/A*
A_exp = (gamma+1)/(2*(gamma-1));
A_ratio_1 = (1/M_1)*((2/(gamma+1))*thing_1)^A_exp;
A_ratio_2 = (1/M_2)*((2/(gamma+1))*thing_2)^A_exp;
A_ratio_3 = (1/M_3)*((2/(gamma+1))*thing_3)^A_exp;
A_star = A_1/A_ratio_1;
A_2_check = A_star*A_ratio_2;
A_3_check = A_star*A_ratio_3;

%fill in table stations 2 and 3
final(1,2) = A_2;   final(3,2) = P_2*psia;
final(4,2) = T_2;   final(5,2) = V_2;

final(1,3) = A_3;   final(3,3) = P_3*psia;
final(4,3) = T_3;   final(5,3) = V_3;

%assume constant total pressure at all stations
P_tot_456 = P_tot_123;

%find temperatures after combustor
T_tot_456 = 4000;
T_4 = T_tot_456*thing_3^T_exp;
T_5 = T_tot_456*thing_2^T_exp;
T_6 = T_tot_456*thing_1^T_exp;

%find pressures and densities after combustor
P_4 = P_tot_456*thing_3^P_exp;
P_5 = P_tot_456*thing_2^P_exp;
P_6 = P_tot_456*thing_1^P_exp;

rho_4 = P_4/(53.34*T_4);
rho_5 = P_5/(53.34*T_5);
rho_6 = P_6/(53.34*T_6);

%find velocities after combustor
V_4 = M_4*sqrt(1.4*32.174*53.34*T_4);
V_5 = M_5*sqrt(1.4*32.174*53.34*T_5);
V_6 = M_6*sqrt(1.4*32.174*53.34*T_6);

%find areas after combustor
A_4 = m_dot_1/(rho_4*V_4);
A_5 = m_dot_1/(rho_5*V_5);
A_6 = m_dot_1/(rho_6*V_6);

%fill in table stations 4, 5 and 6
final(1,4) = A_4;   final(3,4) = P_4*psia;
final(4,4) = T_4;   final(5,4) = V_4;

final(1,5) = A_5;   final(3,5) = P_5*psia;
final(4,5) = T_5;   final(5,5) = V_5;

final(1,6) = A_6;   final(3,6) = P_6*psia;
final(4,6) = T_6;   final(5,6) = V_6;

%find thrust forces
T_diff = (m_dot_1/32.174)*(V_3-V_1) + (P_1 - P_1)*A_1 - (P_3 - P_1)*A_3;
T_comb = (m_dot_1/32.174)*(V_4-V_3) + (P_3 - P_1)*A_3 - (P_4 - P_1)*A_4;
T_nozz = (m_dot_1/32.174)*(V_6-V_4) + (P_4 - P_1)*A_4 - (P_6 - P_1)*A_6;
T_engine = (m_dot_1/32.174)*(V_6-V_1) + (P_1 - P_1)*A_1 - (P_6 - P_1)*A_6;

