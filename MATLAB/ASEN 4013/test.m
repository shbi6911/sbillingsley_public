%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

clear; clc;
%givens
delta = 0.02765;    P_std = 2116.2;     P_0 = P_std*delta;  %lbf/ft^2
omega = 0.7668;     T_std = 518.69;     T_0 = T_std*omega;  %Rankine
rho_std = 0.07647;  %lbm/ft^3
a_std = 1116;       %ft/s
rho_0 = rho_std*(delta/omega);      a_0 = a_std*sqrt(omega);
Isp = 400;      gamma = 1.26;       R = 173.6;      %ft*lbf/lbm*R
T_c = 4840;     A_t = 0.5;  %ft^2
F_i = 100000;   %lbf
g_c = 32.174;   %ft*lbm/lbf*s^2
c_p = (gamma*R)/(gamma-1);

m_dot = F_i/Isp;
Gamma = sqrt(gamma/((gamma+1)/2)^((gamma+1)/(gamma-1)));
P_c = (m_dot*sqrt(T_c)*sqrt(R/g_c))/(A_t*Gamma);
P_e = 4.199*144;
P_ratio = P_e/P_c;
P_a = 14.7*144;

%find effective nozzle exit ratio
stuff = ((P_ratio)^(2/gamma)) - (P_ratio^((gamma+1)/gamma));
epsilon = Gamma/sqrt(((2*gamma)/(gamma-1))*stuff);
%find ideal thrust coefficient 
thing1 = 1 - P_ratio^((gamma-1)/gamma);
thing2 = sqrt(((2*gamma)/(gamma-1))*thing1);
CF_i = Gamma*thing2 + ((P_ratio - (P_a/P_c))*epsilon);

F_i_sl = CF_i*P_c*A_t;