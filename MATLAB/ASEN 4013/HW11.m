%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
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

%find mass flow rate at throat using thrust and specific impulse
m_dot = F_i/Isp;        %g_c and g0 values cancel in this eqn, needed for units
%find characteristic velocity using its definition
Gamma = sqrt(gamma/((gamma+1)/2)^((gamma+1)/(gamma-1)));
c_star = sqrt(R*g_c*T_c)/Gamma;
%find chamber pressure from C* and mass flow and compare
P_c1 = (c_star*(m_dot/g_c))/A_t;
P_c2 = (m_dot*sqrt(T_c)*sqrt(R/g_c))/(A_t*Gamma);
%find effective exhaust velocity at altitude
C = Isp*g_c;        %g_c has same value as g0
%find thrust coefficient and max ideal thrust coefficient
CF_imax = Gamma*sqrt((2*gamma)/(gamma-1));
CF = F_i/(P_c1*A_t);
%find expansion ratio
P_ratio = P_0/P_c1;      %exit to throat ratio assuming ideal expansion
thing = P_ratio^(2/gamma) - P_ratio^((gamma+1)/gamma);
epsilon = Gamma/sqrt((2*gamma)/(gamma-1)*thing);
%find exit diameter
A_e = epsilon*A_t;
D_e = sqrt(A_e*(4/pi));
%find ideal thrust coefficient assuming ideal expansion
thing1 = 1 - P_ratio^((gamma-1)/gamma);
thing2 = sqrt(((2*gamma)/(gamma-1))*thing1);
CF_i = Gamma*thing2;
%find separation point assuming ideal expansion
P_a = 3.5*P_0;
delta_1 = P_a/P_std;
alt = 53642;    %ft
%find molecular weight of exhaust products
R_u = 1545.349;
M = R_u/R;
%find velocity at the exit
%V_e1 = sqrt((2*gamma)/(gamma-1)*R*g_c*T_c*(1 - P_ratio^((gamma-1)/gamma)));
%F_sl1 = ((m_dot*V_e1)/g_c) + ((P_0 - P_std)*A_e);
coeffs = [1,((P_a*A_e*2*c_p)/m_dot),-(2*c_p*R*T_c)];
v_roots = roots(coeffs);
V_e = v_roots(sign(v_roots) > 0);
F_sl = ((m_dot*V_e)/g_c) + ((P_a - P_std)*A_e);
%CEARUN with:
%P_c = 57.85 atm
%liquid hydrogen @ 20K
%LOX @ 90K
%frozen flow at throat, infinite area
%find a f/o ratio of 2.412 gives gamma = 1.26 in chamber - good agreement
%with cstar
% 2.652 gives 1.26 gamma at throat with better agreement with molecular wt
% average gamma of 1.26 is at 3.52, M about 0.11 too high, cstar 25 m/s too high,
%chamber temp 38 K too high

%assuming equilibrium combustion and average gamma approach
%gamma = 1.26 for 3.5 ratio, errors comparable to previous


