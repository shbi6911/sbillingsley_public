%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 2
clear; clc;
M_e = 5.3;  gamma = 1.3;    R = 317; P_c = 225*101325;  T_c = 2800;
D_e = 2.15;     A_e = (pi*D_e^2)/4;     epsilon = 62;
[~,T_ratio,P_ratio,~,A_ratio] = flowisentropic(1.3,5.3);

A_t = A_e/epsilon;
P_e = P_c*P_ratio;
T_e = T_c*T_ratio;

V_e = M_e*sqrt(R*gamma*T_e);

Gamma = sqrt(gamma/(((gamma+1)/2)^((gamma+1)/(gamma-1))));

m_dot = (Gamma*P_c*A_t)/(sqrt(R*T_c));

F_sl = (V_e*m_dot) + (P_e - 101325)*A_e;
F_alt = (V_e*m_dot) + 0;
F_vac = (V_e*m_dot) + (P_e - 0)*A_e;
