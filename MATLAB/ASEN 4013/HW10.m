%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1

%% Problem 2

C_star = 1205.2;        %m/s
gamma = 1.3531;
T_c = 859.12;   %K
M = 10.796;
R = 8.314;
P_c = 20*101325;
Isp = 230;
g0 = 9.80665;
C = Isp*g0;
m_dot = (0.05/C)/0.05;
A_t = (C_star*m_dot)/P_c;


P_e = 0.0087*101325;
[M_e,T_ratio,P_ratio,~,A_ratio] = flowisentropic(gamma,(P_e/P_c),'pres');
T_e = T_c*T_ratio;
A_e = A_t*A_ratio;
V_e = M_e*sqrt(gamma*T_e*(R/M));
C_check = V_e + ((P_e*A_e)/m_dot);





% val = 1;
% while val < 1e-1
