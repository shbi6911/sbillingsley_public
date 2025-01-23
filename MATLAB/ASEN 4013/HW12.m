%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
%givens
Fi = 3000000;   gamma = 1.152;      P_c = 700*144;
Cstar_i = 5912;     Cstar_x = 0.95*Cstar_i;
epsilon = 7.72;     P_e = 2116.2;       P_ratio = P_e/P_c;
g_c = 32.174;   %ft*lbm/lbf*s^2
%find thrust coefficient
Gamma = sqrt(gamma/((gamma+1)/2)^((gamma+1)/(gamma-1)));
thing1 = 1 - P_ratio^((gamma-1)/gamma);
thing2 = sqrt(((2*gamma)/(gamma-1))*thing1);
CF_i = Gamma*thing2;
%find throat area
A_t = Fi/(CF_i*P_c);
D_t = sqrt((4*A_t)/pi);
%find mass flow rate
m_dot = (g_c*P_c*A_t)/Cstar_x;