%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
clear; clc;
%define constants
gamma = 1.4;    theta = deg2rad(7);     M1_1 = 1.6;     M1_2 = 1.5;
%Problem 1a (normal shock, Mach 1.6 incoming)
[~, ~ , ~ , ~ , M2_1a, Pt_ratio_1a, ~] = ...
    flownormalshock(gamma, M1_1);
%Problem 1b (oblique shocks, Mach 1.6 incoming)

%interpolate beta angle between 6 and 8 deg
Beta_1_1 = deg2rad((45.34+48.03)/2);     
%find normal component of Mach
Mn1_1 = M1_1*sin(Beta_1_1);
%find downstream properties based on normal Mach component
[~, ~ , ~ , ~ , Mn2_1b, Pt_ratio_1b, ~] = ...
    flownormalshock(gamma, Mn1_1);
%find total Mach number downstream of first shock
M2_1 = Mn2_1b/sin(Beta_1_1-theta);
%find new Beta angle for the second ramp.  Theta is the same and use M2_1
Beta_2_1 = (deg2rad((58.23+66.91)/2) + deg2rad((57.43+64.29)/2))/2;
%find normal component of Mach
Mn2_1 = M2_1*sin(Beta_2_1);
%find downstream properties based on normal Mach component
[~, ~ , ~ , ~ , Mn3_1c, Pt_ratio_1c, ~] = ...
    flownormalshock(gamma, Mn2_1);
%find total Mach number downstream of second shock
M3_1 = Mn3_1c/sin(Beta_2_1-theta);
% flow is still barely supersonic, so we go through a weak normal shock
[~, ~ , ~ , ~ , M4_1, Pt_ratio_1d, ~] = ...
    flownormalshock(gamma, M3_1);
% find total pressure loss by combining ratios
Pt_ratio_1tot = Pt_ratio_1b*Pt_ratio_1c*Pt_ratio_1d;

%output values for reference
disp("@M1.6, normal shock inlet produces M"+ string(M2_1a) + ...
    " downstream, with pressure ratio " +string(Pt_ratio_1a));
disp("@M1.6, oblique shock inlet produces M"+ string(M4_1) + ...
    " downstream, with pressure ratio " +string(Pt_ratio_1tot));

%Problem 1c, reevaluate ramp inlet for incoming Mach 1.5
%interpolate beta angle between 6 and 8 deg
Beta_1_2 = deg2rad((49.33+52.57)/2);
%find normal component of Mach
Mn1_2 = M1_2*sin(Beta_1_2);
%find downstream properties based on normal Mach component
[~, ~ , ~ , ~ , Mn2_2b, Pt_ratio_2b, ~] = ...
    flownormalshock(gamma, Mn1_2);
%find total Mach number downstream of first shock
M2_2 = Mn2_2b/sin(Beta_1_2-theta);
%no solution for an oblique shock at 7 deg and M1.25.  Use normal shock
[~, ~ , ~ , ~ , M3_2, Pt_ratio_2c, ~] = ...
    flownormalshock(gamma, M2_2);
%find total pressure loss
Pt_ratio_2tot = Pt_ratio_2b*Pt_ratio_2c;
%display output
disp("@M1.5, oblique shock inlet produces M"+ string(M3_2) + ...
    " downstream, with pressure ratio " +string(Pt_ratio_2tot));

%% Problem 2
clear; clc;
M1 = 0.6;   P1 = 10;    T1 = 500;   D=1;    A=1^2;  gamma = 1.4; cf = 0.004; %givens
L = 8;
%upstream properties
[~, T_ratio_1, P_ratio_1, ~, ~, Ptot_ratio_1, fanno_1] = flowfanno(gamma, M1);
[M2,~,~,~,~,~,~] = flowfanno(gamma, (fanno_1 -((4*cf*L)/D)), 'fannosub'); %find M2
%downstream properties
[~, T_ratio_2, P_ratio_2, ~, ~, Ptot_ratio_2, fanno_2] = flowfanno(gamma, M2);
% find final ratios
T2 = T1*(T_ratio_2/T_ratio_1);      P2 = P1*(P_ratio_2/P_ratio_1);
Ptot_ratio_3 = Ptot_ratio_2/Ptot_ratio_1;

%output
disp("Exit Mach = " + string(M2) + " Exit Temp = " + string(T2) + " Exit Pressure " + ...
    string(P2) + " Pt2/Pt1 = " + string(Ptot_ratio_3));

%% Problem 4
clear; clc;
%Givens
pi_c1 = 10; pi_c2 = 20;     M1 = 2;     T_0 = 217;  gamma = 1.4;    cp = 1004;
h_PR = 42800000;
a_0 = sqrt(gamma*286*T_0);
T_t_4 = [1800 2000 2200 2400];

%Assumptions
tau_d = 1;  tau_n = 1;  pi_d = 1;   pi_n = 1;   

tau_c1 = pi_c1^((gamma-1)/gamma);   tau_c2 = pi_c2^((gamma-1)/gamma);

%freestream
tau_r = 1 + (((gamma-1)/2)*M1^2);   pi_r = (1 + (((gamma-1)/2)*M1^2))^(gamma/(gamma-1));

tau_lambda = T_t_4./T_0;

%temp ratios of the turbines
tau_t1 = 1 - (tau_r./tau_lambda)*(tau_c1 -1);   
tau_t2 = 1 - (tau_r./tau_lambda)*(tau_c2 -1);

%fuel-air ratio
f1 = ((cp*T_0)/h_PR).*(tau_lambda - tau_r*tau_c1);  
f2 = ((cp*T_0)/h_PR).*(tau_lambda - tau_r*tau_c2);

%exit Mach
M2_1 = sqrt((2/(gamma-1)).*(tau_r.*tau_c1.*tau_t1 - 1));
M2_2 = sqrt((2/(gamma-1)).*(tau_r.*tau_c2.*tau_t2 - 1));

%exit velocity
V2_1 = a_0.*M2_1.*sqrt(tau_lambda./(tau_r.*tau_c1));
V2_2 = a_0.*M2_2.*sqrt(tau_lambda./(tau_r.*tau_c2));

%specific thrust
F_m0_1 = V2_1 - M1*a_0;
F_m0_2 = V2_2 - M1*a_0;

%thrust specific fuel consumption
S1 = f1./F_m0_1;
S2 = f2./F_m0_2;
