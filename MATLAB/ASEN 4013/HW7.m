%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
clear; clc;
% givens
M_0 = 0.8;      P_0 = 29920;    %Pa
T_0 = 229;      Tt_4 = 1780;    %K
pi_dmax = 0.95;     pi_c = 9;   pi_b = 0.94;    pi_n = 0.98;
e_c = 0.85;     e_t = 0.88;     h_PR = 42800000;    %J/kg
gamma_c = 1.4;      gamma_t = 1.3;  eta_b = 0.96;   eta_m = 0.98;
cp_c = 1004;    cp_t = 1235;    %J/kg*K
P_09 = 0.8;     P_9 = (1/P_09)*P_0;

%find ram pressure and temperature
[~,tau_r,pi_r,~,~] = flowisentropic(gamma_c,M_0);
tau_r = 1/tau_r;    pi_r = 1/pi_r;

%find compressor exit total pressure and temperature
Pt_3 = pi_c*pi_dmax*pi_r*P_0;
tau_c = pi_c^((gamma_c-1)/(gamma_c*e_c));
Tt_3 = tau_c*tau_r*T_0;

%assume ht3 = cp_c*Tt3
ht3 = cp_c*Tt_3;
%find fuel/air ratio assuming ht4 = cp_t*Tt4
f = ((cp_t*Tt_4) - ht3)/((eta_b*h_PR) - (cp_t*Tt_4));

%find turbine temperature and pressure ratio
tau_lambda = (cp_t*Tt_4)/(cp_c*T_0);    
tau_t = 1 - (1/(eta_m*(1+f)))*(tau_r/tau_lambda)*(tau_c - 1);
pi_t = tau_t^((gamma_t)/((gamma_t-1)*e_t));

%find Pt9/P9 and M9
%pi_dmax is used b/c we are subsonic thus eta_d = 1
Pt9_P9 = P_09*pi_r*pi_dmax*pi_c*pi_b*pi_t*pi_n;
M_9 = sqrt((2/(gamma_t-1))*(Pt9_P9^((gamma_t-1)/gamma_t)-1));

%find turbine exit total pressure and temperature
Tt_5 = Tt_4*tau_t;
Pt_5 = pi_t*pi_b*Pt_3;


%% Problem 2
clear; clc;
%givens
gc = 32.174;        %(ft*lbm)/(lbf*s^2)
M_6 = 0.4;          mdot_6 = 230;       %lbm/s
Tt_6 = 1830;        Tt_8 = 3660;    %Rankine
Pt_6 = 38*144;      %lbf/ft^2
gamma_6 = 1.33;     gamma_8 = 1.3;
cp_6 = 0.276;       cp_8 = 0.297;       %BTU/(lbm*R)
R = 53.34;      %(ft*lbf)/(lbm*R)
Pt_86_off = 0.97;       Pt_86_on = 0.94;
eta_AB = 0.95;      h_PR = 18400;   %BTU/lbm

%get static/total temp and pressure ratios
[~,T_ratio_6,P_ratio_6,~,AAstar_off] = flowisentropic(gamma_6,M_6);
%find static temperature and pressure at station 6
P_6 = Pt_6*P_ratio_6;       T_6 = Tt_6*T_ratio_6;
%find area at station 6
A_6 = (mdot_6/(gc*P_6*M_6))*sqrt((gc*R*T_6)/gamma_6);
%find area of the throat with no afterburner
AAstar_real = AAstar_off*Pt_86_off;
A_8_off = A_6/AAstar_real;

%find fuel flow rate
ht6 = cp_6*Tt_6;
ht8 = cp_8*Tt_8;
f = (ht8 - ht6)/((eta_AB*h_PR) - ht8);
mdot_AB = mdot_6*f;
mdot_8 = mdot_6 + mdot_AB;

%find area with afterburner on
M_8 = 1;
[~,T_ratio_8,P_ratio_8,~,~] = flowisentropic(gamma_8,M_8);
P_8 = Pt_6*Pt_86_on*P_ratio_8;
T_8 = Tt_8*T_ratio_8;
A_8_on = ((mdot_8)/(gc*P_8*M_8))*sqrt((gc*R*T_8)/gamma_8);


%% Problem 3
clear; clc;

M_0 = 0.8;  T_0 = 390;  e_f = 0.89;    pi_f = 1.65;     gc = 32.174;
pi_dmax = 0.99;     pi_c = 36;         e_c = 0.9;   BTU = 778.16;
gamma_c = 1.4;      cp_c = 0.24;       alpha = 10;
pi_b = 0.96;        eta_b = 0.99;      h_PR = 18400;
Tt_4 = 3000;        P_09 = 0.9;        eta_m = 0.98;
e_t = 0.89;         gamma_t = 1.33;    cp_t = 0.276;
P_019 = 0.9;        pi_fn = 0.99;      pi_n = 0.99;

%find ram ratios and static temperature and R-values
R_c = ((gamma_c - 1)/gamma_c)*(cp_c*BTU);
R_t = ((gamma_t - 1)/gamma_t)*(cp_t*BTU);
a_0 = sqrt(gamma_c*R_c*T_0*gc);
V_0 = a_0*M_0;
tau_r = 1 + (((gamma_c - 1)/2)*M_0^2);
pi_r = tau_r^(gamma_c/(gamma_c-1));
%find compressor values
tau_c = pi_c^((gamma_c -1)/(gamma_c*e_c));
%eta_c = ((pi_c^((gamma_c - 1)/gamma_c)) - 1)/(tau_c - 1);
%find fan values
tau_f = pi_f^((gamma_c -1)/(gamma_c*e_f));
%eta_f = ((pi_f^((gamma_c - 1)/gamma_c)) - 1)/(tau_f - 1);
%find fuel flow
Tt_3 = T_0*tau_r*tau_c;
ht_3 = cp_c*Tt_3;   
ht_4 = cp_t*Tt_4;
f = (ht_4 - ht_3)/(h_PR*eta_b - ht_4);
%find turbine values
tau_lambda = (cp_t*Tt_4)/(cp_c*T_0);
tau_t = 1 - ((1/(eta_m*(1+f)))*(tau_r/tau_lambda)*(tau_c -1 + (alpha*(tau_f - 1))));
pi_t = tau_t^(gamma_t/((gamma_t -1)*e_t));
%eta_t = (1 - tau_t)/(1 - tau_t^(1/e_t));

%find exit Mach number & velocity
Pt9_P9 = P_09*pi_r*pi_dmax*pi_c*pi_b*pi_t*pi_n;
Tt9_T9 = (cp_c/cp_t)*tau_lambda*tau_t;
T_90 = Tt9_T9/(Pt9_P9^((gamma_t -1)/gamma_t));
T_9 = T_90*T_0;
M_9 = sqrt((2/(gamma_t-1))*(Pt9_P9^((gamma_t -1)/gamma_t) - 1));
V9_a0 = sqrt(((gamma_t*R_t*T_9)/(gamma_c*R_c*T_0))*M_9^2);

%find exit Mach number & velocity for the fan
Tt19_T0 = tau_r*tau_f;
Pt19_P19 = P_019*pi_r*pi_dmax*pi_f*pi_fn;
T19_T0 = (Tt19_T0)/(Pt19_P19^((gamma_c -1)/gamma_c));
M_19 = sqrt((2/(gamma_c-1))*(Pt19_P19^((gamma_c -1)/gamma_c) - 1));
V19_a0 = sqrt(T19_T0*M_19^2);

%find specific thrust using Eq 7.48 and TSFC using Eq 7.49
fan = (alpha/(1+alpha))*(a_0/gc)*(V19_a0 - M_0 + (T19_T0/V19_a0)*((1-P_019)/gamma_c));
jet = (1/(1+alpha))*(a_0/gc)*(((1+f)*(V9_a0)) - M_0 + (1+f)*((R_t*T_90)/(R_c*V9_a0))*((1-P_09)/gamma_c));
F_m0 = jet + fan;
S = f/((1+alpha)*F_m0);
S_hr = S*3600;
