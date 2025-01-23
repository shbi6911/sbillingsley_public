%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
clear; clc;
%givens
M1 = 0.8;   h = 40000; Tt4 = 2520; %degRankine
alpha = 5;  pi_f = 2;   pi_c = 20;  gamma = 1.4;    c_p = 0.24; %BTU/lbm*R
h_PR = 18400;   %BTU/lbm
R = 53.353;     %ft*lbf/lbm (deg R)
gc = 32.174;
delta = 0.1858; P_0 = delta*2116.2;    %lb/ft^2
theta = 0.7519; T_0 = theta*518.69;    %rankine
a_0 = sqrt(theta)*1116;                %ft/s
rho = (delta/theta)*0.07647;            %lbm/ft^3

%turbojet
%Assumptions
tau_d = 1;  tau_n = 1;  pi_d = 1;   pi_n = 1;   

tau_c = pi_c^((gamma-1)/gamma);

%freestream
tau_r = 1 + (((gamma-1)/2)*M1^2);   pi_r = (1 + (((gamma-1)/2)*M1^2))^(gamma/(gamma-1));

tau_lambda = Tt4/T_0;

%temp ratios of the turbines
tau_t = 1 - (tau_r/tau_lambda)*(tau_c -1);   

%fuel-air ratio
f = ((c_p*T_0)/h_PR)*(tau_lambda - tau_r*tau_c);

%exit Mach
M2 = sqrt((2/(gamma-1))*(tau_r*tau_c*tau_t - 1));

%exit velocity
V9 = a_0*M2*sqrt(tau_lambda/(tau_r*tau_c));

%specific thrust
F_m0 = V9 - M1*a_0;

%thrust specific fuel consumption
S = f/F_m0;

%Turbofan
tau_f = pi_f^((gamma-1)/gamma);

M19 = sqrt((2/(gamma-1))*(tau_r*tau_f - 1));

tau_t_2 = 1 - (tau_r/tau_lambda)*(tau_c -1 + (alpha*(tau_f -1)));

coeff = tau_lambda - tau_r*(tau_c -1 + (alpha*(tau_f -1)));
V9_2 = a_0*sqrt(((2/(gamma-1))*coeff) - (tau_lambda/(tau_r*tau_c)));
V19 = a_0*M19;

coeff2 = (V9/a_0) - M1 + alpha*((V19/a_0) - M1);
F_m0_2 = (a_0)*(1/(1+alpha)*coeff2);

S2 = f/((1+alpha)*F_m0_2);


%% Problem 2
clear; clc;
stages = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
ratio = zeros(1,length(stages));
compressor_ratio = ones(1,length(stages));
for i = 1:length(stages)
    ratio(i) = (1+(50/(288.15 + (50*(i-1)))))^((1.4*0.9)/(0.4));
    compressor_ratio(i) = prod(ratio(1:i));
end

