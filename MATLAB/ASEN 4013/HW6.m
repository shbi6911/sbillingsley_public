%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

%% Problem 1
clear; clc;
%givens
e_c = 0.9;
T3 = [900 1300];
T2 = 288.15;
gamma = 1.4;

pi_c = (T3./T2).^((gamma*e_c)/(gamma-1));

%% Problem 2
clear; clc;
%givens
A_ct = 1.2;
%find Mach number at start
Mach = linspace(1,3,100000);    %Mach #s from 1-3
[~,~,~,~,area] = flowisentropic(1.4,Mach);  %A/A* at these Mach #s
[~,~,~,~,~,pressure,~] = flownormalshock(1.4,Mach); %Pty/Ptx at these Mach
area_ratios = area.*pressure;               %Plug into Eqn 11.7
[~,index] = min(abs(area_ratios - A_ct));   %find where it equals given A_ct
M_start = Mach(index);                      %pull starting Mach #
%find Mach number at the throat
area_throat = area(index)/A_ct;
[M_throat,~,~,~,~] = flowisentropic(1.4,area_throat,'sup');
%find Mach number which will unstart
[M_unstart,~,~,~,~] = flowisentropic(1.4,A_ct,'sup');
%find pressure recovery across normal shock at the throat
[~,~,~,~,~,pressure_tot,~] = flownormalshock(1.4,M_throat);

%% Problem 3
clear; clc;
%givens
m_dot = 75;       %kg/s
Pt8 = 350000;    Tt8 = 1600;
A_98 = 1.8;     gamma = 1.33;   R = 287;    P_98 = 0.98;
CD = 0.98;      P_0 = 40000;

[~,T_throat,P_throat,~,~] = flowisentropic(1.33,1);

%find area and radius of throat
A_8e = (m_dot/(Pt8*P_throat))*sqrt((R*Tt8*T_throat)/gamma);
A_8 = A_8e/CD;
r_8 = sqrt(A_8/pi);
%find area and radius of exit
A_9 = A_8*A_98;
r_9 = sqrt(A_9/pi);

%find ideal properties of exhaust
area_ratio_i = A_98/CD;         %ideal A/A*
%ideal exit Mach and P/Pt based on this A/A*
[M9i,~,pres_ratio_i,~,~] = flowisentropic(1.33,area_ratio_i,'sup');
P9i = pres_ratio_i*Pt8;         %ideal exit pressure
thing = ((2*gamma)/(gamma-1))*(1-(P9i/Pt8)^((gamma-1)/gamma));
V9i = sqrt(R*Tt8)*sqrt(thing);  %ideal exit velocity

%find real properties of exhaust
area_ratio = P_98*(A_98/CD);     %real A/A*
%ideal exit Mach and P/Pt based on this A/A*
[M9,~,pres_ratio,~,~] = flowisentropic(1.33,area_ratio,'sup');
P9 = pres_ratio*P_98*Pt8;         %real exit pressure
thing = ((2*gamma)/(gamma-1))*(1-(P9/Pt8)^((gamma-1)/gamma));
V9 = sqrt(R*Tt8)*sqrt(thing);  %real exit velocity

%find velocity coefficient
Cv = V9/V9i;

%find gross thrust coefficient
num1 = 1 - (P9i/Pt8)^(gamma-1)/(gamma);
denom1 = 1 - (P_0/Pt8)^(gamma-1)/(gamma);
coeff1 = sqrt(num1/denom1);
num2 = 1 - P_0/P9;
denom2 = (((Pt8*P_98)/P9)^(gamma-1)/(gamma))-1;
coeff2 = 1 + ((gamma-1)/(2*gamma))*(num2/denom2);
C_fg = CD*Cv*coeff1*coeff2;
%find actual gross thrust
F_g = (m_dot*V9) + (P9 - P_0)*A_9;
%display results
disp("Nozzle throat is " + string(A_8) + " m^2 with a radius of " +...
    string(r_8) + " meters");
disp("Nozzle exit is " + string(A_9) + " m^2 with a radius of " +...
    string(r_9) + " meters");
disp("Velocity coefficient is " + string(Cv));
disp("Gross thrust coefficient is " + string(C_fg));
disp("Actual thrust is " + string(F_g) + " N");


%% Problem 4
clear; clc;
%givens
T_ratio = 3;        CD = 2;     Mi = 0.05;
gamma_i = 1.38;     gamma_e = 1.3;

%find phi (whatever that is)
num1 = Mi^2*(1+ ((gamma_i-1)/2)*Mi^2);
denom1 = 1 + (gamma_i*Mi^2*(1-CD/2)^2);
phi = (gamma_i/gamma_e)*(num1/denom1)*T_ratio;

%find exit Mach using phi
denom2 = 1 - (2*gamma_e*phi) + sqrt(1 - (2*(gamma_e + 1)*phi));
Me = sqrt((2*phi)/denom2);

%find static pressure ratio using exit Mach
num3 = 1 + (gamma_i*Mi^2*(1-(CD/2)));
denom3 = 1 + (gamma_e*Me^2);
P_static_ratio = num3/denom3;

%find total pressure ratio using static pressure ratio and Mach
num4 = (1 + ((gamma_e-1)/2)*Me^2)^(gamma_e/(gamma_e-1));
denom4 = (1 + ((gamma_i-1)/2)*Mi^2)^(gamma_i/(gamma_i-1));
P_tot_ratio = P_static_ratio*(num4/denom4);
%display results
disp("Burner exit Mach number is " + string(Me));
disp("Burner total pressure ratio is " + string(P_tot_ratio));

