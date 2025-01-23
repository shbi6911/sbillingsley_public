%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023


%Problem 20-30
K_conv = 273.15;    %celsius to kelvin
T_film = 32.5;      %given film temperature
T_infty = 20;   %temperature of the air(C)
T_sg = T_film*2 - T_infty;   %assumed surface temp based on given film temp
w = 0.15;       %length of the plate(m)
l = 0.2;        %width of the plate(m)
p = 2*l + 2*w;  %perimeter of the plate(m)
A_s = l*w;      %area of the plate (m^2)
L_c = A_s/p;    %characteristic length for horizontal plate(m)
[rho, cp, k, alpha, mu, nu, Pr] = get_A22_properties(T_film);
beta = 1/(T_film+K_conv);  %beta is 1/T for ideal gas
Q_dot = 8;      %given heat dissipation(W)
g = 9.81;   %acceleration of gravity(m/s)
eps = 0.8;  %emissivity of the plate
sigma = 5.67*10^-8; %Stefan-Boltzmannn constant

%Vertical Plate
%find Rayleigh #
numerator = g*beta*(T_sg - T_infty)*(l^3);
raynum = (numerator/(nu^2))*Pr;

%find nusselt number by stepwise approach
coeff1 = (0.492/Pr)^(9/16);
denominator = (1+coeff1)^(8/27);
numerator = 0.387*raynum^(1/6);
nusselt = (0.825+(numerator/denominator))^2;

%find h based on nusselt number
h = (nusselt*k)/l;

%solve using vpasolve
syms T_s_sym
T_s = double(vpasolve(Q_dot==h*A_s*(T_s_sym-T_infty)+...
    sigma*A_s*eps*((T_s_sym+K_conv)^4-(T_infty+K_conv)^4),T_s_sym,[-273.15,Inf]));
disp("The surface temperature of Case 1 (vertical plate) is "+string(T_s));

%Horizontal Plate(hot side UP)

%find Rayleigh #
numerator = g*beta*(T_sg - T_infty)*(L_c^3);
raynum = (numerator/(nu^2))*Pr;

%find nusselt number based on two cases
if 10^4 < raynum && raynum < 10^7
    nusselt = 0.54*raynum^(1/4);
else
    nusselt = 0.15*raynum^(1/3);
end

%find h based on nusselt number
h = (nusselt*k)/L_c;

%solve using vpasolve
syms T_s_sym
T_s = double(vpasolve(Q_dot==h*A_s*(T_s_sym-T_infty)+...
    sigma*A_s*eps*((T_s_sym+K_conv)^4-(T_infty+K_conv)^4),T_s_sym,[-273.15,Inf]));
disp("The surface temperature of Case 2 (horizontal plate, hot side up) is "+string(T_s));

%Horizontal Plate (hot side DOWN)

%find nusselt number
nusselt = 0.27*raynum^(1/4);

%find h based on nusselt number
h = (nusselt*k)/L_c;

%solve using vpasolve
syms T_s_sym
T_s = double(vpasolve(Q_dot==h*A_s*(T_s_sym-T_infty)+...
    sigma*A_s*eps*((T_s_sym+K_conv)^4-(T_infty+K_conv)^4),T_s_sym,[-273.15,Inf]));
disp("The surface temperature of Case 2 (horizontal plate, hot side down) is "+string(T_s));
