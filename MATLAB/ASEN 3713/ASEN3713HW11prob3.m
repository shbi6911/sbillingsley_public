%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

%Problem 20-60
tic;
K_conv = 273.15;    %celsius to kelvin
Q_dot = 54; %heat flux in watts
T_infty = 25;   %room temperature in celsius
d = 0.08;   %diameter in meters
A_s = 4*pi*(d/2)^2; %surface area of a sphere
eps = 0.9;  %emissivity of glass
sigma = 5.67*10^-8; %Stefan-Boltzmannn constant

guess = 1000;   %initial temperature guess
err = 1;    %initialize an error value
iter = 0;
while abs(err) > 0.000000001    
    h = calcstuff(guess);    %find h-value based on guess
    %find new temperature using vpasolve
    syms T_s_sym
    temp_new = double(vpasolve(Q_dot==h*A_s*(T_s_sym-T_infty)+...
    sigma*A_s*eps*((T_s_sym+K_conv)^4-(T_infty+K_conv)^4),...
    T_s_sym,[-273.15,Inf]));
    err = (temp_new-guess)/guess;   %calculate error between new and old
    guess = temp_new;   %set new guess
    iter = iter+1;
end
disp("The surface temperature is " + string(temp_new));
toc
disp(iter);
function h = calcstuff(T)
    %set constants
    K_conv = 273.15;    %celsius to kelvin
    d = 0.08;   %diameter in meters
    T_infty = 25;    %room temperature in Celsius
    g = 9.81;   %acceleration of gravity
    filmtemp = (T +T_infty)/2;    %film temperature of the water
    beta = 1/(filmtemp+K_conv);  %calculate beta for ideal gas
    %get air properties
    [~, ~, k, ~, ~, nu, Pr] = get_A22_properties(filmtemp);

    %find Rayleigh #
    numerator = g*beta*(T - T_infty)*(d^3);
    raynum = (numerator/(nu^2))*Pr;

    %find nusselt number by stepwise approach
    coeff1 = (0.469/Pr)^(9/16);
    denominator = (1+coeff1)^(4/9);
    numerator = 0.589*raynum^(1/4);
    nusselt = 2+(numerator/denominator);

    h = (nusselt*k)/d;
end
