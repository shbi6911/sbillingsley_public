%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

%Problem 20-16

%initial guess for temperature of cold side of plate, in Kelvin
guess = 53.5;  %given initial temperature guess
err = 1;    %initialize an error value

while abs(err) > 0.000000001    
    temp_new = calcstuff(guess);    %find new temp
    err = (temp_new-guess)/guess;   %calculate error between new and old
    guess = temp_new;   %set new guess
end
disp("The surface temperature is " + string(temp_new));

function temp = calcstuff(T)
    %set constants
    L = 0.2;    %plate length in meters
    t = 0.025;   %plate thickness in meters
    T_infty = 7;    %temperature of cold water in Kelvin
    T_sh = 100;      %temperature of hot side of plate in Kelvin
    k_plate = 15;       %thermal conductivity of plate in W/m*K
    g = 9.81;   %acceleration of gravity
    filmtemp = (T +T_infty)/2;    %film temperature of the water
    [k_water, nu, Pr, beta] = get_A15_properties(filmtemp);   %get water properties

    %find Rayleigh #
    numerator = g*beta*(T - T_infty)*(L^3);
    raynum = (numerator/(nu^2))*Pr;

    %find nusselt number by stepwise approach
    coeff1 = (0.492/Pr)^(9/16);
    denominator = (1+coeff1)^(8/27);
    numerator = 0.387*raynum^(1/6);
    nusselt = (0.825+(numerator/denominator))^2;

    h = (nusselt*k_water)/L;
    R = (t/k_plate)+(1/h);
    Q_dot = (T_sh-T_infty)/R;
    temp = (Q_dot/h)+7;
end

