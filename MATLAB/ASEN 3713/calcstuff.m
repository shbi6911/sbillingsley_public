%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

function [nusselt1, nusselt2] = calcstuff(T)
%this function calculates nusselt numbers using two different methodologies
%for a given plate submerged in water
%
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
    raynum = numerator/(nu^2);

    %find nusselt number #1
    q_dot = (k_plate/t)*(T_sh-T);
    nusselt1 = (q_dot*L)/(k_water*(((T_sh+T)/2)-T_infty));

    %find nusselt number #2
    coeff1 = (0.492/Pr)^(9/16);         %find nusselt number with steps
    denominator = (1+coeff1)^(8/27);
    numerator = 0.387*raynum^(1/6);
    nusselt2 = (0.825+(numerator/denominator))^2;

end