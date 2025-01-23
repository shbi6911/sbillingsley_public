%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

clear;
%define constants
A_s = 8*2.8;    %surface area
Qin = 200;      %radiant heat going in
T_infty = 30;   %ambient air temperature
L = 8;          %length
V = 70*1000*(1/3600); %velocity in m/s
T_s = 40;           %define initial guess
error = 1;          %initialize error value
error_tol = 1E-6;   %define tolerance to be achieved

while abs(error) > error_tol
    [~,~,k,~,~,v,Pr] = get_A22_properties(T_s); %pull properties at current temp
    Re_L = (V*L)/v;                             %calculate current Reynolds #
    %find heat flow out under these conditions
    Qout=((k/L)*(0.037*Re_L^0.8 -871)*(Pr^(1/3)))*(T_s-T_infty);
    error = Qin-Qout;       %compare to radiative heat input
    step = abs(error)*0.01;     %define a step size to change temp guess
    if error >= 0           %increase or decrease temp guess
        T_s = T_s + step;
    else
        T_s = T_s - step;
    end
end
%print results to the screen
disp("The temperature of the train car is " + string(T_s) + ".")
    
