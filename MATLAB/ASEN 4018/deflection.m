%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

clear; clc;

L = 1;      h = 0.1;    %length and deflection height

theta = atan2(h,L);     %deflection angle

x = L/cos(theta) - L;   %reduction in length due to deflection
disp(x);