%By:        Shane Billingsley
%Class:     ASEN 3728 Aircraft Dynamics
%Date:      Spring 2024

Q = [cosd(-60),sind(-60),0; -sind(-60),cosd(-60),0; 0,0,1];
V_B = [100;0;0];
W_E = [0;30;0];
V_E = Q*V_B - W_E;
time = 20000/V_E(2);



