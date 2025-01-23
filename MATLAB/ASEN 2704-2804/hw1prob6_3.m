%By:        Shane Billingsley
%Class:     ASEN 2804 Vehicle Design Lab
%Date:      Spring 2023

% velocity
V = linspace(10,310,100);
% fixed constants
S = 47; AR = 6.5; eo = 0.87; W=103047; C_d0=0.032; T_A=2*40298;
rho_sl=1.225; rho_5km=0.73643; k=1/(pi*eo*AR); T_A5km=T_A*(rho_5km/rho_sl);
% densities
q_sl = 0.5*rho_sl*V.^2; q_5km=0.5*rho_5km*V.^2;
%power @ sea level
P_Asl = T_A*V; P_A5km = T_A5km*V;
T_Rsl = (q_sl*S*C_d0) + ((k*W^2)./(q_sl*S));
%power @ 5 km
T_R5km = (q_5km*S*C_d0) + ((k*W^2)./(q_5km*S));
P_Rsl=T_Rsl.*V; P_R5km=T_R5km.*V;

%calculate maximum excess power
[EP_sl, EP_sl_I] = max(P_Asl-P_Rsl);
[EP_5km, EP_5km_I] = max(P_A5km-P_R5km);
%output max excess power figures
fprintf("Max Excess Power at Sea Level is %8.0f watts \n",EP_sl);
fprintf("Max Excess Power at 5 km is %8.0f watts \n",EP_5km);

% plot sea level
figure (1); hold on;
plot(V,P_Rsl);
plot(V,P_Asl);
xline(V(EP_sl_I));
title("Power Required vs. Power Available at Sea Level");
xlabel("Velocity (m/s)");
ylabel("Power (watts)");
legend("Power Required","Power Available","Max Excess Power");
hold off;

% plot 5 km
figure (2); hold on;
plot(V,P_R5km);
plot(V,P_A5km);
xline(V(EP_5km_I));
title("Power Required vs. Power Available at 5 km Altitude");
xlabel("Velocity (m/s)");
ylabel("Power (watts)");
legend("Power Required","Power Available","Max Excess Power");
hold off;


