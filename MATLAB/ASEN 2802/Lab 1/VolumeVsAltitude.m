%By:        Shane Billingsley
%Class:     ASEN 2802 Aerospace Sciences Lab 1
%Date:      Fall 2022

%% Define constants
height = 0:10:15000;
rho = 0.137;
v0 = 0.105 * (19.5 * 0.0254) * (9.75 * 0.0254);
molMass = 4002.602;
m = rho * v0;
n = m / molMass;
rHe = 2076.9;
%% Allocate temperature and Pressure data
[T,~,P] = atmoscoesa(height);
%% Calculate Volume Array
vol = n * rHe * T ./ P;
%% Plot Volume Array vs Height
plot(height, vol, "LineWidth", 4);
xlabel("Altitude (m)");
ylabel("Volume (m^3)");
title("Volume vs. Altitude");