%By:        Shane Billingsley
%Class:     ASEN 2802 Aerospace Sciences Lab 1
%Date:      Fall 2022

%% Define constants
height = 0:10:30000;
rho = 0.137;
m = 0.175e-2;
rHe = 2077.1;
%% Allocate temperature and Pressure data
[T,~,P] = atmoscoesa(height);
%% Calculate Volume Array
vol = m * rHe * T ./ P;
%% Plot Volume Array vs Height
plot(height, vol, "LineWidth", 4);
xlabel("Altitude (m)");
ylabel("Volume (m^3)");
title("Volume vs. Altitude");