%By:        Shane Billingsley
%Class:     ASEN 2802 Aerospace Sciences Lab 1
%Date:      Fall 2022

%% Define Variables
altitude = 0:10:15000; %m
v = 4.6 * 0.028; %m^3
rAir = 287.0; %J/kg*K
m = 0.03; %kg
%% Use Atmoscoesa to Generate Pressure and Density of Surrounding Air
[~,~,P,rho] = atmoscoesa(altitude);
%% Calculate Temperature Vector
T = (v * P) ./ (rAir * (v * rho - m));
% Plot Temp vs Altitude
a = plot(altitude, T,"LineWidth",4);
hold on;
b = xline(10100, "LineStyle", "--", "Color", "r");
c = yline(522, "LineStyle", "--", "Color", "g");
grid on;
xlim([1 12000]);
ylim([300 550]);
xlabel("Altitude (m)");
ylabel("Temperature (K)");
title("Altitude (m) vs Temperature (K)");
legend([a b c], ["Altitude Vs. Pressure" "Max Altitude" "Max Temperature"], "Location", "southwest");

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
grid on;
xlabel("Altitude (m)");
ylabel("Volume (m^3)");
title("Altitude (m) vs. Volume (m^3)");
