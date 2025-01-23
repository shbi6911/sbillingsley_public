%By:        Shane Billingsley
%Class:     ASEN 2802 Aerospace Sciences Lab 1
%Date:      Fall 2022

%% Import Helium Balloon Data
opts = detectImportOptions("Balloon Prelab Data.xlsx");
opts.DataRange = 'B2';
opts.VariableNamesRange = 'B1';
%opts
%opts.VariableNames;
data = readtable("Balloon Prelab Data.xlsx",opts);
uncert = data(4,1:12);
data = data(1:3,1:12);

%% Move Data to Vectors and Convert Units

%balloon mass without helium in kg
totalmass = (data.BalloonMass_g_+data.StringMass_g_+data.AdditionalMass_g_)/1000;
totalmass_uncert = (uncert.BalloonMass_g_+uncert.BalloonMass_g_+uncert.AdditionalMass_g_)/1000;

%exterior temperature in Kelvin
ext_temp = (((data.PILOTTemp__F_-32)*5)/9)+273.15;
ext_temp_uncert = (((uncert.PILOTTemp__F_-32)*5)/9)+273.15;

%exterior pressure in pascals
ext_pres = str2double(data.PILOTPressure_hPa_)*100;
ext_pres_uncert = str2double(uncert.PILOTPressure_hPa_)*100;

%tank dimensions in meters
height = data.WaterHeightDisplacement_mm_/1000;
height_uncert = uncert.WaterHeightDisplacement_mm_/1000;
length = (str2double(data.TankLength_in_)*2.54)/100;
length_uncert = ((1/16)*2.54)/100;
width = (str2double(data.TankWidth_in_)*2.54)/100;
width_uncert = ((1/16)*2.54)/100;
%calculate balloon volume in meters
volume = length.*width.*height;

%% Define Constants per Engineering Toolbox
sg_He = 0.138; %specific gravity of helium
R_air = 287.05; %in J/kg*K
R_He = 2077.1;  %in J/kg*K
g = 9.81; %gravitational acceleration in m/s^2

%% Calculate Mass of Helium

%mass of helium in each balloon in kg, calculated using ideal gas law
mass_helium = (sg_He*(ext_pres./(R_air*ext_temp))).*volume;
%calculate helium masses using air density and total mass instead
mass_helium_check = ((ext_pres./(R_air*ext_temp)).*volume)-totalmass;

%% Calculate Buoyant Force
rho_air = ext_pres./(R_air*ext_temp);
rho_He = ext_pres./(R_He*ext_temp);
F_buoy_calc = (volume*g).*(rho_air - rho_He);
F_buoy_meas = totalmass*g;
disp("Mean of Measured Values");
disp(mean(F_buoy_meas));
disp("Mean of Calculated Values");
disp(mean(F_buoy_calc));

%% Plot Forces Calculated vs. Measured
scatter([1 2 3],F_buoy_calc,20,'red','filled');
hold on;
grid on;
scatter([1 2 3],F_buoy_meas,20,'blue','filled');
scatter(4,mean(F_buoy_calc),20,'green','filled');
scatter(4,mean(F_buoy_meas),20,'magenta','filled');
xlim([0 5]);
ylim([0.08 0.12]);
xlabel("Balloon Trial");
ylabel("Buoyant Force (N)");
xticks([1 2 3 4]);
title("Calculated Buoyant Force vs. Measured Buoyant Force");
legend("Calculated Values","Measured Values","Calculated Mean","Measured Mean");

%% plot calculated masses of helium
% bar(mass_helium,'red');
% hold on;
% bar(mass_helium_check,'blue');
% xlabel("Balloons");
% ylabel("Mass of Helium (kg)");
% legend;


