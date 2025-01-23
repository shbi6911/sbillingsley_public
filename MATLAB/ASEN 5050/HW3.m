%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1

clear; clc;
const.mu.Saturn = 3.794e7; %grav parameter in km^3/s^2
const.radius.Saturn = 60268; %radius of Mars in km
Rvec_1 = [-720000; 670000; 310000];  %orbital radius in km
Vvec_1 = [2.160; -3.360; 0.620]; %velocity in km/s
R_1 = norm(Rvec_1); V_1 = norm(Vvec_1);     %magnitudes

%define h vector
Hvec_1 = cross(Rvec_1, Vvec_1);  H_1 = norm(Hvec_1);
%define orbital energy
energy_1 = V_1^2/2 - const.mu.Saturn/R_1;
%define semi-major axis
a_1 = -const.mu.Saturn/(2*energy_1);
%find eccentricity vector, unit vector, and scalar eccentricity
evec_1 = cross(Vvec_1,Hvec_1)./const.mu.Saturn - Rvec_1./R_1;
e_1 = norm(evec_1);     e_hat_1 = evec_1./e_1;
%find line of nodes vector and unit vector
nvec_1 = cross([0;0;1], Hvec_1);  n_hat_1 = nvec_1./norm(nvec_1);

%inclination
i_1 = acos(Hvec_1(3)/H_1);  i_1_deg = rad2deg(i_1);

%RAAN
OMEGA_1 = acos(n_hat_1(1));
if n_hat_1(2) < 0
    OMEGA_1 = -OMEGA_1;
end
OMEGA_1_deg = rad2deg(OMEGA_1);

%argument of periapsis
omega_1 = acos(dot(n_hat_1,e_hat_1));
if e_hat_1(3) < 0
    omega_1 = -omega_1;
end
omega_1_deg = rad2deg(omega_1);

%impact point approximated as where orbital radius = planet radius
R_2 = const.radius.Saturn;
%find true anomaly at this point (negative b/c moving toward periapsis)
theta_star_2 = -acos((a_1*(1-e_1^2) - R_2)/(R_2*e_1));

%construct r and v in r/theta/h axes
rvec_2 = [R_2;0;0];
v_r_2 = (const.mu.Saturn/H_1)*e_1*sin(theta_star_2);
v_theta_2 = (const.mu.Saturn/H_1)*(1+(e_1*cos(theta_star_2)));
vvec_2 = [v_r_2;v_theta_2;0];

%construct DCM and rotate r and v into XYZ axes
DCM = AstroRot313(OMEGA_1,i_1,(theta_star_2+omega_1));
Rvec_2 = DCM*rvec_2;
Vvec_2 = DCM*vvec_2;

%% Problem 2
clear; clc;
%givens
a_2 = 6463.8;   e_2 = 0.45454;  i_2 = deg2rad(74.924);
RAAN_2 = deg2rad(1.2410);   omega_2 = deg2rad(353.31);  theta_star_2 = deg2rad(199.38);
const.mu.Mars = 42828.314258067;
const.radius.Mars = 3397.2;
%radius at periapsis
r_p_2 = a_2*(1-e_2);
alt_p_2 = r_p_2 - const.radius.Mars;
%orbital period
P_2 = 2*pi*sqrt(a_2^3/const.mu.Mars);
%convert given Keplerian elements to Cartesian vectors
[Rvec_2,Vvec_2] = KeplerToCartesian(a_2,e_2,i_2,RAAN_2,omega_2,theta_star_2,const.mu.Mars);

%load in and plot data from GMAT, point mass run
ECCdata = readmatrix('MarsPointMassECCPlot_1');
HMAGdata = readmatrix('MarsPointMassHMAGPlot_1');
figure();
plot(ECCdata(:,1),ECCdata(:,2),'r');
ECCavg = mean(ECCdata(:,2));
xlabel("Elapsed Days");
ylabel("Eccentricity");
title("Eccentricity vs Time, Point Mass");
ylim([ECCavg-ECCavg*0.000000001,ECCavg+ECCavg*0.000000001]);
figure();
plot(HMAGdata(:,1),HMAGdata(:,2),'b');
HMAGavg = mean(HMAGdata(:,2));
xlabel("Elapsed Days");
ylabel("Specific Angular Momentum (km^2/s)");
title("Specific Angular Momentum vs Time, Point Mass");
ylim([HMAGavg-HMAGavg*0.000000001,HMAGavg+HMAGavg*0.000000001]);

%load in and plot data from GMAT, high fidelity run
ECCdata_2 = readmatrix('MarsPointMassECCPlot_2');
HMAGdata_2 = readmatrix('MarsPointMassHMAGPlot_2');
figure();
plot(ECCdata_2(:,1),ECCdata_2(:,2),'r');
ECCavg_2 = mean(ECCdata_2(:,2));
xlabel("Elapsed Days");
ylabel("Eccentricity");
title("Eccentricity vs Time, High Fidelity");
%ylim([ECCavg_2-ECCavg_2*0.000000001,ECCavg_2+ECCavg_2*0.000000001]);
figure();
plot(HMAGdata_2(:,1),HMAGdata_2(:,2),'b');
HMAGavg_2 = mean(HMAGdata_2(:,2));
xlabel("Elapsed Days");
ylabel("Specific Angular Momentum (km^2/s)");
title("Specific Angular Momentum vs Time, High Fidelity");
%ylim([HMAGavg_2-HMAGavg_2*0.000000001,HMAGavg_2+HMAGavg_2*0.000000001]);