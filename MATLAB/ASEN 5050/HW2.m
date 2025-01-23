%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc;
const.mu.Mars = 4.305e4; %grav parameter in km^3/s^2
const.radius.Mars = 3397.2; %radius of Mars in km
Rvec_1 = [3.62067e3; -3.19925e2; -4.20645e2];  %orbital radius in km
Vvec_1 = [-4.28843e-1; -3.00176e-2; -3.39801]; %velocity in km/s
R_1 = norm(Rvec_1); V_1 = norm(Vvec_1);     %magnitudes

%define h vector
Hvec_1 = cross(Rvec_1, Vvec_1);  H_1 = norm(Hvec_1);
%define orbital energy
energy_1 = V_1^2/2 - const.mu.Mars/R_1;
%define semi-major axis
a_1 = -const.mu.Mars/(2*energy_1);
%find eccentricity vector, unit vector, and scalar eccentricity
evec_1 = cross(Vvec_1,Hvec_1)./const.mu.Mars - Rvec_1./R_1;
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

%true anomaly
theta_star_1 = acos(dot((Rvec_1./R_1),e_hat_1));
if dot(Vvec_1, Rvec_1) < 0
    theta_star_1 = -theta_star_1;
end
theta_star_1_deg = rad2deg(theta_star_1);

%construct DCM and rotate r and v into r/theta/h axes
DCM = AstroRot313(OMEGA_1,i_1,(theta_star_1+omega_1));
DCMT = DCM';
rvec_1 = DCMT*Rvec_1;
vvec_1 = DCMT*Vvec_1;

%construct r and v in r/theta/h axes by earlier method
rvec_1_check = [R_1;0;0];
v_r_1 = dot((Rvec_1/R_1),Vvec_1);
gamma_1 = asin(v_r_1/V_1);
v_theta_1 = V_1*cos(gamma_1);
vvec_1_check = [v_r_1;v_theta_1;0];

%find position and velocity vectors at periapsis
r_p_1 = a_1*(1-e_1);
v_p_1 = sqrt(2*energy_1 + 2*(const.mu.Mars/r_p_1));
rvec_p_1 = [r_p_1;0;0];
vvec_p_1 = [0;v_p_1;0];
DCM_p = AstroRot313(OMEGA_1,i_1,omega_1);
Rvec_p_1 = DCM_p*rvec_p_1;
Vvec_p_1 = DCM_p*vvec_p_1;

%% Problem 2
n_hat_2 = [0.6428;-0.7660;0];       %given vectors
h_hat_2 = [-0.3237;-0.2717;0.9063];
evec_2 = [0.0475;0.3755;0.1295];

e_2 = norm(evec_2);     %find eccentricity
e_hat_2 = evec_2./e_2;  %eccentricity unit vector

%inclination
i_2 = acos(h_hat_2(3));  i_2_deg = rad2deg(i_2);

%RAAN
OMEGA_2 = acos(n_hat_2(1));
if n_hat_2(2) < 0
    OMEGA_2 = -OMEGA_2;
end
OMEGA_2_deg = rad2deg(OMEGA_2);

%argument of periapsis
omega_2 = acos(dot(n_hat_2,e_hat_2));
if e_hat_2(3) < 0
    omega_2 = -omega_2;
end
omega_2_deg = rad2deg(omega_2);

%given radius at ascending node
r_2 = 4070.6;   %km
%at ascending node, arg of perigee is the negative of true anomaly
theta_star_2 = -omega_2;
%use conic eqn
a_2 = (r_2+(r_2*e_2*cos(theta_star_2)))/(1-e_2^2);