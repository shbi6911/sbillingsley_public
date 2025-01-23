%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

% define banked turn
theta_deg = 60;
theta = theta_deg*(pi/180);
rho = 40;
phi = linspace(0,pi,1000);
X = rho*sin(phi); Y = rho*cos(phi);
Z = zeros(1,1000);

% define other constants
g = 9.81;
h_0 = 125;
h = 25;
v = sqrt(2*g*(h_0-h));

centri = v^2/(g*rho);
Gs = [0 (sin(theta)-(centri*cos(theta))) (cos(theta)+(centri*sin(theta)))];
disp(Gs);