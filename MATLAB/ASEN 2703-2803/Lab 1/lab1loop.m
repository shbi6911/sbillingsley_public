%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

% basis vectors in car-centered frame
car_i = [1 0 0];
car_j = [0 1 0];
car_k = [0 0 1];

% define loop
h_loop = 25;
rho = 40;
phi = linspace(0,2*pi,1000);
Z = (rho*-cos(phi)+rho)+h_loop; Y = rho*sin(phi);
X = ones(1,1000)+25;

% define other constants
m = 1;
g = 9.81;
h_0 = 125;
v = sqrt(2*g*(h_0-Z));

%define G-forces
G = zeros(3,1000);
G(3,:)=(v.^2)./(rho*g)+cos(phi);

%output G-loading
[Mk,Ik] = max(G(3,:));
format = "Maximum G-loading is %2.8f G in the %s direction at theta= %2.2f \n"; 
fprintf(format,Mk,"upward",phi(Ik));
fprintf("Total Length of Track is %3.2f m \n",(pi*rho));

%check surfaces
X_check = -50:50;
Y_check = -50:50;
Z_check1 = ones(101,101)+24;
Z_check2 = ones(101,101)+124;

%plot in 3D
hold on;
plot3(X,Y,Z);
surf(X_check,Y_check,Z_check1);
surf(X_check,Y_check,Z_check2);
grid on;
