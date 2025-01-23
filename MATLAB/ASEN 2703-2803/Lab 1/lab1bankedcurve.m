%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

% basis vectors in car-centered frame
car_i = [1 0 0];
car_j = [0 1 0];
car_k = [0 0 1];

% define banked turn
theta_deg = 60;
theta = theta_deg*(pi/180);
rho = 20;
phi = linspace(0,pi,1000);
X = rho*sin(phi); Y = rho*cos(phi);
Z = zeros(1,1000);


% define other constants
m = 1;
g = 9.81;
h_0 = 125;
h = 107.6795;
v = sqrt(2*g*(h_0-h));

% define unit vectors for forces in car-centered frame
w = [0 sin(theta) cos(theta)];
ac = [0 cos(theta) -sin(theta)];

% define known vectors
weight = m*g*w;
centri = (((v^2)/rho)*m)*ac;
S = (weight(2)+centri(2))*-car_j;
N = (weight(3)+centri(3))*-car_k;
G = abs((S+N)/(m*g));

%output G-loading
format = "G-loading is %2.8f G in the %s direction \n";
fprintf(format,(-G(2)),"leftward");
fprintf(format,(-G(3)),"upward");
fprintf("Total Length of Track is %3.2f m \n",(pi*rho));

%plot in 3D
plot3(X,Y,Z);
grid on;
