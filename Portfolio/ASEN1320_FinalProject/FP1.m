%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

%This script calculates an orbital trajectory for two satellites over the
%period of one day, using the OrbitEOM function.  It then writes the output
%position data to a CSV file for each satellite, and plots the calculated 
%orbits%

%% initial values
sat1 = [1986.21;6388.28;-1237.15;-4.931;0.403;-5.829];
sat2 = [6480.84;1108.20;-2145.48;-0.289;7.071;2.746];

%% calculate orbits
orbitfun = @OrbitEOM;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspan= 0:60:86400;
[~, Sat1Full] = ode45(orbitfun,tspan,sat1,options);
[~, Sat2Full] = ode45(orbitfun,tspan,sat2,options);
%% plotting

[X,Y,Z] = sphere(100);
hold on; grid on;
plot3(-4460.49,2682.08,-3674.26,'.','MarkerSize',30,'Color',[0.9290 0.6940 0.1250]);
plot3(Sat1Full(:,1),Sat1Full(:,2),Sat1Full(:,3),'r','LineWidth',2);
plot3(Sat2Full(:,1),Sat2Full(:,2),Sat2Full(:,3),'k','LineWidth',2);
surf((X*6378),(Y*6378),(Z*6378),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.35);
plot3(0,0,0,'r.','MarkerSize',30);
axis('equal');
legend ('Canberra','ISS','Hubble','Earth Surface','Earth Center');
title('Plot of ISS and Hubble Telescope orbits, with Canberra Ground Station');
xlabel('x-distance from Earth Center (meters)');
ylabel('y-distance from Earth Center (meters)');
zlabel('z-distance from Earth Center (meters)');

%% create position matrices and output to CSV
Sat1Position = [Sat1Full(:,1), Sat1Full(:,2), Sat1Full(:,3)];
Sat2Position = [Sat2Full(:,1), Sat2Full(:,2), Sat2Full(:,3)];

writematrix(Sat1Position,'Sat1Position.csv');
writematrix(Sat2Position,'Sat2Position.csv');