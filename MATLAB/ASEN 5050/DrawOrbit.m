%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

% a script for plotting an orbit using f and g functions
clear; clc; const = getConst();

%define central body parameters
mu = const.mu.Earth; %grav parameter in km^3/s^2
radius = const.radius.Earth; %radius in km
n = 1000;       %number of points on the orbit (resolution)

%define initial state (Keplerian)
a = 10000;   e = 0.1;    i = deg2rad(20);
RAAN = deg2rad(45);   omega = deg2rad(90);  theta_star = deg2rad(90);
[Rvec,Vvec] = KeplerToCartesian(a,e,i,RAAN,omega,theta_star,mu);

% %define initial state (Cartesian)
% Rvec = [5.352950e6; 7.053778e5; -4.059700e5];  %orbital radius in km, inertial
% Vvec = [-4.164248; 1.963690; 0.3191257]; %velocity in km/s

%define 6xn array of n state values along the orbit
orbit_state = FandG_Orbit(Rvec,Vvec,mu,n);

[X,Y,Z] = sphere(100);  %define a sphere for the central body

hold on; grid on;
%plot orbit and central body
plot3(orbit_state(1,:),orbit_state(2,:),orbit_state(3,:),'r','LineWidth',2);    
surf((X*radius),(Y*radius),(Z*radius),...
    'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.15,...
    'EdgeAlpha',0.15);
plot3(0,0,0,'r.','MarkerSize',30);
axis('equal');
ax = gca;
view(20,10);

%plot an xy plane
[x,y] = meshgrid(ax.XLim,ax.YLim); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z,'FaceAlpha',0.35) % Plot the surface

%plot axes
plot3([0;ax.XLim(2)],[0;0],[0;0], 'b','LineWidth',2);
plot3([0;0],[0;ax.YLim(2)],[0;0], 'm','LineWidth',2);
plot3([0;0],[0;0],[0;ax.ZLim(2)], 'g','LineWidth',2);