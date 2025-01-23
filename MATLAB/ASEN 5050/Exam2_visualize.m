%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc; const = getConst();

%characterize initial orbit
r_p1 = 1850;    r_a1 = 5400;    r1 = 3500;      %km
e1 = (r_a1 - r_p1)/(r_a1 + r_p1);   %eccentricity
a1 = (r_a1+r_p1)/2;                 %semimajor axis
p1 = a1*(1 - e1^2);                 %semilatus rectum
h1 = sqrt(const.mu.Moon*p1);        %specific angular momentum
v1 = sqrt(((2*const.mu.Moon)/r1) - (const.mu.Moon/a1)); %velocity
gamma1 = -acos(h1/(r1*v1));         %flight path angle
V1 = [v1*sin(gamma1);v1*cos(gamma1);0];     %velocity vector
theta_star1 = -acos((p1 - r1)/(r1*e1));  %true anomaly

%construct radial vector using rotation matrix
angle1 = (-theta_star1 - (pi/2));
Rot = [cos(angle1),-sin(angle1),0;sin(angle1),cos(angle1),0;0,0,1];
R1_xy = Rot*[r1;0;0];

%characterize velocity in xy axes
V1_xy = Rot*V1;

%characterize final orbit
r_p2 = 3000;    r_a2 = 5400;    r2 = 3500;      %km
e2 = (r_a2 - r_p2)/(r_a2 + r_p2);   %eccentricity
a2 = (r_a2+r_p2)/2;                 %semimajor axis
p2 = a2*(1 - e2^2);                 %semilatus rectum
h2 = sqrt(const.mu.Moon*p2);        %specific angular momentum
v2 = sqrt(((2*const.mu.Moon)/r2) - (const.mu.Moon/a2)); %velocity
gamma2 = -acos(h2/(r2*v2));         %flight path angle
V2 = [v2*sin(gamma2);v2*cos(gamma2);0];     %velocity vector
theta_star2 = -acos((p2 - r2)/(r2*e2));  %true anomaly

%construct radial vector using rotation matrix
angle2 = (-theta_star2 - (pi/2));
Rot = [cos(angle2),-sin(angle2),0;sin(angle2),cos(angle2),0;0,0,1];
R2_xy = Rot*[r2;0;0];

%characterize velocity in xy axes
V2_xy = Rot*V2;


%% Problem 2
clear; clc; const = getConst();

%givens
Vvec_in = [2.2789;5.8841;0];    %km/s
a1 = 1.8e6;
gamma = atan2(Vvec_in(1),Vvec_in(2));
denom = dot(Vvec_in,Vvec_in)/2 + const.mu.Jupiter/(2*a1);
r1 = const.mu.Jupiter/denom;
h1 = norm(cross([r1;0;0],Vvec_in));
e1 = sqrt(1 - (h1^2/(const.mu.Jupiter*a1)));
p1 = a1*(1 - e1^2);
theta_star1 = acos((p1 - r1)/(r1*e1));

%construct radial vector using rotation matrix
angle1 = -(theta_star1 + (pi/2));
Rot = [cos(angle1),-sin(angle1),0;sin(angle1),cos(angle1),0;0,0,1];
R1_xy = Rot*[r1;0;0];   R2_xy = R1_xy;

%characterize velocity in xy axes
V1_xy = Rot*Vvec_in;

v_out = 7.6759;     gamma_out = deg2rad(20.91);
V2 = [v_out*sin(gamma_out);v_out*cos(gamma_out);0];
V2_xy = Rot*V2;

%velocity of Titan relative to Saturn
v_moon = sqrt(const.mu.Jupiter/r1);
V_moon = [0;v_moon;0];
R_moon = [r1;0;0];

%% plotting

%define 6xn array of n state values along the orbit
orbit_state_moon = FandG_Orbit(R_moon,V_moon,const.mu.Jupiter,1000);

%define 6xn array of n state values along the orbit
orbit_state_1 = FandG_Orbit(R1_xy,V1_xy,const.mu.Jupiter,1000);

%define 6xn array of n state values along the orbit
orbit_state_2 = FandG_Orbit(R2_xy,V2_xy,const.mu.Jupiter,1000);


[X,Y,Z] = sphere(100);  %define a sphere for the central body
radius = const.radius.Moon;
figure();
hold on; grid on;
%plot orbit and central body
plot3(orbit_state_moon(1,:),orbit_state_moon(2,:),orbit_state_moon(3,:),'r','LineWidth',2);    
surf((X*radius),(Y*radius),(Z*radius),...
    'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.15,...
    'EdgeAlpha',0.15);
plot3(0,0,0,'r.','MarkerSize',30);
axis('equal');
ax = gca;
view(20,10);

plot3(orbit_state_1(1,:),orbit_state_1(2,:),orbit_state_1(3,:),'b','LineWidth',2);

plot3(orbit_state_2(1,:),orbit_state_2(2,:),orbit_state_2(3,:),'g','LineWidth',2);

%plot an xy plane
[x,y] = meshgrid(ax.XLim,ax.YLim); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z,'FaceAlpha',0.35) % Plot the surface

%plot axes
plot3([0;ax.XLim(2)],[0;0],[0;0], 'b','LineWidth',2);
plot3([0;0],[0;ax.YLim(2)],[0;0], 'm','LineWidth',2);
plot3([0;0],[0;0],[0;ax.ZLim(2)], 'g','LineWidth',2);