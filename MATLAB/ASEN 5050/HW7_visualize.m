%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

clear; clc; const = getConst();
%givens
r_p = 600000;           %km
r_a = 1800000;
a_titan = 1221830;
r_titan = 2575;
radius = const.radius.Saturn;

m_titan = 1.3455e23;    %kg

%characterize spacecraft orbit around Saturn
a_1 = (r_p+r_a)/2;
e_1 = 1 - (r_p/a_1);
p_1 = a_1*(1 - e_1^2);
h_1 = sqrt(p_1*const.mu.Saturn);
v_p = sqrt((2*const.mu.Saturn/r_p) - (const.mu.Saturn/a_1));
R_p_sc1 = [r_p;0;0];
V_p_sc1 = [0;v_p;0];
%find true anomaly at Titan orbit intersection
theta_star_1 = -acos((p_1 - a_titan)/(a_titan*e_1));
%construct radial vector using rotation matrix
angle1 = -(theta_star_1 + (pi/2));
Rot = [cos(-angle1),-sin(-angle1),0;sin(-angle1),cos(-angle1),0;0,0,1];
R1_xy = Rot*[a_titan;0;0];

%find velocity and flight path angle
v_1 = sqrt(const.mu.Saturn*((2/a_titan)-(1/a_1)));
gamma = -acos(h_1/(a_titan*v_1));
%find components of velocity and express as a vector
v_r1 = v_1*sin(gamma);   v_theta1 = v_1*cos(gamma);
V1_roh = [v_r1;v_theta1;0];
%characterize velocity in xy axes
V1_xy = Rot*V1_roh;


%velocity of Titan relative to Saturn
v_titan = sqrt(const.mu.Saturn/a_titan);
V_titan = [0;v_titan;0];
R_titan = [a_titan;0;0];

V_infty_in = V1_roh - V_titan;  %velocity of s/c relative to Titan

%characterize hyperbolic trajectory relative to Titan
r_p_h = 3000;       %km
mu_titan =  m_titan*const.G;
a_h = -mu_titan/norm(V_infty_in)^2;
e_h = 1 - (r_p_h/a_h);
delta = 2*asin(1/e_h);

%find outgoing velocity relative to Saturn
beta_1 = pi - acos(dot(V_infty_in,V_titan)/(norm(V_infty_in)*norm(v_titan)));
beta_2 = beta_1 - delta;
V_infty_out = [-norm(V_infty_in)*sin(beta_2);-norm(V_infty_in)*cos(beta_2);0];
V_out = V_titan + V_infty_out;

%characterize new orbit
h_2 = norm(cross([a_titan;0;0],V_out));
a_2 = -const.mu.Saturn/(dot(V_out,V_out) - (2*const.mu.Saturn)/a_titan);
e_2 = sqrt(1 - (h_2^2/(const.mu.Saturn*a_2)));
p_2 = a_2*(1-e_2^2);
theta_star_2 = -acos((p_2 - a_titan)/(a_titan*e_2));
%find equivalent impulsive maneuver
delta_V = V_out - V1_roh;

%find new orbit vectors
r_p2 = p_2/(1 - e_2);
v_p2 = h_2/r_p2;
V_p_sc2 = [0;v_p2;0];
R_p_sc2 = [r_p2;0;0];
V2_xy = Rot*V_out;

%define 6xn array of n state values along the orbit
orbit_state_titan = FandG_Orbit(R_titan,V_titan,const.mu.Saturn,1000);

%define 6xn array of n state values along the orbit
orbit_state_1 = FandG_Orbit(R1_xy,V1_xy,const.mu.Saturn,1000);

%define 6xn array of n state values along the orbit
orbit_state_2 = FandG_Orbit(R1_xy,V2_xy,const.mu.Saturn,1000);

%% plotting
[X,Y,Z] = sphere(100);  %define a sphere for the central body

figure();
hold on; grid on;
%plot orbit and central body
plot3(orbit_state_titan(1,:),orbit_state_titan(2,:),orbit_state_titan(3,:),'r','LineWidth',2);    
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