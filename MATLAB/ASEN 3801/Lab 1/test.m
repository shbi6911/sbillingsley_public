[dens,~,~,~,~]=stdatmo(1655,0);  
g=9.81;             %acceleration of gravity in m/s
m = 0.05;           %mass of the object in kg
coeff_drag = 0.6;           %coefficient of drag of the object
cross_sec_a = pi*0.01^2;      %cross-sectional area of the object in m^2
tol = 1E-8;         %tolerance value for ode45
wind_vel = [0;0;0]; %initial wind velocity

x_ic = [0;0;0;0;20;-20];    %given initial state values (pos xyz, vel xyz)
                            % in meters and m/s

tspan = [0 20];             %time span for integration

%set options for ode45, including tolerances and stop event function
opts = odeset('Events',@HitGround,'RelTol',tol,'AbsTol',tol);
%run initial trajectory with no wind
[t_out_i,x_out_i,~,~,~] = ode45(@(t,state_vector) objectEOM(t,state_vector, dens, coeff_drag, cross_sec_a, m, g, wind_vel)...
    ,tspan,x_ic,opts);
%plot initial trajectory in 3D space
figure(); view(3);
hold on; grid on;
plot3(x_out_i(:,1),x_out_i(:,2),x_out_i(:,3),'LineWidth',2);
%plot3(x_out_basic(:,1),x_out_basic(:,2),x_out_basic(:,3),'LineWidth',2);
set(gca,'ZDir','reverse');  %reverse axes so that -z is on top
title("Initial Trajectory with Zero Wind Speed");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
%legend("Gravity and Drag","Gravity Only",'location','best');
hold off;

function [x_dot] = objectEOM(t,state_vector, dens, coeff_drag, cross_sec_a, m, g, wind_vel)
% Inputs:   
%           t = time
%           state_vector = state vector [pos ; vel]
%               pos = need inertial position [x,y,z]'
%               vel = and inertial velocity [vx,vy,vz]'
%           dens = density
%           coeff_drag = coefficient of drag
%           cross_sec_a = cross sectional area
%           m = mass
%           g = gravity
%           wind_vel = wind velocity vector [x,y,z]'           
% Outputs:  
%           x_dot = state vector [ V_e ; A_e ] 
%           V_e = [vx,vy,vz]' inertial velocity vector
%           A_e = [ax,ay,az]' inertial acceleration vector        
% Methodology: im working on it, okay. 
%              but above is what we need to put in, and get out.

%%
% velocity vector is the 3 velocity parts from the state vector.
vel = state_vector(4:6,1);
wind_vel = wind_vel';
% deal with wind vector
V_e = vel + wind_vel;  
% Airspeed
vel_air = norm(V_e);

%if state_vector(3,1)

% drag vector *negative, opposite direction of motion*
drag = -0.5* dens * vel_air * cross_sec_a * coeff_drag * V_e;

% newtons 2nd law F=ma so a=F/m
grav = [0 , 0, g]';
A_e = grav + (drag ./ m);

 x_dot = [ V_e ; A_e ];
 
end

function [zpos,isterminal,direction] = HitGround(t,x)
%INPUTS     t       scalar time value
%           x       6-element state vector
%OUTPUTS    zpos    z-coordinate of position vector
%           isterminal  Boolean flag of whether to stop integration or not
%           direction   allowable direction to approach the stop point from
%METHODOLOGY    HitGround is an Event function as required by ode45.  It is
%   checking for when zpos = 0.  Isterminal is set to 1 so integration will
%   stop when this condition is detected.  Direction is set to 0 so that
%   the function is agnostic of direction.  Essentially, HitGround will
%   detect when zpos=0 (i.e. when the projectile has hit the ground) and
%   will stop integration at that point.
zpos = x(3);
isterminal = 1;
direction = 0;
end