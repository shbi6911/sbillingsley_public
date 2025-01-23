%% assign constant values and plot initial trajectory
clear;  close all; tic;
%air density in kg/m^3 at 1655 m geopotential altitude
[rho,~,~,~,~]=stdatmo(1655,0);  
g=9.81;             %acceleration of gravity in m/s
m = 0.05;           %mass of the object in kg
Cd = 0.6;           %coefficient of drag of the object
A = pi*0.01^2;      %cross-sectional area of the object in m^2
tol = 1E-8;         %tolerance value for ode45
c = 1;              %plot counter

wind_vel = [0;0;0]; %initial wind velocity

x_ic = [0;0;0;0;20;-20];    %given initial state values (pos xyz, vel xyz)
                            % in meters and m/s

tspan = [0 20];             %time span for integration

%set options for ode45, including tolerances and stop event function
opts = odeset('Events',@HitGround,'RelTol',tol,'AbsTol',tol);

%run initial trajectory with no wind
[t_out_i,x_out_i,~,~,~] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind_vel)...
    ,tspan,x_ic,opts);

%plot initial trajectory in 3D space
figure(c); view(3);
hold on; grid on;
plot3(x_out_i(:,1),x_out_i(:,2),x_out_i(:,3),'LineWidth',2);
set(gca,'ZDir','reverse');  %reverse axes so that -z is on top
title("Initial Trajectory with Zero Wind Speed");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
hold off;   c = c+1;

%% investigate variation in wind speed
wind_vel_matrix = zeros(3,11);  %preallocate matrix of wind velocities
%establish a linear variation in horizontal (x) direction wind speed from
%0 m/s to 50 m/s
wind_vel_matrix(1,:) = linspace(0,50,11);

%loop through wind speeds and get a trajectory for each
for ii = 1:length(wind_vel_matrix)
    %trajectories and time values are stored in separate structs
    %field names for struct are t_out.trajectory1, etc.
    %trajectory1 = 0 m/s wind speed; trajectory11 = 50 m/s wind speed
    fieldname = "trajectory" + string(ii);
    wind_vel = wind_vel_matrix(:,ii);
    [t_out.(fieldname),x_out.(fieldname),~,~,~] = ode45(@(t,x) ...
        objectEOM(t,x,rho,Cd,A,m,g,wind_vel),tspan,x_ic,opts);
end

%plot all trajectories on the same 3D plot
figure(c); view(3);
set(gca,'ZDir','reverse');
hold on; grid on;
for jj = 1:length(wind_vel_matrix)
    fieldname = "trajectory" + string(jj);
    plot3(x_out.(fieldname)(:,1),x_out.(fieldname)(:,2),...
        x_out.(fieldname)(:,3));
end
title("Trajectories with Wind Speed Variation");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
hold off;   c = c+1;

%% evaluate wind speed effects on landing positions

%preallocate matrix of landing position vectors
landing_pos = zeros(3,length(wind_vel_matrix));
%extract final landing position of each trajectory
for kk = 1:length(wind_vel_matrix)
    fieldname = "trajectory" + string(kk);
    landing_pos(1:3,kk) = x_out.(fieldname)(end,1:3);
end

%determine horizontal displacement of landing and relate to wind speed
horz_disp = landing_pos(1,:)./wind_vel_matrix(1,:);
%plot 
figure(c); hold on; grid on;
plot(wind_vel_matrix(1,:),horz_disp);
title("Horizontal Displacement with Variation in Wind Speed");
xlabel ("Wind Speed (m/s)");
ylabel ("Displacement in X-Direction per m/s of Wind Speed (s)");
hold off;   c = c+1;

%determine total displacement and relate to wind speed
total_disp = vecnorm(landing_pos)./wind_vel_matrix(1,:);
%plot
figure(c); hold on; grid on;
plot(wind_vel_matrix(1,:),total_disp);
title("Total Displacement with Variation in Wind Speed");
xlabel ("Wind Speed (m/s)");
ylabel ("Total Displacement per m/s of Wind Speed (s)");
hold off;   c = c+1;

%% investigate variation in altitude and wind speed

%set a vector of altitudes from sea level to 4500 m in 500 m increments
altitude_vec = 0:500:4500;

%get air density at these altitudes
[rho_vec,~,~,~] = stdatmo(altitude_vec,0);

%preallocate a matrix to store results
wind_altitude_disp = zeros(length(wind_vel_matrix),length(altitude_vec));

%loop through vectors of altitude and wind velocity in nested for loops.
%Generate a trajectory for each combination, and determine the final
%displacement by norming the final position vector
for mm = 1:length(rho_vec)
    rho = rho_vec(mm);
    for nn = 1:length(wind_vel_matrix)
        wind_vel = wind_vel_matrix(:,nn);
        [~,traj_temp,~,~,~] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind_vel)...
        ,tspan,x_ic,opts);
        wind_altitude_disp(nn,mm) = norm(traj_temp(end,1:3));
    end
end

%plot results
figure(c);  hold on;    grid on;
for oo = 1:length(altitude_vec)
    plot(wind_vel_matrix(1,:),wind_altitude_disp(:,oo));
end
title("Displacement with Variation in Wind Speed and Altitude");
xlabel("Wind Speed (m/s)");
ylabel("Total Displacement at a Certain Altitude (m)");
legend(string(altitude_vec),'Location','eastoutside');
hold off; c = c+1;

%% investigate variations in kinetic energy

%find kinetic energy of given initial velocity vector, in joules
ke_initial = 0.5*m*(norm(x_ic(4:6)))^2;

%create vector of linear variation in mass values between 10g and 1kg
m_vec = linspace(0.01,1,100);

h_vec = x_ic(4:6)./norm(x_ic(4:6)); %unit heading vector at launch
%preallocate matrix of initial velocity vectors
ke_variation = zeros(3,length(m_vec));
ke_variation = h_vec + ke_variation;  %all vectors have same heading
%create vector of velocity magnitudes corresponding to given kinetic energy
%and previously created vector of masses
v_mag = sqrt((2*ke_initial)./m_vec);
%create matrix of initial velocity vectors
ke_variation = v_mag.*ke_variation;

%preallocate a vector of landing distances corresponding to masses
ke_dist = zeros(1,length(m_vec));

wind_vel = [0;0;0]; %reset wind velocity to zero

%loop through velocity vectors, find a trajectory for each one, determine
%the total displacement at landing, and store it
for ll = 1:length(m_vec)
    x_ic = [0;0;0;ke_variation(:,ll)];  %initial conditions
    [~,trajectory,~,~,~] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m_vec(ll),...
        g,wind_vel),tspan,x_ic,opts);  %get trajectory
    ke_dist(ll) = norm(trajectory(end,1:3));    %store landing displacement
end

%plot mass vs. landing displacement
figure(c);
hold on; grid on;
plot(m_vec,ke_dist);
title("Total Displacement with Variation in Mass and Constant KE");
xlabel("Mass (kg)");
ylabel("Total Displacement (m)");
hold off; c = c+1;

disp("Elapsed Time = " + string(toc));
%% functions

function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind_vel)
%
%NOTE:  All vector inputs/outputs should use NED axes
%
%INPUTS     t   scalar value of current integration time
%           x   six-element state vector (position x;y;z;velocity x;y;z)
%           rho scalar value of local air density
%           Cd  scalar value of object drag coefficient
%           A   scalar value of object cross-sectional area
%           m   scalar value of object mass
%           g   scalar value of local gravitational acceleration
%           wind_vel    three element wind velocity vector (x;y;z)
%
%OUTPUTS    xdot    six-element derivative state vector
%
%METHODOLOGY    objectEOM creates a derivative state vector for an object
%   moving through space under Newton's 2nd Law, assuming gravity and drag
%   as the only accelerations.  Wind, air density, gravity and the physical
%   properties of the object are all assumed constant.

a_grav = [0;0;g];   %define a gravitational acceleration vector

v_e = x(4:6) - wind_vel;    %find air relative velocity vector
v_a = norm(v_e);            %find airspeed
D = 0.5*Cd*A*rho*(v_a^2);   %calculate magnitude of drag force
a_drag = (-D.*(v_e./v_a))./m;   %define a drag acceleration vector

%create derivative state vector.  Velocity is derivative of position, and
%calculated acceleration is derivative of velocity
xdot = [x(4:6);(a_grav+a_drag)];

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