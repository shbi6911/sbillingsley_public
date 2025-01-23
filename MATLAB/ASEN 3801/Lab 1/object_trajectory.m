% Contributors: Shane Billingsley and Gabriel Law
% Course number: ASEN 3801
% File name: object_trajectory.m
% Created: 1/17/24

clc;
clear;
close all;
tic;

%% Question 1
%% Defining Constants
tspan = [0 20];         % Time span to be integrated over
init_con = [1;1;1;1];   % Initial conditions

%% Solving 1a
opts = odeset('RelTol',1E-8,'AbsTol',1E-8);         % Setting tolerances
[t_out,x_out] = ode45(@odefun,tspan,init_con,opts); % Running ode45 

%% Plotting 1a
f = figure();
subplot(4,1,1)
plot(t_out,x_out(:,1))
title("w vs Time")
xlabel("Time (n.d.)")
ylabel("w (n.d.)")
ylim([min(x_out(:,1))-0.1 1])

subplot(4,1,2)
plot(t_out,x_out(:,2))
title("x vs Time")
xlabel("Time (n.d.)")
ylabel("x (n.d.)")
ylim([min(x_out(:,2))-0.25 max(x_out(:,2))+0.25])

subplot(4,1,3)
plot(t_out,x_out(:,3))
title("y vs Time")
xlabel("Time (n.d.)")
ylabel("y (n.d.)")
ylim([min(x_out(:,3))-0.25 max(x_out(:,3))+0.25])

subplot(4,1,4)
plot(t_out,x_out(:,4))
title("z vs Time")
xlabel("Time (n.d.)")
ylabel("z (n.d.)")
ylim([min(x_out(:,4))-0.25 max(x_out(:,4))+0.25])

saveas(f,'P1a','png')
%% Problem 1b
tol = 1E-12;                                                % Creating tolerance variable
optsR = odeset('RelTol',tol,'AbsTol',tol);                  % Creating options for reference values
[t_outR,x_outR] = ode45(@odefun,tspan,init_con,optsR);      % Running reference ode45
final_values = zeros(4, 5);                                 % Pre-allocating the final output matrix

for ii = 1:5                                                    % Looping through the other cases
    tol = tol/((1E-2) * ii);                                    % Getting tolerance value
    opts = odeset('RelTol',tol,'AbsTol',tol);                   % Creating options
    [t_out,x_out] = ode45(@odefun,tspan,init_con,opts);         % running ode45
    final_values(:,6-ii) = (abs(x_out(end,:)-x_outR(end,:)))';  % Calculation of difference and placement into final matrix
end

format shortE
disp(final_values)

%% Question 2
%% assign constant values and plot initial trajectory
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
%run trajectory with gravity only for comparison
[t_out_basic,x_out_basic,~,~,~] = ode45(@basic_EOM,tspan,x_ic,opts);

%plot initial trajectory in 3D space
figure(); view(3);
hold on; grid on;
plot3(x_out_i(:,1),x_out_i(:,2),x_out_i(:,3),'LineWidth',2);
plot3(x_out_basic(:,1),x_out_basic(:,2),x_out_basic(:,3),'LineWidth',2);
set(gca,'ZDir','reverse');  %reverse axes so that -z is on top
title("Initial Trajectory with Zero Wind Speed");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
legend("Gravity and Drag","Gravity Only",'location','best');
hold off;

%% investigate variation in wind speed
wind_vel_matrix = zeros(3,101);  %preallocate matrix of wind velocities
%establish a linear variation in horizontal (x) direction wind speed from
%0 m/s to 50 m/s
wind_vel_matrix(1,:) = linspace(0,50,length(wind_vel_matrix));

%loop through wind speeds and get a trajectory for each
for ii = 1:length(wind_vel_matrix)
    %trajectories and time values are stored in separate structs
    %field names for struct are t_out_w.trajectory1, etc.
    %trajectory1 = 0 m/s wind speed; trajectory101 = 50 m/s wind speed
    fieldname = "trajectory" + string(ii);
    wind_vel = wind_vel_matrix(:,ii);
    [t_out_w.(fieldname),x_out_w.(fieldname),~,~,~] = ode45(@(t,x) ...
        objectEOM(t,x,rho,Cd,A,m,g,wind_vel),tspan,x_ic,opts);
end

%plot all trajectories on the same 3D plot
figure(); view(3);
set(gca,'ZDir','reverse');
hold on; grid on;
for jj = 1:length(wind_vel_matrix)
    if mod((jj-1),10) == 0
        fieldname = "trajectory" + string(jj);
        plot3(x_out_w.(fieldname)(:,1),x_out_w.(fieldname)(:,2),...
            x_out_w.(fieldname)(:,3),'DisplayName',...
            ("Wind Speed " + string(wind_vel_matrix(1,jj)) + " m/s"));
    end
end
title("Trajectories with Wind Speed Variation");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
legend show
hold off;

%% evaluate wind speed effects on landing positions

%preallocate matrix of landing position vectors
landing_pos = zeros(3,length(wind_vel_matrix));
%extract final landing position of each trajectory
for kk = 1:length(wind_vel_matrix)
    fieldname = "trajectory" + string(kk);
    landing_pos(1:3,kk) = x_out_w.(fieldname)(end,1:3);
end

%determine horizontal displacement of landing and relate to wind speed
horz_disp = landing_pos(1,:);
%plot 
figure(); hold on; grid on;
plot(wind_vel_matrix(1,:),horz_disp);
title("Horizontal Displacement with Variation in Wind Speed");
xlabel ("Wind Speed (m/s)");
ylabel ("Displacement in X-Direction (m)");
hold off;

%determine total displacement and relate to wind speed
total_disp = vecnorm(landing_pos);
%plot
figure(); hold on; grid on;
plot(wind_vel_matrix(1,:),total_disp);
title("Total Displacement with Variation in Wind Speed");
xlabel ("Wind Speed (m/s)");
ylabel ("Total Displacement (m)");
hold off;

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

%find minimum distance for a given altitude
disp_minimum = min(wind_altitude_disp);

%plot results
figure();  hold on;    grid on;
for oo = 1:length(altitude_vec)
    plot(wind_vel_matrix(1,:),wind_altitude_disp(:,oo),'DisplayName',...
        (string(altitude_vec(oo)) + " m Altitude"));
end
title("Displacement with Variation in Wind Speed and Altitude");
xlabel("Wind Speed (m/s)");
ylabel("Total Displacement at a Certain Altitude (m)");
legend('Location','eastoutside');
hold off;

figure();  hold on;    grid on;
plot(altitude_vec,disp_minimum)
title("Minimum Displacement (for given wind speed) vs. Altitude");
xlabel("Altitude (m)");
ylabel("Minimum Displacement at Given Wind Speed (m)");
hold off;

%% investigate variations in kinetic energy

%find kinetic energy of given initial velocity vector, in joules
ke_initial = 0.5*m*(norm(x_ic(4:6)))^2;

%create vector of linear variation in mass values between 0.1g and 50 g
m_vec = linspace(0.001,0.05,100);

h_vec = x_ic(4:6)./norm(x_ic(4:6)); %unit heading vector at launch
%preallocate matrix of initial velocity vectors
ke_variation = zeros(3,length(m_vec));
ke_variation = h_vec + ke_variation;  %all vectors have same heading
%create vector of velocity magnitudes corresponding to given kinetic energy
%and previously created vector of masses
v_mag = sqrt((2*ke_initial)./m_vec);
%create matrix of initial velocity vectors
ke_variation = v_mag.*ke_variation;

ke_wind_vel_matrix = zeros(3,11);  %preallocate matrix of wind velocities
%establish a linear variation in horizontal (x) direction wind speed from
%0 m/s to 50 m/s
ke_wind_vel_matrix(1,:) = linspace(0,50,length(ke_wind_vel_matrix));

%preallocate a vector of landing distances corresponding to masses
ke_dist = zeros(length(ke_wind_vel_matrix),length(m_vec));

%loop through velocity vectors and wind speeds, find a trajectory for each 
% one, determine the total displacement at landing, and store it
for bb = 1:length(ke_wind_vel_matrix)
    wind_vel = ke_wind_vel_matrix(:,bb);
    for qq = 1:length(m_vec)
        x_ic_ke = [0;0;0;ke_variation(:,qq)];  %initial conditions
        [~,ke_traj,~,~,~] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m_vec(qq),...
            g,wind_vel),tspan,x_ic_ke,opts);  %get trajectory
        ke_dist(bb,qq) = norm(ke_traj(end,1:3));    %store landing displacement
    end
end

%plot mass vs. landing displacement
figure();
hold on; grid on;
for cc = 1:length(ke_wind_vel_matrix)
    plot(m_vec,ke_dist(cc,:),'DisplayName',...
        ("Wind Speed " + string(ke_wind_vel_matrix(1,cc)) + " m/s"));
end
title("Total Displacement with Variation in Mass and Constant KE");
xlabel("Mass (kg)");
ylabel("Total Displacement (m)");
legend ('Location','eastoutside');
hold off;

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

function [y_prime] = odefun(t,y_vec,var)
%INPUTS     t   scalar time value given by ode45
%           y   4-element state vector
%
%OUTPUTS    y_prime     4-element derivative state vector
%
%METHODOLOGY    odefun takes in a state vector and outputs the derivate 
% of that state vector according to given equations

w = y_vec(1);   x=y_vec(2);     y=y_vec(3);     z=y_vec(4);

y_prime = [0;0;0;0];
y_prime(1)= (-9*w) + y;
y_prime(2)= (4*w*x*y) - x^2;
y_prime(3)= (2*w) - x -(2*z);
y_prime(4)= (x*y) - y^2 - (3*z^3);

end

function x_out = basic_EOM(t,x)
%NOTE NED AXES
% INPUTS:       t   scalar time value (from ode45)
%               x   input state vector (position x;y;z; velocity x;y;z)
%OUTPUTS:       x_out   derivative state vector (velocity x;y;z; accel
%                       x;y;z)
%METHODOLOGY    This is an EOM function to feed ode45 that solely takes
%               into account acceleration due to gravity

g= 9.81;    a_grav = [0;0;g];
p_prime = x(4:6);
v_prime = a_grav;
x_out = [p_prime;v_prime];
end
