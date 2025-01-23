%% ASEN 2012: Project 2: Bottle Rocket Trajectory
% by:               [Shane Billingsley & Margot Daberkow]
% SID:              [110231742]
% date created:     [11/11/22]
% last modified:    [12/07/22]
%
% PURPOSE: 
% The purpose of this script is to integrate the trajectory of a water
% bottle rocket, using Newton's Laws of Motion, and equations governing the
% pressure and volume of air within the rocket.  It can also define a
% parameter space and examine all trajectories within that space.
%
% To accomplish this, the script defines a function to feed into ode45 for
% integration, and then plots the integrated data vs. a given verification
% case.
% 
% INPUTS:  
% The code requires a verification case for comparison.  All other input
% values are set using the getConst and getParams functions.
% OUTPUTS:
% The code outputs plots of a set of bottle rocket trajectories, defined
% in terms of position and velocity in the x and z directions.  Rocket
% thrust vs. time can also be plotted.  The script contains an n-value that
% defines the parameter space.  n^4 trajectories are created in each run.
% ASSUMPTIONS:
% We assume that atmospheric conditions at launch can be treated as
% constant (i.e. the rocket does not go high enough for these values to
% change.)
% We assume that forces in the y-direction are negligible for this
% calculation.
% We assume that the drag coefficient of the rocket is constant (i.e. it
% does not change with orientation of the rocket.)
% We assume that the expansion of air in the rocket, and its escape through
% the throat, is isentropic.

%% NOTE:  For simplicity in coding, in order to examine the variation of an
%individual parameter with respect to the others, certain lines of code
%must be commented in and out, and certain values in the code must be changed.
%This will be noted in the code.


%% housekeeping I
clc;
clear all;
tic                                                         %start time

%% Main script to calculate trajectory

% import comparison data
% load ("project2verification.mat");

% set number of data points for varying parameters
% if n=1, will run a default trajectory (currently set to the 85m
% trajectory
%% NOTE: n-values greater than 4 will take considerable time to run
n=1;
count = 0; %counting variable to keep track of loops
%set time span for integration
tspan = [0 10];

%preallocate arrays for results
% params_reference = zeros(n,5);           %comment in for individual parameters
% distances = zeros(n,3);                     
distances = zeros((n^4),6);                %use for all parameters
params_reference = zeros(n^4,5);
a=1; b=1; c=1; d=1;     %initialize loop variables (needed for ind. params.)

%% Run all trajectories in a nested for loop structure
%NOTE: if it is desired to examine an individual parameter while holding
%others constant, comment out all but one for loop.  Remember to comment
%out the end statements as well.
% for a = 1:n                                             %air pressure
    for b = 1:n                                         %water volume
%         for c = 1:n                                     %drag coefficient
%             for d = 1:n                                 %launch angle
                ind = [a,b,c,d]; %define index array for this trajectory
                count = count +1; %iterate count variable

                %define initial conditions for integration
                const = getConst(ind,n); %define constants for main script
                %set initial conditions
                IC = [const.p0;const.v0;const.V0;const.m_water;const.m_air];
                %state vector contains 9 variables [position(x,y,z),velocity(x,y,z),air
                %volume, water mass, air mass]
                
                %integrate trajectory
                [t,y,data] = getTrajectory(tspan,IC,ind,n);

                %store data from this trajectory in output structures
                trajectories.("trajectory"+string(count))=y;
                trajectories.("time"+string(count))=t;
                trajectories.("data"+string(count))=data;
                %record the parameters of this trajectory for reference
                params_reference(count,:) = [count,(const.P0-const.Pa),const.Vw,...
                    const.C_drag,const.theta];
%             end
%         end
    end
% end

%% for individual parameters, build array of results

% for i = 1:n
%     distances(i,1) = params_reference(i,2);   %alter second index per parameter
%     distances(i,2) = max(trajectories.("trajectory"+string(i))(:,1));
%     distances(i,3) = max(trajectories.("trajectory"+string(i))(:,3));
% end

%% for all parameters, extract relevant data

%loop through stored trajectories and extract relevant values
% for i = 1:(n^4)
%     distances(i,1) = i;     %index
%     distances(i,2) = max(trajectories.("trajectory"+string(i))(:,1));
%     %maximum distance reached in trajectory (needed to find 85m target
%     %trajectory)
%     distances(i,3:6) = params_reference(i,2:5); 
%     %parameter values which created this trajectory
% end

%sort resulting array by distance
% dist_sort = sortrows(distances,2);

%create array of all trajectories within +/- 0.5 m of 85 m
% for i = 1:(n^4)
%     if dist_sort(i,2) >= 84.5 & dist_sort(i,2) <= 85.5
%         dist_targ(i,:) = dist_sort(i,:);
%     end
% end

target_trajectories = reshape(nonzeros(dist_targ),[],6);
%convert pressure data back to psi
target_trajectories(:,3) = target_trajectories(:,3) ./ (101325/14.7);


%% Plot results
%NOTE:  This section contains many different plots for many situations.  It
%is intended that plots not being used in the current run be commented out.

%NOTE:  The single trajectory plots plot a given trajectory from the final
% "trajectories" data structure.  Change this number to plot a different
% trajectory.  If the n-value is low enough that the trajectory being plotted is
% not being created, this will cause an error.

% figure(1);                                  % plot x-position
% hold on; grid on;                           % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,1),"LineWidth",2);       
% scatter(verification.time,verification.distance,5,"red","filled","o");
% title ("Bottle Rocket Trajectory (X)");
% xlabel ("Time (s)");
% ylabel ("X-Position (m)");
% xlim([0 5]);
% legend ("My Data","Verification");
% hold off;

% figure(2);                                % plot y-position
% hold on; grid on;                         % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,2),"LineWidth",2);       
% title ("Bottle Rocket Trajectory (Y)");
% xlabel ("Time (s)");
% ylabel ("Y-Position (m)");
% hold off;

% figure(3);                                  % plot z-position
% hold on; grid on;                           % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,3),"LineWidth",2);       
% scatter(verification.time,verification.height,5,"red","filled","o");
% title ("Bottle Rocket Trajectory (Z)");
% xlabel ("Time (s)");
% ylabel ("Z-Position (m)");
% xlim([0 5]);
% legend ("My Data","Verification");
% hold off;

% figure(4);                                  % plot x-velocity
% hold on; grid on;
% plot(trajectories.time1,trajectories.trajectory1(:,4),"LineWidth",2);       
% scatter(verification.time,verification.xVelocity,5,"red","filled","o");
% title ("Bottle Rocket Velocity (X)");
% xlabel ("Time (s)");
% ylabel ("X-Velocity (m/s)");
% xlim([0 4]);
% legend ("My Data","Verification");
% hold off;

% figure(5);                                % plot y-velocity
% hold on; grid on;                         % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,5),"LineWidth",2);       
% title ("Bottle Rocket Velocity (Y)");
% xlabel ("Time (s)");
% ylabel ("Y-Velocity (m/s)");
% hold off;

% figure(6);                                  % plot z-velocity
% hold on; grid on;                           % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,6),"LineWidth",2);       
% scatter(verification.time,verification.zVelocity,5,"red","filled","o");
% title ("Bottle Rocket Velocity (Z)");
% xlabel ("Time (s)");
% ylabel ("Z-Velocity (m/s)");
% xlim([0 4]);
% legend ("My Data","Verification");
% hold off;

% figure(7);                                  % plot air volume
% hold on; grid on;                           % single trajectory
% plot(trajectories.time1,trajectories.trajectory1(:,7),"LineWidth",2);      
% scatter(verification.time,verification.airVolume,5,"red","filled","o");
% title ("Bottle Rocket Air Volume (m^3)");
% xlabel ("Time (s)");
% ylabel ("Volume (kg/m^3)");
% xlim([0 0.3]);
% legend ("My Data","Verification");
% hold off;

figure(8);                                  % plot thrust
hold on; grid on;                           % single trajectory
plot(trajectories.time1,trajectories.data1(:,1),"LineWidth",2);      
scatter(verification.time,verification.thrust,5,"red","filled","o");
title ("Bottle Rocket Thrust");
xlabel ("Time (s)");
ylabel ("Thrust (N)");
xlim([0 0.3]);
legend ("My Data","Verification");
hold off;

% thrust = trajectories.data1(:,1);        % sort thrust data by phase
% flag = trajectories.data1(:,4);
% phase1 = thrust(flag == 10); t1 = trajectories.time1(flag == 10);
% phase21 = thrust(flag == 21); t21 = trajectories.time1(flag == 21);
% phase22 = thrust(flag == 22); t22 = trajectories.time1(flag == 22);
% phase3 = thrust(flag == 30); t3 = trajectories.time1(flag == 30);
% phase4 = thrust(flag == 40); t4 = trajectories.time1(flag == 40);
% 
% figure(9);                                  % plot thrust by phase
% hold on; grid on;                           % single trajectory
% plot(t1,phase1,"LineWidth",2,"Color",'r');
% plot(t21,phase21,"LineWidth",2,"Color",'g');
% plot(t22,phase22,"LineWidth",2,"Color",'b');
% plot(t3,phase3,"LineWidth",2,"Color",'y');
% plot(t4,phase4,"LineWidth",2,"Color",'k');
% title ("Bottle Rocket Thrust");
% xlabel ("Time (s)");
% ylabel ("Thrust (N)");
% xlim([0 0.3]);
% legend ("Phase 1","Phase 2 (choked)","Phase 2 (unchoked)", ...
%     "Phase 3 (ballistic)","Phase 4 (landed)");
% hold off;

figure(10);                                 % plot pressure
hold on; grid on;                           % single trajectory
plot(trajectories.time1,trajectories.data1(:,2),'LineWidth',2);       
title ("Bottle Rocket Air Pressure");
xlabel ("Time (s)");
ylabel ("Pressure (Pa)");
xlim([0 0.3]);
hold off;

% figure (11);                                  %plot launch angle vs dist
% hold on; grid on;             %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,2),'LineWidth',2);
% title ("Launch Angle vs. Distance");
% xlabel ("Launch Angle (radians)");
% ylabel ("Max Distance (m)");
% xticks([(pi/6),(pi/4),(pi/3)]);
% xticklabels({'pi/6','pi/4','pi/3'});
% 
% figure (12);                                 %plot launch angle vs height
% hold on; grid on;              %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,3),'LineWidth',2);
% title ("Launch Angle vs. Height");
% xlabel ("Launch Angle (radians)");
% ylabel ("Max Height (m)");
% xticks([(pi/6),(pi/4),(pi/3)]);
% xticklabels({'pi/6','pi/4','pi/3'});

% figure (13);                                  %plot drag coeff vs dist
% hold on; grid on;             %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,2),'LineWidth',2);
% title ("Drag Coefficient vs. Distance");
% xlabel ("Drag Coefficient");
% ylabel ("Max Distance (m)");
% 
% figure (14);                                  %plot drag coeff vs height
% hold on; grid on;             %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,3),'LineWidth',2);
% title ("Drag Coefficient vs. Height");
% xlabel ("Drag Coefficient");
% ylabel ("Max Height (m)");

% figure (15);                                  %plot water volume vs dist
% hold on; grid on;             %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,2),'LineWidth',2);
% title ("Water Volume vs. Distance");
% xlabel ("Water Volume (m^3)");
% ylabel ("Max Distance (m)");
% 
% figure (16);                                  %plot water vol vs height
% hold on; grid on;             %single parameter, multiple trajectories
% plot(distances(:,1),distances(:,3),'LineWidth',2);
% title ("Water Volume vs. Height");
% xlabel ("Water Volume (m^3)");
% ylabel ("Max Height (m)");

% figure (17);                                  %plot pressure vs dist
% hold on; grid on;             %single parameter, multiple trajectories
% plot((distances(:,1)./(101325/14.7)),distances(:,2),'LineWidth',2);
% title ("Air Pressure vs. Distance");
% xlabel ("Air Pressure (psi)");
% ylabel ("Max Distance (m)");
% 
% figure (18);                                  %plot pressure vs height
% hold on; grid on;             %single parameter, multiple trajectories
% plot((distances(:,1)./(101325/14.7)),distances(:,3),'LineWidth',2);
% title ("Air Pressure vs. Height");
% xlabel ("Air Pressure (psi)");
% ylabel ("Max Height (m)");

%output final comparison values to the command window
% disp("Max Altitude (My Data)");
% disp(max(trajectories.trajectory1(:,3)));
% disp("Max Altitude (Verification)");
% disp(max(verification.height))
% disp ("Max Distance (My Data)");
% disp (max(trajectories.trajectory1(:,1)));
% disp ("Max Distance (Verification)");
% disp (max(verification.distance))
%% housekeeping II
elapsedTime = toc                                              %end time
save("trajectories");           %save trajectory data to avoid rerunning

%% Subfunction getTrajectory.  This integrates a trajectory.
function [t,y,data] = getTrajectory(tspan,IC,ind,n)
% 
% Inputs:   tspan = span of time for integration
%           IC = initial state vector
%           ind = an index corresponding to constant values for this
%           particular trajectory, defined by getParams
%           n = input parameter for getParams defining number of data
%           points for each parameter range
%                
% Outputs:  t = time values from ode45 integration
%           y = function values from ode45 integration
%           data = misc data values not returned by ode45
%
% Methodology: this function calls ode45 to integrate a trajectory using
% the rocketEOM function.  By calling ode45 from this function, we can
% integrate multiple trajectories with changes in initial values or other
% constants.
%
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);  %set tolerances for ode45

%integrate trajectory using ode45
[t,y] = ode45(@(t,X) rocketEOM(t,X,ind,n),tspan,IC,opts);

% extract additional data from rocketEOM function
% most importantly thrust, but also other values for debugging purposes
data = zeros(length(t),4);
for i = 1:length(t)
    [~,data(i,:)]=rocketEOM(t(i),y(i,:),ind,n);
end

end

%% Subfunction rocketEOM. This is the subfunction called by ODE45. 
function [dX,data] = rocketEOM(t,X,ind,n)
% 
% Inputs:   t = anonymous time variable
%           X = anonymous state vector containing
%             = [p_x;p_y;p_z;v_x;v_y;v_z;air volume;mass of water, mass of air]
%           ind = a 4-element vector for parameters matrix (see getParams)
%           n = a scalar value for number of data points in each parameter
%           range (see getParams) 
%                
% Outputs: dX = derivative state vector containing
%             = [dp/dx;dp/dy;dp/dz;dv/dx;dv/dy;dv/dz;dV/dt;dm_w/dt;dm_a/dt]
%           data = a vector of additional data for debugging
%
% Methodology: function used for integrating the flight of a water bottle
% rocket
%

% Extract state vector
p = [X(1);X(2);X(3)];   % From the variable X coming in, extract the indices 
                        % corresponding to position
v = [X(4);X(5);X(6)];   % From the variable X coming in, extract the indices 
                        % corresponding to velocity
V = X(7);   % Extract the index corresponding to volume of air
m_w = X(8);   % Extract the index corresponding to mass of water
m_a = X(9);   % Extract the index corresponding to mass of air


const = getConst(ind,n);     %define constants within the function

% set heading vector and adjust for wind
v_rel_vec = v-const.w0;
v_rel_mag = norm(v_rel_vec);
if norm(p) <= const.Ls          %rocket still on launch stand
    h = const.h0;
else
    h = v_rel_vec/v_rel_mag;    %heading vector free to change      %EQ 2
end

%determine air pressure at two different phases of flight
if V <= const.V     % before water is exhausted
    P = const.P0*((const.V0/V)^const.gamma);                        %EQ 5
else                % after water is exhausted
    P = const.P_end*((m_a)/const.m_air)^const.gamma;                %EQ 16
end

% Determine the phase of flight of the system to write the proper equations of motion 
% using if statements to check the relevant state variables

if (V <= const.V) && (P >= const.Pa)  % check that volume of air is less than volume of the
                  % bottle rocket and pressure is greater than atmospheric
    % phase 1: expelling water

    % Calculate forces for Newton's Second Law
    fGrav = [0;0;-((m_w+m_a+const.m_dry)*const.g0)];                %EQ 1
    fThrust = h*(2*const.C_dis*const.At*(P-const.Pa));              %EQ 9
    fDrag = -h*(0.5*const.Ac*const.C_drag*const.rho_air*(v_rel_mag^2)); %EQ 4
    % define change in mass
    dm_wdt = -(const.C_dis*const.At*sqrt(2*const.rho_water*(P-const.Pa))); %EQ 6
    dm_adt = 0;
    % define change in volume
    dVdt = const.C_dis*const.At*...                                 %EQ 10
        sqrt((2/const.rho_water)*(const.P0*((const.V0/V)^const.gamma)-const.Pa));
    % define values for debugging purposes
    M_e = 0;
    flag = 10;

elseif P >= const.Pa % check that air pressure is greater than ambient
                     % air pressure
    % phase 2: expelling air
  
    % Calculate forces for Newton's Second Law
    fGrav = [0;0;-((m_w+m_a+const.m_dry)*const.g0)];                %EQ 1
    fDrag = -h*(0.5*const.Ac*const.C_drag*const.rho_air*(v_rel_mag^2)); %EQ 4

    %define critical pressure
    P_crit = P*((2/(const.gamma+1))^(const.gamma/(const.gamma-1))); %EQ 17
    % Check for choked flow through the throat
    if P_crit > const.Pa     %choked flow
        rho = m_a/const.V;  %find current air density in bottle     %EQ 18
        T = P/(rho*const.R_air);    %find current temperature in bottle %EQ 18
        P_e = P_crit;    %define pressure at the exit               %EQ 18
        M_e = 1;         %define Mach number at the exit
        T_e = (2/(const.gamma+1)*T); %define temperature at the exit    %EQ 18
        rho_e = P_e/(const.R_air*T_e);  %define density at the exit %EQ 18
        V_e = M_e*sqrt(const.gamma*const.R_air*T_e); %define exit velocity %EQ 18
        % define values for debugging purposes
        flag = 21;

    else                    %not choked flow
        rho = m_a/const.V;  %find current air density in bottle     %EQ 19
        T = P/(rho*const.R_air);    %find current temperature in bottle %EQ 19
        P_e = const.Pa;    %define pressure at the exit             %EQ 19
                           %define Mach number at the exit
        M_e = sqrt((((P/const.Pa)^((const.gamma-1)/const.gamma))-1)*...
            (2/(const.gamma-1)));                                   %EQ 19
                           %define temperature at the exit
        T_e = T/((((const.gamma-1)/2)*M_e^2)+1);                    %EQ 19
        rho_e = const.Pa/(const.R_air*T_e);  %define density at the exit %EQ 19
        V_e = M_e*sqrt(const.gamma*const.R_air*T_e); %define exit velocity %EQ 19
        % define values for debugging purposes
        flag = 22;
    end
    % define change in mass
    dm_adt = -const.C_dis*rho_e*const.At*V_e;                        %EQ 21
    dm_wdt = 0;

    fThrust = h*(((-dm_adt)*V_e)+((P_e-const.Pa)*const.At));        %EQ 20

    % define change in volume
    dVdt = 0;

elseif p(3) >= 0   % check that the rocket is above the ground 
    % phase 3: ballistic phase

    % Calculate forces for Newton's Second Law
    fGrav = [0;0;-((m_w+m_a+const.m_dry)*const.g0)];                %EQ 1
    fThrust = [0;0;0];
    fDrag = -h*(0.5*const.Ac*const.C_drag*const.rho_air*(v_rel_mag^2)); %EQ 4
    
    % define change in mass
    dm_wdt = 0;
    dm_adt = 0;

    % define change in volume
    dVdt = 0;
    
    % define values for debugging purposes
    M_e = 0;
    flag = 30;

else                % after rocket has landed, stop any changes
    fThrust = [0;0;0];
    fGrav = [0;0;0];
    fDrag = [0;0;0];
    v = [0;0;0];
    dm_adt = 0;
    dm_wdt = 0;
    dVdt = 0;
    M_e = 0;
    flag = 40;
end

% expression for net force and acceleration per Newton's Second Law
fnet = fThrust + fGrav + fDrag;                                  %EQ 1
%forces are signed properly in previous steps

a = fnet/(m_a+m_w+const.m_dry); % calculate acceleration

dX = [v;a;dVdt;dm_wdt;dm_adt]; %create output derivative vector using changed variables
data = [norm(fThrust),P,M_e,flag]; %output thrust data and debugging values
end
%% Subfunction getConst which sets up the constants needed to solve this problem. 
function const = getConst(ind,n)
% 
% Inputs:   ind = a 4-element vector for parameters matrix (see getParams)
%           n = a scalar value for number of data points in each parameter
%           range (see getParams)
%                
% Outputs: const = structure containing initial conditions, physical
% constants, given rocket characteristics, and other values not changing
% with time.
%
% Methodology: function used to define a constant structure for ode45 integration
% 

% physical constants
const.g0 = 9.81; % gravitational acceleration (m/s^2)
const.Pa = 83403.57; % ambient air pressure at launch altitude (Pa) (N/m^2)
const.rho_air = 0.961; % density of ambient air at launch altitude (kg/m^3)
const.rho_water = 1000; % density of water (kg/m^3)
const.gamma = 1.4; % ratio of specific heats for air
const.R_air = 287; % R-value for air (J/kg*K)

% define parameters of this individual trajectory
params = getParams(n);                  % define parameter matrix
const.P0 = params(1,ind(1))+const.Pa; % gage pressure+ambient pressure (Pa) (N/m^2)
const.Vw = params(2,ind(2)); %initial volume of water (m^3)
const.C_drag = params(3,ind(3)); % drag coefficient of the rocket
const.theta = params(4,ind(4)); % define launch angle (rad)
const.h0 = [cos(const.theta);0;sin(const.theta)]; % initial heading vector (m)

% initial conditions
const.p0 = [0;0;0.25]; % initial position vector (m)
const.v0 = [0;0;0]; % initial velocity vector (m/s)
const.w0 = [0;0;0]; % initial wind vector (m/s)
const.T0 =  300; % initial temperature (K)
const.Ls = 0.5; % length of launch stand (m)

% properties of the bottle rocket
const.m_dry =  0.15; % dry mass of bottle rocket (kg)
const.d =  0.105; % diameter of the bottle rocket (m)
const.d_t = 0.021; % diameter of the rocket throat (m)
const.V = 0.002; % volume of empty bottle rocket (m^3)
const.C_dis = 0.8; % discharge coefficient of the rocket

% some derived properties
const.V0 = const.V-const.Vw; %initial volume of air
const.m_air = ...   % initial mass of air in the bottle rocket (kg)
    ((const.P0)*const.V0)/(const.R_air*const.T0);
const.m_water = ... % initial mass of water in the bottle rocket (kg)
    (const.rho_water*const.Vw);
const.m_wet = ...   % initial wet mass of the bottle rocket (kg)
    const.m_dry + const.m_air + const.m_water;                      %EQ 12
const.P_end = ...   % air pressure when water is exhausted (Pa) (N/m^2)
    const.P0*((const.V0/const.V)^const.gamma);                      %EQ 14
const.T_end = ...   % air temperature when water is exhausted (K)
    const.T0*((const.V0/const.V)^(const.gamma-1));                  %EQ 15
const.Ac = pi * const.d^2 * (1/4);  % cross-sectional area of rocket (m^2)
const.At = pi * const.d_t^2 * (1/4); %cross-sectional area of the throat (m^2)

end

%% Subfunction getParams.  This defines a paramater matrix for varying trajectories.
function params = getParams(n)
% 
% Inputs:  n = a scalar value for number of data points in each parameter
% range
%                
% Outputs: params = matrix containing a range of parameters for a trajectory
%
% Methodology: this function defines a matrix containing a range of values
% for each of four parameters for a trajectory.  These parameters include
% initial air pressure, initial volume of water, drag coefficient, and
% launch angle.
%
if n == 1
    params = [61.25*(101325/14.7);0.0005;0.425;((3*pi)/16)]; %define default case
else
    %set parameter ranges according to n-value
    %NOTE: to examine one parameter with the others held constant, change
    %the range of the constant parameters to have the same start and end
    %value.  It is recommended to comment out some for loops when doing
    %this.  See the section "Run all trajectories in a nested for loop
    %structure."
    airpressure = linspace(40,40,n);
    airpressure = airpressure .* (101325/14.7);
    watervolume = linspace(0.0005,0.0015,n);
    dragcoefficient = linspace(0.5,0.5,n);
    launchangle = linspace((pi/4),(pi/4),n);
    %combine parameters into a single 4 x n matrix
    params = [airpressure;watervolume;dragcoefficient;launchangle];
end

end