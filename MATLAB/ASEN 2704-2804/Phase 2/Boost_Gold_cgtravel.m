%% ASEN 2804: Vehicle Design Lab
%% Boost Glider Model
% Authors: Alex Virga, Christian Bowman. Ethan Davis, Murilo Tibana, 
% Nick Ratajczyk, Sam Allen, Sebastian Escobar Gonzalez, Shane Billingsley

% housekeeping
clear; clc; close all

b = 0.53109; % [m] Wing span
S = 0.111; % [m^2] Planform area
AR = 1.875; % Aspect Ratio
Re_x_cr = 5*10^5; % Critical Reynold's Number - Assume flat plate
M = 0.3; % Mach number
cr = 35.719 * 10^-2; % Root chord
ct = 7.138 * 10^-2; % Tip chord
c_mean = (cr+ct)/2; % Mean chord
rho = 1.225; % [kg/m^3] Assume sea level, STD ATM
mu = 1.79e-5; % [kg/m*s] Assume sea level, STD ATM
T = 288.16; % [K] Assume sea level, STD ATM
a = sqrt(1.4*287*T); % [m/s] Assume sea level, STD ATM
V = M * a; % [m/s] V infinity
x_cr = (mu* Re_x_cr) / (rho * V); % [m] Critical point in the wing
Re_L = c_mean*rho*V/mu; % Reynold's Number
turb_skin_fric = 0.074 / (Re_L ^ 0.2); % Laminar skin friction
lam_skin_fric = 1.328 / sqrt(Re_L); % Turbulent skin friction

% Body
Height_max = 10.769 * 10 ^ -2; % [m] Max Height
Width_max = 9.955 * 10 ^ -2; % [m] Max Width
Length_body = 62.611 * 10 ^ -2; % [m] Length of fuselage
Height_base = 0.172 * 62.611 * 10 ^ -2; % [m] Height of fuselage
Width_base = 0.058 * 2 * 62.611 * 10^ -2; % [m] Width of fuselage
f = (Length_body) / sqrt((4/pi) * Height_max * Width_max); % Length divided by diameter of body
cf_body = 0.455 / ( (log10(Re_L)) ^ (2.58) * (1 + 0.144 * M ^ 2) ^ (0.65)); % Skin friction of body
FF_body = (0.9 + (5 / (f ^ (1.5))) + f / 400); % Free form of body
Q_body = 1; % Interference drag
S_wet_body = 2* Height_max * Width_max + 2 * Width_max * Length_body + 2 * Length_body * Height_max; % Wetted area of body

% Wing
x_c = 0.3; % x/c for big wing
t_c = 0.08; % t/c for big wing
swing_theta = 35; % Swing angle for big wing
cf_wing = lam_skin_fric * x_cr / c_mean + turb_skin_fric * (c_mean - x_cr) / c_mean; % Skin friction for big wing
FF_wing = (1 + (0.6 / (x_c) ) * (t_c) + 100 * (t_c)^4) * (1.35 * M ^ (0.18) * (cos(deg2rad(swing_theta))) ^ (0.28)); % Free form of big wing
Q_wing = 1.1; % Interference drag for big wing
S_wet_wing = 0.1864; % Wetted area of big wing
 
% Tails
cr_tail = 10.955 * 10^-2; % Root chord for tail wing
ct_tail = 3.81 * 10^-2; % Tip chord for tail wing
c_mean_tail = (ct_tail + cr_tail) / 2; % Mean chord of tail wing
Re_l = c_mean_tail*rho*V/mu; % Reynold's number of tail wing
turb_skin_fric_tail = 0.074 / (Re_l ^ 0.2); % Turbulent skin friction of tail wing
lam_skin_fric_tail = 1.328 / sqrt(Re_l); % Laminar skin friction of tail wing
x_c_tail = 0.3; % x/c for tail wing
t_c_tail = 0.12; % t/c for tail wing
swing_theta_tail = 45; % Swing angle for tail wing
cf_tail = lam_skin_fric_tail * (x_cr / c_mean_tail) + turb_skin_fric_tail * ((c_mean_tail - x_cr) / c_mean_tail); % Skin friction for tail wing
FF_tail = (1 + (0.6 / (x_c_tail) ) * (t_c_tail) + 100 * (t_c_tail)^4) * (1.36 * M ^ (0.18) * (cos(deg2rad(swing_theta_tail))) ^ (0.28)); % Free form for tail wing
Q_tail = 1.1; % Interference drag for tail wing (randomly chosen)
S_wet_tail = 0.0115; % Wetted area of tail wing

% Extras
wing = cf_wing*FF_wing*Q_wing*S_wet_wing; % Total wing drag
body = cf_body*FF_body*Q_body*S_wet_body; % Total body drag
tails = 2*cf_tail*FF_tail*Q_tail*S_wet_tail; % Total tail drag
CD_misc = ((0.139 + 0.419 * (M - 0.161) ^ 2) * Height_base * Width_base) * (1 / (S)); % Miscelaneous drag coefficient
C = ( wing+body+tails) / (S) + CD_misc;
CD_lp = C*0.02; % Add leakage drag coefficient
CD0 = C + CD_lp; % Parasite drag

%% Drag model validation and testing
% Raymer's method for calculating oswald factor
e1 = 4.61 * (1 - 0.045 * AR ^ 0.68) * (cos(deg2rad(swing_theta))).^(0.15) - 3.1;

% % Obert's method for calculating oswald factor
e2 = 1 / (1.05 + 0.007 * pi * AR);

% Kroo's method for calculating oswald factor
u = 0.99; % Equal to Oswald's theoretical factor
d_F = (Height_max * Width_max); % Diameter of Fuselage
s = 1 - 2 * (d_F / b) ^ 2; % Lift dependent drag
K = 0.38; % Determined experimentally
e3 = 1 / ((1/ u * s) + K * CD0 * pi * AR);

% Stinton Method
e_inv = 0.95; % Inviscid coefficient
e4 = 1 / ((1/e_inv) + 0.45);
% import given shuttle data
Shuttle_Data = readtable('ShuttleTruth.xlsx');

% polyfit to derive function relating CD to CL
poly_Data = polyfit(Shuttle_Data.CL, Shuttle_Data.CD,5); %5th order polynomial
y_Fit = polyval(poly_Data,Shuttle_Data.CL);

% list of CL values to test
CL_span = 0:0.1:1;

% Coefficient of drag as a function of coefficient of lift
CD = @(CL,e) CD0 + (CL.^2 / (pi * e * AR));

% list of drag coefficients derived from the above function
CD_span1 = CD(CL_span,e1);
CD_span2 = CD(CL_span,e2);
CD_span3 = CD(CL_span,e3);
CD_span4 = CD(CL_span,e4);

% Plotting all relevant data
figure(1)
hold on;
plot(Shuttle_Data.CL,Shuttle_Data.CD, 'b--.', 'MarkerSize',20,'LineWidth',2)
plot(Shuttle_Data.CL, y_Fit, 'LineWidth',2,'Color','r');
plot(CL_span, CD_span1, 'b--.', 'MarkerSize',20,'LineWidth',2, 'Color', [0,255/255,0])
plot(CL_span, CD_span2, 'b--.', 'MarkerSize',20,'LineWidth',2, 'Color', [255/255,100/255,0])
plot(CL_span, CD_span3, 'b--.', 'MarkerSize',20,'LineWidth',2, 'Color', [255/255,255/255,0])
plot(CL_span, CD_span4, 'b--.', 'MarkerSize',20,'LineWidth',2, 'Color', [0,0.5,1])
xlabel('C_L','FontSize',15);
ylabel('C_D','FontSize',15);
title('C_D vs. C_L For the Whole Shuttle','FontSize',20);
legend('Given Data', 'Line of Best Fit', "Raymer's method", ...
    "Obert's method", "Kroo's method", "Stinton Method" , ...
    'FontSize',20,'Location','northwest');
hold off;

%% Glide Flight modelling
% parameters
ground_elevation = 5200; %[ft] ground height at which glider is tested
launch_height = 10; %[m] height above ground at which glider is launched
W = 1; % [lb] Aircraft weight
e = e4; % Sets oswald efficiency factor

% retrieves standard atmosphere conditions
[atm.temp, ~, atm.pres, atm.rho] = atmosisa((ground_elevation * 0.3048));

% coefficient of lift and coefficient of drag calculations
CL = (CD0 * pi * e * AR) ^ (1/2); % CL and CD assume that the aircraft
CD = 2 * CD0;                  % flies at max L/D
disp("CL/CD max: " + (CL/CD))

% calculate glide range
glide_range = launch_height * (CL / CD); %[m] range across the ground of the glider

% calculate glide angle
theta = (pi / 2) - atan(CL/CD);
theta = rad2deg(theta);

v_glide = sqrt((2*W)/(atm.rho*S*CL)); % glide velocity [m/s]
sink_rate = v_glide*sind(theta); % sink rate [m/s]

%% Boost Team Code


 


%% List initial conditions

const = getConst(); % call function to access structure of constants

x0 = const.x0; %initial launch distance (horizontal) [m]
y0 = const.y0; % intial wind drift position [m]
z0 = const.z0; % inintial launch height [m]
vx0 = const.vx0; % initial launch velocity in the x-direction(North) [m/s]
vy0 = const.vy0; % initial wind velocity in the y-direction(East) [m/s]
vz0 = const.vz0; % initial launch velocity in the z-direction(Down) [m/s]
V0_air = const.V_air; % initial air volume [m^3]
m0_air = const.m_air; % initial mass of air [kg]
m0_rocket = const.m_rocket; % initial mass of rocket with air and water [kg]

% Set initial condition vector
X0 = [x0 y0 z0 vx0 vy0 vz0 V0_air m0_air m0_rocket];

% Set integration time span
tspan = [0 10];

% Set integration parameters
odepar = odeset('RelTol',(1e-7),'AbsTol',1e-9);

% Call ode45 with function bottle_rocketEOM
[t,X] = ode45(@(t,X) bottle_rocketEOM(t,X),tspan,X0,odepar);

% Pre-allocate vectors for additional variable collection
thrustCorrection = zeros(height(t),1);
Thrust = zeros(height(t),1);
V_AIR = zeros(height(t),1);

% Determine apogee and collect index and apogee height
[h_max, index] = max(X(:,3));

% Construct array of x, y, and z coordinates of apogee for glide function
array = X(index,1:3);

% Loop to extract thrust, air volume, and change in mass data from function bottle_rocketEOM
for i = 1:height(t)
    [DX, Thrust(i,1), V_AIR(i,1)] = bottle_rocketEOM(t(i),X(i,:));
    thrustCorrection(i,1) = (const.m_rocket*const.g0)-(X(i,7)*const.g0);
end

% Calculate glide phase
G = glide(array,sink_rate,const.wind(2),theta);

% Plotting boost and glide trajectory
figure(2)
hold on
grid on
plot3(X(1,1),X(1,2),X(1,3),'^',"MarkerEdgeColor","r","MarkerFaceColor","k")
plot3(X(2:index,1),X(2:index,2),X(2:index,3),"LineWidth",2)
plot3(G(1,1),G(1,2),G(1,3),'>',"MarkerEdgeColor","g","MarkerFaceColor","k")
plot3(G(2:end-1,1),G(2:end-1,2),G(2:end-1,3),'g-.',"LineWidth",2)
plot3(G(end,1),G(end,2),G(end,3),'v',"MarkerEdgeColor","c","MarkerFaceColor","k")
title('Boost-Glide Trajectory')
legend('Launch','Boost Phase','Vehicle Transition','Glide Phase','Landing')
xlabel('Horizontal Distance [m]')
ylabel('Wind-Drift Distance [m]')
zlabel('Height [m]')
view(-176,20)
hold off


cg_rocket = zeros(height(X),1);

 for i = 1:height(X)

    cg_rocket(i) = centerGrav(const.V_bottle - X(i,7));

 end 

plot(t(1:300,:), cg_rocket(1:300,:))
title('Boost-Glide Center of Gravity Evolution')
xlabel('Time [s]')
ylabel('Distance From Vehicle Front [m]')

bob = "Final Center of Gravity, (measured from the front of the vehicle) = ";
cg_glide = cg_rocket(end); % [m] measured from the front of the vehicle.
disp(bob + cg_glide + " meters")


dole = "Total Vehicle length = ";
vehicle_length = const.vehicle_length; % [m]
disp(dole + vehicle_length + " meters")



%% Subfunction for atmospheric constants

%% Subfunction getConst which sets up the constants
% Date Created: 04/15/2022
% Date modified: 02/01/2023
% Inputs:  N/A
% Outputs: const = structure containing constant parameters to use with ODE45 integration
function const = getConst()
% Certain constant value were obtained from the following website
% http://www.aircommandrockets.com/labs/pressure_chambers/pc_1250_001/pc_1250_001.htm

[temp, ~, pres, rho] = atmosisa((5200 * 0.3048)); % inputs: location elevation in ft converted to m

% environmental constants
const.g0 = 9.81; % gravity vector [x,z] [m/s^2]
const.P_atm = pres; % atmospheric pressure [Pa]
const.T_air = temp; % atmospheric temperature [K]
const.rho_air = rho; % atmospheric density [kg/m^3]
const.R = 287; % air constant [J/kg*K]
const.gamma = 1.4; % air specific heat ratio
const.wind = [0 5 0]; % wind velocity vector [m/s]
const.uK = 0.6; % kinetic friction coefficient

% water constants
const.V_water = 0.000625; % initial volume of water [m^3]
const.rho_water = 1000; % water density [kg/m^3]
const.m_water = const.V_water*const.rho_water; % initial mass of water [kg]
% bottle parameters

const.m_bottle = 0.039; % mass of bottle and fins [kg]
const.V_bottle = 0.00125; % bottle volume [m^3]
const.P_bottle = 40*6895; % bottle pressure [Pa]
const.P_gage = const.P_atm+const.P_bottle; % initial gauge pressure [Pa]
const.V_air = const.V_bottle-const.V_water; % initial air volume [m^3]
const.r_bottle = 0.0425; % bottle radius [m]
const.r_throat = 0.0105; % bottle throat radius [m]
const.Cd = 0.75; % discharge coefficient
const.At = pi*(const.r_throat^2); % bottle throat area [m^2]
const.Ac = pi*(const.r_bottle^2); % bottle cross sectional area [m^2]
const.m_air = (const.P_gage*const.V_air)/(const.R*const.T_air); % initial mass of air [kg]
const.m_rocket = const.m_water+const.m_bottle+const.m_air; % initial mass of rocket, air, and water [kg]
const.l_bottle = const.V_bottle / (pi * const.r_bottle^2);  % [m]
const.height_water = const.V_water / (pi * const.r_bottle^2); % [m]

%% Aerospike
const.m_spike = 0; % [kg] 
const.l_spike = 0; % [m]
const.x_cg_spike = const.l_spike / 2; % [m]

%% Front Cone 
const.m_cone_front = 70 * 0.001; % [kg]
const.cg_cone_front = 150 / 1000; % [m]
const.x_cg_cone_front = const.l_spike + const.cg_cone_front; % [m]

%% Rear Cone
const.m_cone_rear = 20 * 0.001; % [kg]
const.cg_cone_rear = 20.95 / 1000; % [m]
const.x_cg_cone_rear = const.l_spike + const.l_bottle + const.cg_cone_rear; % [m]

%% Wings
const.m_wing = 22.85 * 0.001; % [kg]
const.l_wing = 523 / 1000; % [m]
const.cg_wing = (const.l_wing * cosd(30)) / 2; % [m]
const.x_cg_wing = const.l_spike + const.cg_wing; % [m] 

%% Horizontal Stabilizer
const.m_stab = 15.45 * 0.001; % [kg]
const.cg_stab = 20 / 1000; % [m]
const.x_cg_stab = const.l_spike + const.l_bottle + const.cg_stab; % [m]

%% Vertical Stabilizer
const.m_tail = 4.08 * 0.001; % [kg]
const.cg_tail = 52.32 / 1000; % [m]
const.x_cg_tail = const.l_spike + const.l_bottle + const.cg_tail; % [m]

const.vehicle_length = const.l_spike + const.l_bottle + .04;

% flight parameters
const.CD = 0.0648; % drag coefficient
const.z0 = 0.25; % initial vertical height [m]
const.vz0 = 0; % initial vertical velocity [m/s]
const.y0 = 0; % initial orthogonal distance [m]
const.vy0 = 0; % initial orthogonal wind velocity [m/s]
const.x0 = 0; % initial horizontal distance [m]
const.vx0 = 0; % initial horizontal velocity [m/s]
const.ls = 0.5; % launch stand length [m]
const.theta0 = 60; % initial launch angle [deg]
const.h0 = [cosd(const.theta0),0,sind(const.theta0)];
end


%% Boost Function
% Date Created: 04/15/2022
% Date Modified: 02/01/2023
% Inputs:   t = anonymous time variable
%           X = anonymous state vector containing [x;z;vx;vz;]         
% Outputs: dX = derivative state vector containing [dx;dz;dvx;dvz] = [vx;vz;ax;az]
% Methodology: function used for integrating the flight of a bottle rocket
% using equations of motion and mass flow rates

function [dX, f_thrust, V_air] = bottle_rocketEOM(t,X)
% Call function to pass the structure of constants into this EOM function
const = getConst();
% Extract state vector
x = X(1); % From the variable X coming in, extract the index corresponding to horizontal position
y = X(2); % From the variable X coming in, extract the index corresponding to orthogonal position
z = X(3); % From the variable X coming in, extract the index corresponding to vertical position
vx = X(4); % From the variable X coming in, extract the index corresponding to velocity in the x-direction
vy = X(5); % From the variable X coming in, extract the index corresponding to velocity in the y-direction
vz = X(6); % From the variable X coming in, extract the index corresponding to velocity in the z-direction
V_air = X(7); % From the variable X coming in, extract the index corresponding to volume of air
m_air = X(8); % From the variable X coming in, extract the index corresponding to mass of air
m_rocket = X(9); % From the variable X coming in, extract the index corresponding to mass of rocket
% declare relationships that are independent of phase
V = [vx vy vz]-const.wind; % Subtract the wind velocity from the rocket velocity
Vmag = norm(V); % magnitude of the rockets total velocity
q = 0.5*const.rho_air*Vmag^2; % dynamic pressure
fDrag = q*const.CD*const.Ac; % drag force
fGrav = -const.g0*m_rocket; % force of gravity

% phase 1: powered flight with water
if V_air <= const.V_bottle % check that there is still water in the rocket by using the conditional statement 
    % that the volume of air is less than the total volume of the bottle
    p_air = const.P_gage*((const.V_air/V_air)^(const.gamma)); % air pressure inside of bottle
    f_thrust = 2*const.Cd*const.At*(p_air-const.P_atm); % thrust for water powered flight phase
    dV_air = const.Cd*const.At*sqrt((2*(p_air-const.P_atm))/const.rho_water); % change in air volume
    dm_rocket = -const.Cd*const.At*sqrt((2*const.rho_water)*(p_air-const.P_atm)); % change in mass of rocket
    dm_air = 0; % change in mass of air
else
    p_end = const.P_gage*(const.V_air/const.V_bottle)^(const.gamma);
    p_air = p_end*(m_air/const.m_air)^(const.gamma);
end

% phase 2: powered flight phase with air after the water has been exhausted
if V_air >= const.V_bottle && p_air > const.P_atm % check that there is no water in the rocket by using the conditional statement that 
    % the air pressure inside the bottle is greater than atospheric pressure
    p_end = const.P_gage*((const.V_air/const.V_bottle)^(const.gamma)); % pressure at the end of the bottle
    p_air = p_end*((m_air/const.m_air)^(const.gamma)); % air pressure inside of bottle
    p_critical = p_air*(2/(const.gamma+1))^(const.gamma/(const.gamma-1)); % critical pressure
    rho_air = m_air/const.V_bottle; % air density inside bottle
    T = p_air/(rho_air*const.R); % air temperature inside bottle
        if p_critical > const.P_atm % Conditiona statement checking if critical pressure is greater than atmospheric pressure
            p_exit = p_critical; % exit pressure
            T_exit = T*(2/(const.gamma+1)); % exit temperature
            V_exit = sqrt(const.gamma*const.R*T_exit); % exit velocity
            rho_exit = p_exit/(const.R*T_exit); % exit density
        else
        p_exit = const.P_atm; % exit pressure
        M_exit = sqrt(((((p_air/const.P_atm)^((const.gamma-1)/const.gamma))-1))/((const.gamma-1)/2)); % exit mach number
        T_exit = T/(1+(((const.gamma-1)/2)*(M_exit^2))); % exit temperature
        rho_exit = const.P_atm/(const.R*T_exit); % exit density
        V_exit = M_exit*sqrt(const.gamma*const.R*T_exit); % exit velocity
        end
    
    dm_air = -const.Cd*rho_exit*const.At*V_exit; % change in mass of air
    dV_air = 0; % change in volume of air
    dm_rocket = dm_air; % change in mass of rocket
    f_thrust = (-dm_air*V_exit)+(const.At*(const.P_atm-p_exit)); % force of thrust delivered by air pressure flight phase
end

% phase 3: ballistic phase after water and air are exhausted
if V_air >= const.V_bottle && p_air <= const.P_atm % check that there is no water or air left in the rocket 
    % all thrust and changes in mass and volume are zero in the ballistic
    % phase
    f_thrust = 0;
    dm_rocket = 0;
    dV_air = 0;
    dm_air = 0;
  
end

if z <= 0 % checking to see if the rocket has hit the ground
vx = 0;
vy = 0;
vz = 0;
end

% conditional statements that apply during all phases
if z < ((1/2)*sind(const.theta0))+0.25 % checking to see if the rocket has left the test stand
xh = const.h0(1);
yh = const.h0(2);
zh = const.h0(3);
% Implement guide rail friction
f_friction = const.uK*m_rocket*const.g0*xh;
else % heading for all other phases of flight
f_friction = 0;
xh = V(1)/Vmag;
yh = V(2)/Vmag;
zh = V(3)/Vmag;
end

% expressions assigning variables into output vector
dx = vx; % horizontal position 
dy = vy; % orthogonal position
dz = vz; % vertical position
dvx = (f_thrust-fDrag-f_friction)*(xh/m_rocket); % horizontal veloity
dvy = (f_thrust-fDrag-f_friction)*(yh/m_rocket); % wind drift velocity
dvz = (f_thrust-fDrag-f_friction)*(zh/m_rocket)+(fGrav/m_rocket); % vertical velocity

% Write the output vector dX containing the equations to be integrated
dX = [dx;dy;dz;dvx;dvy;dvz;dV_air;dm_air;dm_rocket];
end


%% Glide Function
% Date Created: 02/22/2023
% Inputs: X - 1x3 array containing initial x, y, and z positions [m]
%         Vs - Sink rate of the glider [m/s]
%         Vw - Wind velocity value [m/s]
%         theta - glide decsent angle [degrees]
% Outputs: Xmat - a matrix of x, y, and z positions as a function of the
% linearly spaced glide time vector

function [Xmat] = glide(X,Vs,Vw,theta)
x0 = X(1); % x position (horizontal) at apogee [m]
y0 = X(2); % y position (orthogonal) at apogee [m]
z0 = X(3); % z position (vertical) at apogee [m]

t = z0/Vs; % total glide time [s]

t = linspace(0,t); % linear spaced time vector for additional data points

V = Vs/sind(theta); % velocity magnitude [m/s]

vx = V*cosd(theta); % velocity in x-direction (constant) [m/s]

vy = Vw; % velocity in y-direction [m/s]

vz = Vs; % velocity in z-direction (constant) [m/s]

Xmat = zeros(length(t),3); % Pre-allocating vector for xyz outputs

for i = 1:length(t)

Xmat(i,1) = x0 + vx*t(i); % x position (horizontal) as a function of time [m]
Xmat(i,2) = y0 - vy*t(i); % y position (orthogonal) as a function of time [m]
Xmat(i,3) = z0 - vz*t(i); % z position (vertical) as a function of time [m]

end

end

%% Changing Center of Gravity Function
% Date Created: 03/08/2023
% NOTE: x locations measured from the FRONT of the vehicle
% Inputs: Volume of water
% Outputs: 

function [cg_int] = centerGrav(v_water)

const = getConst();
      
%Calculating Center of Gravity of the Water
x_water = const.l_spike + (const.l_bottle - (v_water / (2 * pi * (const.r_bottle^2))));

%Calculating the Numerator of CG
num = ((const.m_spike * const.x_cg_spike) + ...
    (const.m_cone_front * const.x_cg_cone_front) + ...
    (const.m_cone_rear * const.x_cg_cone_rear) + ...
    (2 * const.m_wing * const.x_cg_wing) + ...
    (2 * const.m_stab * const.x_cg_stab) + ...
    (2 * const.m_tail * const.x_cg_tail) + ...
    ((v_water * const.rho_water) * x_water));
        
%Calculating the Denominator of CG
denom = (const.m_spike + const.m_cone_front + const.m_cone_rear + ...
    (2 * const.m_wing) + (2 * const.m_tail) + (2 * const.m_stab) + ...
    (v_water * const.rho_water));

%Final CG calculation
cg_int = num / denom;



end
