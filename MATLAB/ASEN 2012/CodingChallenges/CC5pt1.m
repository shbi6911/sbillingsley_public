%% ASEN 2012: Coding Challenge #5: Integrating a 1st Order ODE with ODE45 
% by: Shane Billingsley
% SID: 110231742
% last modified: 10/31/22
%
% OUTLINE: 
% The purpose of this coding challenge is to practice using ode45 on a simplified first 
% order system. You will create an expression to model the acceleration of a weather balloon,
% eventually finding his terminal velocity.

% housekeeping
clc
close all

%% Main Script

% givens
r = 17; %radius of the balloon
g0 = 9.81; % gravitational acceleration [m/s^2]
CD = 0.5; % coefficient of drag
rho_air = 1.225; % density of air (outside) [kg/m^3]
rho_gas = 0.1786; % density of helium (inside) [kg/m^3]

% intermediate calculations
V = (4/3)*pi*r^3;% volume [m3]
Ac = pi*r^2; % cross-sectional area [m2]

% masses
m_payload = 470; % mass of balloon + payload (not including mass of the gas) [kg]
m_gas = V*rho_gas; % recall mass = density*volume [kg] 
m_sys = m_payload+m_gas; % total mass of the system

% set up ode45 integration
tspan = [0 10]; % in [s]
IC = 0;% Intial conditions (starting velocity) in [m/s]
const = struct('radius',r,'gravity',g0,'drag',CD,'air_density',rho_air,...
    'gas_density',rho_gas,'volume',V,'area',Ac,'mass',m_sys); % Vector or structure containing important geometric and constant values. Note the index (if using a vector) ...
        % or the name (if using a structure) so you can call these values within your subfunction

% plotting
figure(); hold on
    % call ode45
    [t,v] = ode45(@(t,v) calcAccel(t,v,const),tspan,IC); % Using the syntax required for ODE45 call the function here. 
    plot(t,v)
hold off

% You will see that the balloon reaches terminal velocity. Write an expression to determine what the terminal velocity is.
vterm = v(end);


%% End of Main Function

%% Start of Subfunction

function accel = calcAccel(t,v,const)
    % calcAccel: a function used for integration of a first order system,
    % taking into account the forces impacting a balloon in steady flight
    %
    % INPUTS: t,v,const
    %  t is a scalar value of time, v is a scalar velocity
    %  const is a structure of constant values
    % OUTPUTS: fnet, accel
    % fnet is the net force on the balloon, accel is its current
    % acceleration
    % Note: your inputs and outputs of your subfunction must be the same variables named below
    
    % compute individual forces
    fGrav =  const.mass*const.gravity;
    fBuoyancy = const.volume*const.air_density*const.gravity;
    fDrag = 0.5*const.area*const.drag*const.air_density*(v^2);
    
    
    % expression for net force and acceleration
    fnet = fBuoyancy-fDrag-fGrav;
    accel = fnet/const.mass;
    %disp("fnet "); disp(fnet);
    %disp("accel "); disp(accel);
 
end
