% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025

%% Problem 1
% nondimensionalize the Earth-Moon and Sun-Earth systems

clear; clc; const = getConst();
% find nondimensional mass quantities
m_star_em = (const.mu.Earth/const.G) + (const.mu.Moon/const.G);
m_star_se = (const.mu.Sun/const.G) + (const.mu.Earth/const.G);
% find mass ratios
mu_em = (const.mu.Moon/const.G)/m_star_em;
mu_se = (const.mu.Earth/const.G)/m_star_se;
% find nondimensional length quantities
l_star_em = const.a.Moon;
l_star_se = const.a.Earth;
% find nondimensional time quantities
t_star_em = sqrt((l_star_em^3)/(const.G*m_star_em));
t_star_se = sqrt((l_star_se^3)/(const.G*m_star_se));

%% Problem 2
%integrate and plot four specified trajectories in the Earth-Moon circular
%restricted 3 body problem

%clear; clc; const = getConst();
prompt = {'Enter initial state vector as space-separated numbers',...
    'Enter integration time span as space-separated numbers'};
dlgtitle = 'CR3BP Integration';
fieldsize = [1 60; 1 25];
definput = {'0 0 0 0 0 0','0 1'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
initial_state = str2num(answer{1});
tspan = str2num(answer{2});

opts = odeset('RelTol',10e-12,'AbsTol',10e-12);
EOMfun = @(t,r)threeBP_EOM(t,r,mu_em);
[t,state_matrix] = ode89(EOMfun, tspan, initial_state, opts);

grid on;
plot(state_matrix(:,1),state_matrix(:,2));
%earth = circle(-mu_em,0,const.radius.Earth/l_star_em);
moon = circle((1-mu_em),0,const.radius.Moon/l_star_em);
axis equal

%% Problem 3
%clear; clc; const = getConst();
initial_state = [0.12 0 0 0 3.45 0];
tspan = [0 25];

opts = odeset('RelTol',10e-12,'AbsTol',10e-12,'Events',@eventFun);
EOMfun = @(t,r)threeBP_EOM(t,r,mu_em);
[t,state_matrix] = ode45(EOMfun, tspan, initial_state, opts);

grid on;
plot(state_matrix(:,1),state_matrix(:,2));
earth = circle(-mu_em,0,const.radius.Earth/l_star_em);
moon = circle((1-mu_em),0,const.radius.Moon/l_star_em);
axis equal

%% Problem 4
clear; clc; const = getConst();

%% Problem 5
clear; clc; const = getConst();

%% Functions
function drdt = threeBP_EOM(t,r,mu)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025

% threeBP_EOM calculates the equations of motion for the third body within
% a Circular Restricted Three Body Problem, using nondimensional quantities
% for time and state vector, and a rotating coordinate frame.  It is
% intended for use with ODE45.
%
% INPUTS:   t       nondimensional time unit (from ODE45)
%           r       nondimensional state vector, as a column vector, of the
%                   form [x;y;z;x_dot;y_dot;z_dot]
%           mu      mass ratio parameter
%
% OUTPUTS:  drdt    nondim derivative state vector, as a column vector, 
%                   of the form [x_dot;y_dot;z_dot;x_ddot;y_ddot;z_ddot]

    drdt = zeros(6,1);      %preallocate output vector
    drdt(1:3) = r(4:6);     % derivative of position is velocity
    %calculate respective radii
    r1 = sqrt((r(1) + mu)^2 + r(2)^2 + r(3)^2);
    r2 = sqrt((r(1) - 1 + mu)^2 + r(2)^2 + r(3)^2);
    %calculate accelerations per CR3BP equations of motion
    %note x,y,z = r(1,2,3)  x_dot,y_dot,z_dot = r(4,5,6)
    drdt(4) = 2*r(5) + r(1) - ((1-mu)*(r(1)+mu))/r1^3 - (mu*(r(1)-1+mu))/r2^3;
    drdt(5) = -2*r(4) + r(2) - ((1-mu)*r(2))/r1^3 - (mu*r(2))/r2^3;
    drdt(6) = -((1-mu)*r(3))/r1^3 - (mu*r(3))/r2^3;
end

function [value, isterminal, direction] = eventFun(t,y)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025
%
%eventFun is an event function for use with threeBP_EOM and an ode solver.
% It stops integration at the first point where y = 0 and the function y is
% increasing
%
%INPUTS:    t,y     time and state from the ode solver
%
%OUTPUTS    value       the state at which y=0
%           isterminal  a value of 1 indicates stop integration
%           direction   a value of 1 indicates the function is increasing
    value = y(2);
    isterminal = 1;
    direction = 1;
end

function h = circle(x,y,r)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025
%
% circle is a simple function to plot a circle of defined center and radius
%
%INPUTS:    x,y     coordinates of the circle center
%           r       radius of the circle
%
%OUTPUTS    h       a Line object with properties as above
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
    hold off
end