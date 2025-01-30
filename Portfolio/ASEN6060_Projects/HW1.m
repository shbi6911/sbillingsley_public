% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025

%% Problem 1
% nondimensionalize the Earth-Moon and Sun-Earth systems

clear; clc; const = getConst(0);
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

clear; clc; 
const = getConst(1);

%code for manual input of initial states and integration times
% prompt = {'Enter initial state vector as space-separated numbers',...
%     'Enter integration time span as space-separated numbers'};
% dlgtitle = 'CR3BP Integration';
% fieldsize = [1 60; 1 25];
% definput = {'0 0 0 0 0 0','0 1'};
% answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
% initial_state = str2num(answer{1});
% tspan = str2num(answer{2});

%preset initial states and time spans
initial_states = [0.98 0 0 0 1.2 0;0.98 0 0 0 1.7 0;0.12 0 0 0 3.45 0;...
    0.12 0 0 0 3.48 0];
tspans = [0 2;0 8;0 25;0 25];
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);

for i = 1:length(tspans)
    [results.("t_" + string(i)),results.("state_matrix_" + string(i))] = ...
        ode89(EOMfun, tspans(i,:), initial_states(i,:), opts);
    results.("jc_" + string(i)) = ...
        jacobi(results.("state_matrix_" + string(i)),const.mu.em);
end

if const.vrb
    grid on;
    x = results.state_matrix_1(:,1);    y = results.state_matrix_1(:,2);
    p1 = plot(x,y);
    axis equal
    frac = round(length(x)/2);
    p2 = arrow([x(1),y(1)],[x(1+1),y(1+1)]);        %plot dir of motion
    arrow([x(frac),y(frac)],[x(1+frac),y(1+frac)]);
    p3 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
    xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    title("CR3BP Trajectory 1");
    legend([p1,p2,p3],"Trajectory","Direction of Travel","Lunar Diameter",...
        'Location','southoutside');
    %saveas(gcf,"ASEN6060_HW1_Fig1",'png');
    
    figure(); grid on;
    x = results.state_matrix_2(:,1);    y = results.state_matrix_2(:,2);
    p1 = plot(x,y);
    axis equal
    frac = round(length(x)/2);
    p2 = arrow([x(1),y(1)],[x(1+1),y(1+1)]);
    arrow([x(frac),y(frac)],[x(1+frac),y(1+frac)]);
    p3 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
    xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    title("CR3BP Trajectory 2");
    legend([p1,p2,p3],"Trajectory","Direction of Travel","Lunar Diameter",...
        'Location','southoutside');
    %saveas(gcf,"ASEN6060_HW1_Fig2",'png');
    
    figure(); grid on;
    x = results.state_matrix_3(:,1);    y = results.state_matrix_3(:,2);
    p1 = plot(x,y);
    axis equal
    frac = round(length(x)/3);
    frac = round(length(x)/2);
    p2 = arrow([x(1),y(1)],[x(1+1),y(1+1)]);
    arrow([x(frac),y(frac)],[x(1+frac),y(1+frac)]);
    p3 = circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r','');
    p4 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
    xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    title("CR3BP Trajectory 3");
    legend([p1,p2,p3,p4],"Trajectory","Direction of Travel","Earth Diameter",...
        "Lunar Diameter",'Location','southoutside');
    %saveas(gcf,"ASEN6060_HW1_Fig3",'png');
    
    figure(); grid on;
    x = results.state_matrix_4(:,1);    y = results.state_matrix_4(:,2);
    p1 = plot(x,y);
    axis equal
    frac = round(length(x)/8);
    p2 = arrow([x(1),y(1)],[x(1+1),y(1+1)]);
    arrow([x(7*frac),y(7*frac)],[x(1+7*frac),y(1+7*frac)]);
    p3 = circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r','');
    p4 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
    xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    title("CR3BP Trajectory 4");
    legend([p1,p2,p3,p4],"Trajectory","Direction of Travel","Earth Diameter",...
        "Lunar Diameter",'Location','southoutside');
    %saveas(gcf,"ASEN6060_HW1_Fig4",'png');
end

%% Problem 2d
clear; clc;
const = getConst(1);
initial_states = [0.98 0 0 0 1.2 0;0.98 0 0 0 1.7 0;0.12 0 0 0 3.45 0;...
    0.12 0 0 0 3.48 0];
tspans = [0 2;0 8;0 25;0 25];
%dimensionalize initial state 3
init = initial_states(3,:);
init_dim = zeros(1,length(init));
tim = tspans(3,:);
init_dim(1:3) = init(1:3).*const.l_star.em;
init_dim(4:6) = init(4:6).*(const.l_star.em/const.t_star.em);
tim_dim = tim(2)*const.t_star.em;
dist = norm(init_dim(1:3)) + (const.mu.em*const.l_star.em);
tim_dim_period = tim(2)/(2*pi);
tim_period_check = tim_dim/(2*pi*const.t_star.em);

%% Problem 3
% write an event function with specified parameters to stop integration of
% an ode solver, run a specified trajectory, and plot it

%clear; clc; 
const = getConst(1);
initial_state = [0.12 0 0 0 3.45 0];    %trajectory 3
tspan = [0 25];

%calculate trajectory using specified event function
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14,'Events',@eventFun);
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
[t,state_matrix] = ode89(EOMfun, tspan, initial_state, opts);
final_state = state_matrix(end,:);
final_time = t(end);

if const.vrb
    %plot trajectory
    grid on;
    x = state_matrix(:,1);    y = state_matrix(:,2);
    p1 = plot(x,y);
    axis equal
    frac = round(length(x)/3);
    p2 = arrow([x(frac),y(frac)],[x(frac+1),y(frac+1)]);
    p3 = circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r','');
    %p4 = circle((1-mu_em),0,const.radius.Moon/l_star_em,'k');
    xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
    title("CR3BP Trajectory 3 with Stopping Condition");
    legend([p1,p2,p3],"Trajectory","Direction of Travel","Earth Diameter",...
        'Location','southoutside');
end

%% Problem 4
% write code to find zero-velocity curves, plot some specified curves, and
% use zero-velocity curves to estimate Jacobi constant values for the
% Lagrange points of the Earth-Moon system

clear; clc; 
const = getConst(1);
x = linspace(-1.5,1.5,1000);     %define nondimensional coordinates for Earth-Moon system
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);         
%calculate a matrix of Jacobi constants at zero velocity
r1 = sqrt((X + const.mu.em).^2 + Y.^2);
r2 = sqrt((X - 1 + const.mu.em).^2 + Y.^2);
Z = (X.^2 + Y.^2) + (2*(1-const.mu.em))./r1 + (2*const.mu.em)./r2;

%set specified values of Jacobi constant
jconst = [3.189,3.173,3.013,2.995];

if const.vrb
    %plot zero-velocity curves for specified values
    for i = 1:length(jconst)
        figure();
        colormap(gray)
        contourf(X,Y,-(Z - jconst(i)),[0 0]);
        hold on;
        p1 =scatter(nan,nan,[],[0.5,0.5,0.5],'filled',DisplayName="Forbidden Areas");
        p2 =scatter(nan,nan,[],[1,1,1],'filled',DisplayName="Allowable Areas");
        p3 = circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r',"Earth Diameter");
        p4 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k',"Moon Diameter");
        xlabel("X-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
        ylabel("Y-Position ($l^*$ = 1 nondimensional unit)",'interpreter','latex');
        title("Earth-Moon CR3BP Zero-Velocity Curve for Jacobi Constant " + string(jconst(i)));
        legend([p1,p2,p3,p4],'Location','southoutside');
        axis equal
    end
end

%% Problem 4d

clear; clc; 
const = getConst(0);

x = linspace(-1.5,1.5,1000);     %define nondimensional coordinates for Earth-Moon system
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);         
%calculate a matrix of Jacobi constants at zero velocity
r1 = sqrt((X + const.mu.em).^2 + Y.^2);
r2 = sqrt((X - 1 + const.mu.em).^2 + Y.^2);
Z = (X.^2 + Y.^2) + (2*(1-const.mu.em))./r1 + (2*const.mu.em)./r2;

if const.vrb
    figure();
    %colormap('gray');
    contourf(X,Y,Z,linspace(2.98,3.00,100));
    circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'r','');
    axis equal
end

%% Functions
function drdt = threeBP_EOM(t,r,mu)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025

% threeBP_EOM calculates the equations of motion for the third body within
% a Circular Restricted Three Body Problem, using nondimensional quantities
% for time and state vector, and a rotating coordinate frame.  It is
% intended for use with ODE solvers such as ODE45 or ODE89.
%
% INPUTS:   t       nondimensional time unit (from ODE solver)
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
%eventFun is an event function for use with threeBP_EOM and an ODE solver.
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

function h = circle(x,y,r,color,name)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025
%
% circle is a simple function to plot a circle of defined center and radius
%
%INPUTS:    x,y     coordinates of the circle center
%           r       radius of the circle
%           color   color value of the plotted circular line
%           name    string or char array of the DisplayName of the line obj
%
%OUTPUTS    h       a Line object with properties as above
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
    h.Color = color;
    h.DisplayName = name;
    hold off
end