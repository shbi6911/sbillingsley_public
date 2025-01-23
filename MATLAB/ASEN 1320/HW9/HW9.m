%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

%This script numerically integrates the flight of a projectile in the x and
%y dimension.  For a set of given constants and initial velocity and launch
%angle, it uses a function EOM to calculate rates of change and uses ode45
%to integrate this function.  It then plots the output using color for a
%third dimension, using the color_line3d function.

%%          %given constants
mass = 65.0;
radius = 0.38;
Cdrag = 0.5;
rho = 1.275;
g = 9.81;
parameters = [Cdrag,radius,mass,rho,g];
%%         %initial values
tspan = [0,10];
V0 = 50;
theta = 45;
state0 = [V0*cosd(theta);V0*sind(theta);0;100];

 %%         %numerical integration
EOMFun = @(time,state)EOM(time, state, parameters);
[timeVector, stateMatrix] = ode45(EOMfun,tspan,state0);

%%          %drag calculation
Drag = zeros(length(timeVector),1);
for ii = 1:length(timeVector)
    Drag(ii,1) = 0.5*rho*Cdrag*(pi*radius*radius)*(stateMatrix(ii,1)^2+stateMatrix(ii,2)^2);
end

%%          %plotting data
%creates a plot with four subplots, all with x,y data plotted.  Color is
%used to indicate Time, Velocity in X direction, Drag force, and Velocity
%in the Y direction, respectively, clockwise starting from upper left.

ColorPlotFun = @(color)color_line3d(color, stateMatrix(:,3), stateMatrix(:,4));
figure (1);
subplot(2,2,1),plot (stateMatrix(:,3), stateMatrix(:,4));
subPlot1 = ColorPlotFun(timeVector);
subplot (2,2,1),title('Time (s)'),xlabel('x (meters)'),ylabel('y (meters)');
subplot (2,2,2),plot (stateMatrix(:,3), stateMatrix(:,4));
subPlot2 = ColorPlotFun(stateMatrix(:,1));
subplot (2,2,2),title('X velocity (m/s)'),xlabel('x (meters)'),ylabel('y (meters)');
subplot (2,2,3),plot (stateMatrix(:,3), stateMatrix(:,4));
subPlot3 = ColorPlotFun(stateMatrix(:,2));
subplot (2,2,3),title('Y velocity (m/s)'),xlabel('x (meters)'),ylabel('y (meters)');
subplot (2,2,4),plot (stateMatrix(:,3), stateMatrix(:,4));
subPlot4 = ColorPlotFun(Drag);
subplot (2,2,4),title('Drag Force (N)'),xlabel('x (meters)'),ylabel('y (meters)');