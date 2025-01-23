%% ASEN 2012: Coding Challenge #4: Runge-Kutta and Euler Integration 
% by: [Shane Billingsley]
% SID: [110231742]
% last modified: [10/18/22]
%
% OUTLINE: 
% The purpose of this coding challenge is to familiarize yourself with some
% techniques for numerical integration, using both a first-order euler
% integration as well as a 4th-order Runge-Kutta expansion.
%
% To accomplish this, you are tasked with creating two functions to perform
% these methods, and apply them with variables step sizes to test their
% effectiveness.
% 

% housekeeping
clc; close all

% we wish to integrate with the following step sizes: pi, pi/2, pi/4, pi/8
% the function in question is given by the differential equation:
%     g(t,y) = y*(sin(t))^2
%
% with the following initial condition:
%     y(0) = pi
%
% which has the explicit solution below:
%     y(t) = C*e^(t/2 - sin(2t/4))
%
% using integration techniques, adjust the value for "C" until your curve
% looks appropriate

step = [pi,pi/2,pi/4,pi/8]; % step sizes to iterate through
g = @(t,y) y.*(sin(t)).^2; % anonymous function handle for differential equation 
IC = [0,pi]; % initial conditions [t0,y0]
tf = 3*pi; % set final time

% true solution is defined here, you'll need to try different values of C
C = pi;
y = @(t) C*exp(t/2 - sin(2*t)/4); % true solution

% now, run your functions in a FOR loop to get numerical expressions that
% we can fit the curve y(t) to

f = figure();
f.Position = [800,65,650,850];
for i = 1:length(step)
    n = tf/step(i); % set number of integration intervals
    [t_e,y_e] = euler(g,step(i),n,IC);% call euler function
    [t_rk,y_rk] = rk4(g,step(i),n,IC);% call rk function
    
    % plot the results in a m x 1 subplot, where m is the number of step
    % sizes you used 
    subplot(length(step),1,i); hold on
        plot(t_e,y_e,'*-') % plot euler outputs
        plot(t_rk,y_rk,'o--') % plot rk4 ouputs
        plot(t_e,y(t_e),':','LineWidth',1.5) % plot true solution
        
        title("Comparison of Integration Techniques: \Deltat = \pi/" + string(pi/step(i)))
        xlabel('Time (t)')
        ylabel('y-axis (y)')
        legend('Euler Approximation','Runge-Kutta Approximation','True Solution','location','northwest')
        grid on; grid minor
    hold off
end


%% Function definitions

% euler integration function
function [t,y] = euler(g,dt,n,IC)
    
    % pre-allocate vectors for t and y and set initial conditions
    t = zeros(1,n); t(1) = IC(1);
    y = zeros(1,n); y(1) = IC(2); 

    % run euler integration in a FOR loop
    for i = 1:n
        t(i+1) = dt*i;
        y(i+1) = y(i)+dt*g(t(i),y(i)); 
    end
    
end

% runge-kutta integration function
function [t,y] = rk4(g,dt,n,IC)

    % pre-allocate vectors for t and y and set initial conditions
    t = zeros(1,n); t(1) = IC(1);
    y = zeros(1,n); y(1) = IC(2); 

    % run runge-kutta integration in a FOR loop
    for i = 1:n
        t(i+1) = dt*i;
        
        % expressions for intermediate slopes (k1-k4)
        k1 = g(t(i),y(i));
        k2 = g((t(i)+(dt/2)),(y(i)+(k1*(dt/2))));
        k3 = g((t(i)+(dt/2)),(y(i)+(k2*(dt/2))));
        k4 = g((t(i)+dt),(y(i)+(k3*dt)));
        
        % combine to forward-integrate using RK expansion
        y(i+1) = y(i) + (1/6)*(dt)*(k1 + 2*k2 + 2*k3 +k4);        
    end

end