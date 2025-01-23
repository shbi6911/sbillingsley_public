% Contributors: Shane Billingsley and Gabriel Law
% Course number: ASEN 3801
% File name: 
% Created: 1/17/24

clc;
clear;
close all;

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

%% Functions
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