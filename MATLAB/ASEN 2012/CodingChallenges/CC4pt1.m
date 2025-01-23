%% [Shane Billingsley] - SID: [110231742]
% ASEN 2012: Coding Challenge 4
% Last modified: [10/17/22]


%% Housekeeping

clc         % Clear command window
close all   % Close out of any open figures

%% Define important variables for Euler Integration

% Anonymous function of differential equation
g = @(x) exp(-x^2);

% All the parameters required to solve Euler's method
x_init = 0;
h = 0.1;
n_steps = 20;
x = x_init:h:(n_steps*h);
y=ones(1,n_steps);

%% Solve the Integration

% Hint: you will want to use a for-loop for this one, since it cannot be
% done using vector/matrix calculations. This is because one step 
% will depend on the previous one, which usually prohibits vectorisation
for i = 1:n_steps
    y(i+1) = y(i) + h*g(x(i));
end


%% Plot Data

% Load in accurate data
load("AccurateData.mat");
% Plot both solutions
figure;
hold on;
plot(x_accurate,y_accurate);
plot(x,y,'red');


% Format the figure. Pay attention to the order of the legend
title("Numerical Integration using Euler's Method");
legend('Accurate', 'Euler', 'Location', 'nw');
xlabel('x');
ylabel('y');


