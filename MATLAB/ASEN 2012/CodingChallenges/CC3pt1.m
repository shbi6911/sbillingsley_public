%% ASEN 2012: Coding Challenge #3 - Linear Least Squares Regression
% by: Shane Billingsley
% SID: 110231742
% last modified: 10/5/22
%
% OUTLINE: 
% When fitting a line of best fit, we utilize a method called the Least
% Squares Estimation: minimizing the squared error between a set of data
% and some line of "best fit." For this coding challenge, you'll be
% provided a set of data and be asked to perform a Least Squares regression
% by using the matrix method (since this is MATLAB, after all)
% 


% clear the command window and close all open figures
clc; close all

%% Problem 1: Compute Linear Least Squares Regression

% read in and extract data (mind your units: 1 Volt = 1 Ohm * 1 Ampere)
data = readmatrix('IV_Characteristic.xlsx');
V = data(:,2); % voltage values from 0V to 3.4V [V]
I = data(:,1); % corresponding measured current values [mA]
    I = I/1000; % convert current to [A]
N = length(data); % compute the length of the data stream

% for our LLS fit, we need to solve "d = A*x_hat"
% start with defining the solution matrix, d (depedant variable)
d = V;

% derive an expression for the LS approximation matrix: A
% recall that the voltage for this circuit can be expressed as 
%               V = R*I + V0
%         or    y = m*x + b

Acol1 = I; % column 1 represents "m" or the x^1 term
Acol2 = ones(N,1); % column 2 represents "b" or the x^0 term
A = [Acol1 Acol2]; % concatenate columns into one matrix

% compute the coefficient matrix: d = A*x_hat
x_hat = inv(A'*A)*(A'*d);

% compute the best estimate for the resistance and reference voltage
R_best = x_hat(1);
V0_best = x_hat(2);

% write an anonymous function handle
voltageLeastSquares = @(m,b,x) m*x+b;

% implement step 9 of the 10 step method by checking the solution with an 
% alternate approach, such as solving the system with the \ operator
x_hat = (A'*A)\(A'*d);

% plot the input data as points, as well as the LS regression
figure(1); hold on;
scatter (I,V,10,'blue','filled');
xlabel ("Voltage (V)");
ylabel ("Measured Current (A)");
title ("Finding Resistance of Unknown Resistor")

plot(I,voltageLeastSquares(R_best,V0_best,I),'red');






hold off

%% Problem 2: Applying the Covariance Matrix 
% given that the error in the current measurement has a SYSTEMATIC ERROR of
% +/- 10μA, find an expression for the covariance matrix and associated LS 
% error


% compute the RESIDUAL SUM OF SQUARES of the line of best fit
RSS = sum((V-voltageLeastSquares(R_best,V0_best,I)).^2);
err_RSS = sqrt((1/(N-2))*RSS);

% write an expression for the SYSTEMATIC ERROR
err_sys = 1e-5; % [A]

% combine these errors in quadrature
sigmaI = norm([err_sys err_RSS]);

% Create W matrix (hint: type >> help diag or help eye in command line)
W = (1/(sigmaI)^2)*eye(N);

% Solve for covariance matrix: P = inv(A'*W*A)
P = inv(A'*W*A);

% extract errors from the covariance matrix 
sigmaR = sqrt(P(1,1));
sigmaV0 = sqrt(P(2,2));

% add errorbars with the total uncertainty to the figure
totalErr = diag(sqrt(A*P*A'));

figure(1); hold on
    % use the errorbar() function
errorbar(I,voltageLeastSquares(R_best,V0_best,I),totalErr,'red')    
    
hold off

