%% ASEN 2012: Coding Challenge #7
% NAME: [Shane Billingsley]
% SID: [110231742]
% Last modified: [11/28/22]
%

% The purpose of this coding challenge it to familiarize yourself with root
% finding methods for use in your later aerospace curriculum. You will be
% using methods discussing in class as well as built-in MATALB functions to
% complete these tasks.

% housekeeping
clc; close all

%% Define error tolerance
s = 4; % number of decimals we want precision to
epsilon = 10^-(s+1); % tolerance (same for all methods)

%% Define constants
alpha = 0.008;
beta = 0.0002;
gamma = -0.26;

% set up inputs for our rootfunc function
f = @(x) ((gamma*x)./(1-(beta*x)))-(alpha*x.^2); % original function of energy density rho
df = @(x) (gamma/(1-(beta*x)).^2)-(2*alpha*x); % Derivative of function with respect to rho
g1 = @(x) ((alpha*x.^2)*(1-(beta*x)))/gamma;
g2 = @(x) (gamma/(alpha*(1-(beta*x))));
%g3 = @(x) (sqrt((gamma*x)/(alpha-(alpha*beta*x))));

% determine a reasonable starting point for your recursive formulas
guess = [-32 1]; % Initial guess(es) [hint: you should have as many distinct initial guesses you have roots - how can you approximate the values of these roots?

%% Call function to calculate roots once per guess
roots = zeros(5,length(guess));

for i = 1:length(guess)
     % each column of "roots" represents all five solutions for one guess
    switch i
        case 1
            roots(:,i) = rootsfunc(guess(i),f,df,g1,epsilon);
        case 2
            roots(:,i) = rootsfunc(guess(i),f,df,g2,epsilon);
        case 3
            roots(:,i) = rootsfunc(guess(i),f,df,g3,epsilon);
    end
    % HINT: the input function you use for g will depend on which root
    % you're solving for, try using a switch command to set these
    %   i.e.   "switch i"
    
end

% printing results: don't change anything here
method = categorical({'Bisection Method';'Successive Substitution';'Newtons Method';'Secant Method';'MATLAB Built-ins'});
tab = table(method,roots);






function z = rootsfunc(x0,f,df,g,epsilon) % Define function to calculate roots for all methods
% ROOTSFUNC(x0,f,df,g,epsilon) evaluates the convergence of some function,
% f, using a variety of methods as outlined for ASEN 2012.
% 
% INPUTS:
%   x0      : initial guess
%   f       : original function
%   df      : derivative function
%   g       : back-solved function
%   epsilon : error tolerances
%
% OUPUTS:
%   z       : a 5 x 1 vector containing the roots for each method
%

% initialize output
z = zeros(5,1);
fprintf("For initial guess: %.3g\n",x0)

%% Bisection method - the general process for this can be seen in slide 15 of lecture 20
% Calculate how long this method takes to complete (do this for every single method you use)
tic

% hint: for the bisection method, see the function sign()

% describe initial params
a = x0 - 10;
b = x0 + 10;
c = x0;

% set initial guesses s.t. a and b have different signs
while sign(f(a)) == sign(f(b))
    a = a-1;
    b = b+1;
end

% modify and iterate c until f(c) is within tolerance
while abs(f(c)) >= epsilon
    if sign(f(a)) ~= sign(f(c))
        b = c;
    else 
        a = c;
    end
    c = (a+b)/2;
end

% output into z vector
z(1) = c;
fprintf("Bisection method: %0.6f s\n",toc)

%% Successive substitution - the general process for this can be seen in slide 6 of lecture 21
tic

% manipulate original f(x) function to get x = g(x) [x, in this case, is rho (independent variable)]

% set up inputs
x1 = x0;
x2 = g(x1);
%iterate until |x2 - x1| is within some tolerance
while abs(x2-x1) >= epsilon
    x1 = x2;
    x2 = g(x1);
end

% output into z vector
z(2) = x2;
fprintf("Successive Substitution: %0.6f s\n",toc)

%% Newton's method - the general process for this can be seen on slide 9 of lecture 21
tic

% set up inputs
x1 = x0;
x2 = x0 - (f(x0)/df(x0));

% iterate until |x2 - x1| is within some tolerance
while abs(x2-x1) >= epsilon
    x1 = x2;
    x2 = x1 - (f(x1)/df(x1));
end

% output into z vector
z(3) = x2;
fprintf("Newton's method: %0.6f s\n",toc)

%% Secant method - the general process for this can be seen on slide 14 of lecture 21
tic

% set up inputs
xn = x0-1;
x1 = x0;
x2 = x1 - ((f(x1)*(x1-xn))/(f(x1)-f(xn)));
% iterate until |x3 - x2| is within some tolerance
while abs(x2-x1) >= epsilon
    xn = x1;
    x1 = x2;
    x2 = x1 - ((f(x1)*(x1-xn))/(f(x1)-f(xn)));
end

% output into z vector
z(4) = x2;
fprintf("Secant method: %0.6f s\n",toc)

%% Check using MATLAB built-in functions: fzero(), roots(), or fsolve()
% Don't forget "doc" and "help" are useful tools to find syntax of built-in functions
tic

z(5) = fzero(f,x0);

fprintf("MATLAB built-ins: %0.6f s\n",toc)
fprintf("\n\n")

end