%% ASEN 2012: Coding Challenge #2
% by: [Shane Billingsley]
% SID: [110231742]
% last modified: [9/16/22]
%

%% Problem 1
% You should develop one function to calculate the standard error of the
% mean, or SEM, and one function to calculate the weighted average and its
% associated uncertainty. Both of these functions should be called from the
% main script. If you wish to, you are allowed to use built-in MATLAB
% functions for intermediate calculations like the mean and standard
% deviation.

%% Functions: hints for syntax and formulation
% All functions should be defined at the very bottom of your main script,
% but it is recommended that you fill in the function before completing the
% code that calls it from the script.

% The number of outputs, number of inputs, and general syntax of any
% function call should correspond directly to the function definition
% itself

% Function syntax for 1 input, 1 output: output = function_name(input)
% Function syntax for multiple inputs, multiple outputs: [output1,output2]
% = function_name(input1,input2)
    % You can also mix and match - number of inputs does not have to equal
    % number of outputs


%% Housekeeping
clc; close all

%% Read In and Clean Data
% thrust data collected by two different measurement devices during a
% static-fire rocket test

% read in data - xlsread, readmatrix, or readtable might be of use here
data = readtable("Static_Thrust_Data.xlsx"); 

t = data.Time_s_; % extract time data in [s]
f = [data.Sensor1Force_N_,data.Sensor2Force_N_]; % extract thrust data (should be a 373x2 matrix) in [N]

% clean both time and force data - think about what values should be 
% included in the weighted average and which values should not

t_clean = t(f(:,1)~=0);
f_clean = f(f(:,1)~=0,:);


% extrapolate size values from post processed matrix          
[n_rows,n_cols] = size(f_clean); % should NOW be a 217x2 matrix

% optional: plot the data to get a sense of what you're working with
figure(); hold on
    plot(t_clean,f_clean)
hold off

% % Compute Standard Error of the Mean
% call SEM function here, from the main part of your code
% 
% compute the SEM for each dataset (once for each column of f)
% hint: a FOR loop will make this easier
% 
% pre-allocate vectors to save space
sigma_f_bar = zeros(n_cols,1);
for i = 1:n_cols
    sigma_f_bar(i) = SEM(f_clean(:,i));
end
%disp(sigma_f_bar)
% 
% % Compute Weighted Average Between The Datasets
% call wAvg function here, from the main part of your code
% 
% we want the weighted average of the mean of the two datasets
% use the standard deviation for each dataset's uncertainty
[f_wav,sigma_wav] = wAvg(mean(f_clean),std(f_clean));
% 
% 
% % Functions
% define your SEM and wAvg functions here

function sigma_x_bar = SEM(x)
   sigma_x_bar = std(x)/sqrt(length(x));
end

% 
function [x_wav,sigma_wav] = wAvg(x,sigma_x) 
%     WAVG(x,sigma_x) computes the weighted average of x with uncertainties
%     defined in the vector, sigma_x, where sigma_x(i) is the uncertainty
%     of x(i)
%     
    if length(x) ~= length(sigma_x)
        warning('Your inputs should be vectors of the same length!')
    end
    
%    compute weighted average
    w = 1./sigma_x.^2;
    disp(w);
    x_wav = (sum(x'.*w)/sum(w)); % sum() and dot() might be useful here
    
%    compute uncertainty in the weighted average
    sigma_wav = 1/sqrt(sum(w));
end