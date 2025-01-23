%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%{
Retrieves the appropriate variables from a given file, and returns 

file {str] - file name for the data to be retrieved

dat [struct{double}} - data
%}
function [theta, w, v, t] = LCSDATA(file)
% import data
rawData = importdata(file);

% assign variables
t = rawData.data(:,1);
theta = rawData.data(:,2);
% position_slide = rawData.data(:,3);
w = rawData.data(:,4);
v = rawData.data(:,5);
% sample_time = rawData.data(:,6);

% fix degrees to start between 0 and 360
f = floor(theta(1) / 360);
theta = theta - (360.*f);

end