%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%% Housekeeping
tic
clc
clear
close all

% file setup
outfilename = "table.txt";
outfile = fopen(outfilename, "w");
fprintf(outfile, "Experiment number: Mean: Absolute Value Mean: Standard Deviation: Absolute Value Std Dev \n");

% for loop performs calculations on all experiments
for i = 1:6
    %% Import data
    tnum = i; % number of the test being analyzed
    str = "Locomotive_Data_2020\Test1_" + (tnum+4) + "pt5V";
    [theta_exp, w_exp, V_exp, time] = LCSDATA(str);
    
    
    %% Slide velocity
    % component dimensions
    r = (7.5 * 10); % [mm]
    d = (15.5 * 10); % [mm]
    l = (26 * 10); % [mm]
    theta = theta_exp; % [deg]
    w = w_exp; % [deg/s]
    
    % run caluculations
    [V_slide] = LCSMODEL(r, d, l, theta_exp, w_exp);
    
    % plotting
    fig = figure (tnum);
    hold on
    plot(theta_exp, V_exp/10)  %convert to cm and plot
    plot(theta, V_slide/10)
    titlestr = "Angular position vs. slide velocity (Experiment " + tnum + "):";
    title(titlestr)
    xlabel("Degrees")
    ylabel("cm/s")
    legend("Experimental Data", "Model Data")
    hold off

    %% residuals
    % calculations
    residual = (V_exp - V_slide);
    mean_exp = mean(residual);
    mean_exp_abs=mean(abs(residual));
    deviation = std(residual);
    deviation_exp_abs=std(abs(residual));

    % write to table:
    fprintf(outfile, (tnum + " " + mean_exp + " " + mean_exp_abs + " " + ...
        deviation + " " + deviation_exp_abs + "\n"));
    hold off

    % plot residuals
    fig2 = figure(tnum + 6);
    hold on
    plot(time, residual/10)  %convert to cm and plot
    title("Residual vs. time (Experiment " + tnum + "):")
    xlabel("Time [s]")
    ylabel("Residual [cm]")
    hold off
    
end
fclose(outfile);

toc

%{
Calculates the vertical speed of the sleeve

Inputs :
    r {float} - radius of the rotating disk [mm]
    d {float} - distance between the center of the disk and the wall [mm]
    l {float} - length of rod connecting disk to collar [mm]
    theta {float} - angle of the disk relative to the vertical [deg]
    w {float} - angular velocity of the disk (derivative of theta) [deg/s]

Outputs:
    V_slide {float} - vertical speed of the collar [mm/s]
%}
function [V_slide] = LCSMODEL(r, d, l, theta, w)

% convert theta and w to radians and rads/s repsectively
theta = deg2rad(theta);
w = deg2rad(w);
                                     
beta = asin((d - (r .* sin(theta))) ./ l); % [ rad] angle between rod and collar

V_slide = -(w .* r .* sin(theta)) - (w .* r .* cos(theta) .* tan(beta));
end

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