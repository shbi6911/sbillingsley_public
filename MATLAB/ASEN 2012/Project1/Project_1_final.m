%% ASEN 2012 Project One Calorimetry 
%By Shane Billingsley (110231742) & Zach Sollins (110675868)
% date created: [10/3/22]
% last modified: [10/18/22]

% OUTLINE: This code reads in given calorimeter data (time and
% temperature.)  It uses this data, along with given values of mass and
% specific heat, to find the specific heat of an unknown sample material.
% This specific heat can be compared to a list of possiblities to identify
% the sample.
%INPUTS: A whitespace-delimited data file of time and temperature data 
% from a calorimeter
%OUTPUTS: Specific heat with error value of unknown sample, outputted to
%the display window, and four figures showing data and best fit lines.

%% read in and clean data
data = readmatrix("SampleB");
time = data(:,2);
thermo1 = data(:,3); % calorimeter 1 thermocouple
thermo2 = data(:,4); % sample thermocouple (boiling water) 
thermo3 = data(:,5); % room temperature thermocouple
thermo4 = data(:,6); % calorimeter 2 thermocouple

%% average calorimeter thermocouples and determine error
thermo_cal = mean([thermo1 thermo4],2); %cal1 and cal2
%combine the error of the mean in quadrature
thermo_cal_sigma_x = sqrt(((thermo1 - thermo_cal)+(thermo4 - thermo_cal).^2));
thermo_cal_err = thermo_cal_sigma_x./sqrt(2); %standard error of the mean


%find values for endpoints of our fit lines by observing plots
%determine indices of these values in temperature data
%disp(find(time == 442.461)); %endpoint of line 1
%disp(find(time == 681.44));  %start point of line 2
%disp(find(time == 500.739)); %end point of line 3
%disp(find(time == 439.625)); %t0

% split arrays to calculate three best fit lines
% note that line numbers do not go strictly left to right across graph
% the line numbers also do not match what was used in our report
time_line1 = time(1:312);             %line 1 is before t0
thermo_cal_line1 = thermo_cal(1:312);

time_line2 = time(480:746);         %line 2 is after calorimeter equalizes
thermo_cal_line2 = thermo_cal(480:746);

time_line3 = time(315:350);        %line 3 is during equalization
thermo_cal_line3 = thermo_cal(315:350);


%% find least-squares fit lines using split arrays
line1_fit = leastSquaresLine(time_line1,thermo_cal_line1);
line2_fit = leastSquaresLine(time_line2,thermo_cal_line2);
line3_fit = leastSquaresLine(time_line3,thermo_cal_line3);
%anonymous function to calculate fit lines
line = @(fit,x) fit(2)*x + fit(1);

%% find values for equation temperatures
t0 = 439;  %time when sample is added
T_L = line(line1_fit,439);  %intersection of line 1 and t0
T_H = line(line2_fit,439);  %intersection of line 2 and t0
T_avg = (T_L+T_H)/2;         %avg temperature
t_avg = (1/line3_fit(2)*T_avg)-(line3_fit(1)/line3_fit(2)); %time @ T_avg

%% define all values for calorimetry equation
T0 = T_L;            %temp of calorimeter at t0
T1 = thermo2(310);    %temp of sample at t0
T2 = line(line2_fit,t_avg);  %temp of calorimeter at equilibrium
m_c = 0.51*1000; % mass of calorimeter in g
c_c = 0.895; %specific heat of calorimeter in J/(g*degC)
m_s = 39.306; %mass of sample in g

%% calculate specific heat of sample and output
c_s = (m_c*c_c*(T2-T0))/(m_s*(T1-T2));
disp("C_s");
disp(c_s);

%% compute partial derivatives
m_c_deriv = (c_c*(T2-T0))/(m_s*(T1-T2));
m_s_deriv = -(m_c*c_c*(T2-T0))/((m_s^2)*(T1-T2));
T0_deriv = -(m_c*c_c)/(m_s*(T1-T2));
T1_deriv = -(m_c*c_c*(T2-T0))/(m_s*(T1-T2)^2);
T2_deriv = (m_c*c_c*(T1-T0))/(m_s*(T1-T2)^2);

%% define error values of given masses
m_c_err = 0.05; %error value of calorimter mass in grams
m_s_err = 0.001; %error value of sample mass in grams

%% calculate sigma_y, W, and P for each fit line
[sigma_y1,W1] = sigma_y(time_line1,thermo_cal_line1,line1_fit);
[sigma_y2,W2] = sigma_y(time_line2,thermo_cal_line2,line2_fit);
[sigma_y3,W3] = sigma_y(time_line3,thermo_cal_line3,line3_fit);
P1 = covariance(time_line1,W1);
P2 = covariance(time_line2,W2);
P3 = covariance(time_line3,W3);

%% calculate sigma y of each line using alternate method
%line 1
AAA= line(line1_fit,time_line1); % arbitray variable for math
sigY1 = sum((thermo_cal_line1-AAA).^2); % sum of difference
err_siggY1 = sqrt(sigY1/(312-2)); % (error) sigma y
%line 2
BBB= line(line2_fit,time_line2);
sigY2 = sum((thermo_cal_line2-BBB).^2);
err_siggY2 = sqrt(sigY1/(267-2));
%line 3
CCC= line(line3_fit,time_line3);
sigY3 = sum((thermo_cal_line3-CCC).^2);
err_siggY3 = sqrt(sigY3/(length(thermo_cal_line3)-2));


%% create error bounds on each line of fit
line_one_err = line(line1_fit,0:444); %matrix size work around
line_two_err = line(line2_fit,419:1070); %matrix size work around
line_three_err = line(line3_fit,440:515); %matrix size work around
L1errplus = line_one_err + sigY1; % error up
L1errmin = line_one_err - sigY1;  % error down
L2errplus = line_two_err + sigY2; % err up
L2errmin = line_two_err - sigY2;  % err d
line_three_errp = line(line3_fit,440:503);
L3errplus = line_three_errp + sigY3; % err up 
L3errmin = line_three_err - sigY3;  % err d


%% find sigma_q for each line
time_line1_extrap = (0:450)';
sigma_q1 = sigma_q(time_line1_extrap,P1);
time_line2_extrap = (400:1200)';
sigma_q2 = sigma_q(time_line2_extrap,P2);
sigma_q3 = sigma_q(time_line3,P3);

%%  find error of T_L and T_H
% should use sigma_y or sigma_q, whichever is bigger
if (2*sigma_y1) >= (2*sigma_q1(t0)) %compare for T_L on line 1
    T_L_err = sigma_y1;
else
    T_L_err = sigma_q1(t0);
end

if (2*sigma_y2) >= (2*sigma_q2(39)) %compare for T_H on line 2
    T_H_err = sigma_y2;
else
    T_H_err = sigma_q2(39);
end

%%  find error of T0 and T2
T0_err = T_L_err;  %T_L is used for T0, so their error is equal
T_avg_err = norm([T_L_err T_H_err]);%combine T_L and T_H error in quadrature
% use the General Method applied to line 3 to determine error in t_avg
% which is the time at which we presume equalization was complete
t_avg_err = norm([((1/line3_fit(2))*T_avg_err),... 
   ((-1/line3_fit(2))*sqrt(P3(1,1))),...
   (((T_avg-line3_fit(1))/line3_fit(2)^2)*sqrt(P3(2,2)))]);
%use the General Method to combine error in t_avg with errors of line 2
%coefficients to arrive at an error value for T2
T2_err = norm([(line2_fit(2)*t_avg_err),(t_avg*sqrt(P2(2,2))),(sqrt(P2(1,1)))]);

%% find error of T1
%fit a line to the values of thermocouple 2 before sample was removed
%and placed in calorimeter
thermo2fit = leastSquaresLine(time(1:310),thermo2(1:310));
%calculate sigma_y of this line
[sigma_y_thermo2,~] = sigma_y(time(1:310),thermo2(1:310),thermo2fit);
%use a 2 sigma value for error in T1
T1_err = 2*sigma_y_thermo2;

%% calculate final error value and output
%use error values and previously calculated derivatives in the General
%Method
c_s_err = norm([(m_c_err*m_c_deriv),(T2_err*T2_deriv),(T0_err*T0_deriv),...
    (m_s_err*m_s_deriv),(T1_err*T1_deriv)]);
disp("error in C_s")
disp(c_s_err);

%% Plots

%plot with fit lines-------------------------------------------------------
figure (1); hold on; 
plot(0:444,line(line1_fit,0:444),'red','LineWidth',.2);
plot(443:508,line(line3_fit,443:508),'red','LineWidth',.2);
plot(419:1070,line(line2_fit,419:1070),'red','LineWidth',.2);
plot (time,thermo_cal,"g");
legend("Fit Line 1","Fit Line 2","Fit Line 3","Sample B");
title("Sample B with lines of fit");
xlabel("Time (S)");
ylabel("Tempurature (C)");

%plot with fit & error & marked temps-------------------------------------
%line one-----------------------------------------------------------------
figure (2); hold on; 
plot(0:444,line(line1_fit,0:444),'m','LineWidth',.2);
% line one error bars 
plot(L1errplus,'r:','LineWidth',.1);
plot(L1errmin,'r:','LineWidth',.1);
plot (time,thermo_cal,"g");
%markers for T values
plot(439,25.203,'or',"markersize",5); % marker for T0
% indicating lines
xline(439);  % sample added
legend("Fit Line","sigma y","sigma y","Sample B","T0"); %this will be annoying
title("Line One data");
xlabel("Time (S)");
ylabel("Tempurature (C)");



%line three---------------------------------------------------------------
figure(3); hold on; 
plot(443:508,line(line3_fit,443:508),'m','LineWidth',.2);
% line three error bars
plot(440:503,L3errplus,'r:','LineWidth',.1);
plot(440:515,L3errmin,'r:','LineWidth',.1);
plot (time,thermo_cal,"g"); % plotting data last so its shown on top of other lines
xline(t_avg); % avg time
legend("Fit Line","sigma y","sigma y");
title("Line two data");
xlabel("Time (S)");
ylabel("Tempurature (C)");


%line two-----------------------------------------------------------------
figure(4); hold on; 
plot(419:1070,line(line2_fit,419:1070),'m','LineWidth',.2);
% line two error bars
plot(419:1070,L2errplus,'r:','LineWidth',.1);
plot(419:1070,L2errmin,'r:','LineWidth',.1);
plot (time,thermo_cal,"g");
% indicating lines
xline(439);  % sample added
xline(t_avg); % avg time
%markers for T values
plot(439,29.77,'ok',"markersize",5); % marker for TH
plot(t_avg,T2,'or',"markersize",5); % marker for T2
legend("Fit Line","sigma y","sigma y","Sample B","combined","midpoint","TH","T2",'Location','southeast');
title("Line Three data");
xlabel("Time (S)");
ylabel("Tempurature (C)");









%% functions

%This function fits a least squares line to a given set of x and y data
%INPUTS: column vectors of x and y data
%OUTPUTS: A column vector of line coefficients in the order [b m]
function fit = leastSquaresLine(x,y)
N = length(x);
H = [ones(N,1) x];
fit = inv(H'*H)*(H'*y);
end

%This function calculates the covariance matrix for a given line fit
%INPUTS: a column vector of x data and a W matrix
%OUTPUTS: a 2x2 covariance matrix of error values
function P = covariance(x,W)
N = length(x);
H = [ones(N,1) x];
P = inv(H'*W*H);
end

%This function calculates sigma_y for a given line fit
%INPUTS: column vectors of x and y data, and a fit matrix with coefficients
%OUTPUTS: a scalar value for sigma_y, and a weighting matrix with 
% 1/sigma_y^2 on the diagonal
function [sigma_y, W] = sigma_y(x,y,fit)
sigma_y = sqrt(1/(length(x)-2) * sum(((fit(2)*x+fit(1))-y).^2));
W = (1/sigma_y^2)*eye(length(x));
end

%This function calculates sigma_q for a given line fit, on a given set of
%x-data, which may or may not be the original x-data
%INPUTS: a column vector of x-data and a 2x2 covariance matrix
%OUTPUTS: a column vector of sigma_q values for each x in the input vector
function sigma_q = sigma_q(x,P)
N = length(x);
H = [ones(N,1) x];
sigma_q = diag(sqrt(H * P * H'));
end