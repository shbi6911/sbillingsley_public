%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%define random arrays for values subject to error
r = (75-0.5)+((75+0.5)-(75-0.5))*randn(1,500);
l = (260-0.5)+((260+0.5)-(260-0.5))*randn(1,500);
d = (155-0.5)+((155+0.5)-(155-0.5))*randn(1,500);

%pull experimental data
tnum = 5; % number of the test being analyzed
str = "Locomotive_Data_2020\Test1_" + (tnum+4) + "pt5V";
[theta_exp, w_exp, V_exp, time] = LCSDATA(str);
theta = theta_exp; % [deg]
w = w_exp; % [deg/s]

%preallocate
residual_monte = zeros(252,length(r));

%run Monte Carlo
for i=1:length(r)
    [V_slide] = LCSMODEL(r(i), d(i), l(i), theta_exp, w_exp);
    residual = (V_exp - V_slide);
    residual_monte(:,i)=residual;
end

%calculate mean and std dev
residual_monte_mean = mean(residual_monte);
residual_monte_mean = sort(residual_monte_mean);
residual_monte_stddev = std(residual_monte);
residual_monte_stddev = sort(residual_monte_stddev);

%plotting
figure(1);
plot(residual_monte_mean);
title("Residual Mean of Monte Carlo");
ylabel("Mean of Residuals (mm)");
xlabel("Monte Carlo Trial");

figure(2);
plot(residual_monte_stddev);
title("Standard Deviation of Monte Carlo");
ylabel("Std Deviation of Residuals (mm)");
xlabel("Monte Carlo Trial");

