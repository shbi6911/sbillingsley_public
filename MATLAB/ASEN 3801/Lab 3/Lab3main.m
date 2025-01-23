%% Bias and Calibration

clear; close all;   %housekeeping

%read in data files
data_1 = readin("24_02_21_001_MEMS_0.2_0.5");   %first dataset
data_2 = readin("24_02_21_001_MEMS_0.5_0.25");  %second dataset
data_3 = readin("24_02_21_001_MEMS_0.5_0.75");  %third dataset
data_4 = readin("24_02_21_001_MEMS_0.75_0.75");
%calculate calibration values using biasCalc function, omitting fourth
%dataset as an outlier
[K1,b1] = biasCalc("24_02_21_001_MEMS_0.2_0.5");
[K2,b2] = biasCalc("24_02_21_001_MEMS_0.5_0.25");
[K3,b3] = biasCalc("24_02_21_001_MEMS_0.5_0.75");

%find mean and std deviation of bias and scaling factor
mean_scale_factor = mean([K1 K2 K3]);
mean_bias = mean([b1 b2 b3]);
std_scale_factor = std([K1 K2 K3]);
std_bias = std([b1 b2 b3]);

%disp([mean_scale_factor std_scale_factor]);
%disp([mean_bias std_bias]);

%plot first dataset (calibrated)
figure();
hold on; grid on;
plot(data_1(:,1),data_1(:,3));
plot(data_1(:,1),mean_scale_factor.*data_1(:,2)+mean_bias);
title("Time History for 0.2 Hz and 0.5 Amplitude (Calibrated)");
xlabel("Time (s)");
ylabel("Angular Rate (rad/s)");
legend("Encoder Data","Calibrated Gyro Data");
hold off;

%plot second dataset (calibrated)
figure();
hold on; grid on;
plot(data_2(:,1),data_2(:,3));
plot(data_2(:,1),mean_scale_factor.*data_2(:,2)+mean_bias);
title("Time History for 0.5 Hz and 0.25 Amplitude (Calibrated)");
xlabel("Time (s)");
ylabel("Angular Rate (rad/s)");
legend("Encoder Data","Calibrated Gyro Data");
hold off;

%plot third dataset (calibrated)
figure();
hold on; grid on;
plot(data_3(:,1),data_3(:,3));
plot(data_3(:,1),mean_scale_factor.*data_3(:,2)+mean_bias);
title("Time History for 0.5 Hz and 0.75 Amplitude (Calibrated)");
xlabel("Time (s)");
ylabel("Angular Rate (rad/s)");
legend("Encoder Data","Calibrated Gyro Data");
hold off;

%calculate error values (measured - truth)
error_1 = -data_1(:,2) - data_1(:,3);
error_2 = -data_2(:,2) - data_2(:,3);
error_3 = -data_3(:,2) - data_3(:,3);
calibrated_error_1 = (mean_scale_factor.*data_1(:,2)+mean_bias) - data_1(:,3);
calibrated_error_2 = (mean_scale_factor.*data_2(:,2)+mean_bias) - data_2(:,3);
calibrated_error_3 = (mean_scale_factor.*data_3(:,2)+mean_bias) - data_3(:,3);

%plot measured error first dataset
figure();
hold on; grid on;
plot(data_1(:,1),error_1);
title("Time History of Measured Error (O.2 Hz, 0.5 A)");
xlabel("Time (s)");
ylabel("Error (Measured vs. True)");
hold off;

%plot measured error second dataset
figure();
hold on; grid on;
plot(data_2(:,1),error_2);
title("Time History of Measured Error (0.5 Hz, 0.25 A)");
xlabel("Time (s)");
ylabel("Error (Measured vs. True)");
hold off;

%plot measured error third dataset
figure();
hold on; grid on;
plot(data_3(:,1),error_3);
title("Time History of Measured Error (0.5 Hz, 0.75 A)");
xlabel("Time (s)");
ylabel("Error (Measured vs. True)");
hold off;

%plot calibrated error first dataset
figure();
hold on; grid on;
plot(data_1(:,1),calibrated_error_1);
title("Time History of Calibrated Error (0.2 Hz, 0.5 A)");
xlabel("Time (s)");
ylabel("Error (Calibrated vs. True)");
hold off;

%plot calibrated error second dataset
figure();
hold on; grid on;
plot(data_2(:,1),calibrated_error_2);
title("Time History of Calibrated Error (0.5 Hz, 0.25 A)");
xlabel("Time (s)");
ylabel("Error (Calibrated vs. True)");
hold off;

%plot calibrated error third dataset
figure();
hold on; grid on;
plot(data_3(:,1),calibrated_error_3);
title("Time History of Calibrated Error (0.5 Hz, 0.75 A)");
xlabel("Time (s)");
ylabel("Error (Calibrated vs. True)");
hold off;

%plot position error vs. truth data angular rate (first dataset)
figure();
hold on; grid on;
scatter(data_1(:,3),error_1,15,'blue','.');
title("Angular Position Error vs. Encoder Rate (0.2 Hz, 0.5 A)");
xlabel("Encoder Rate Measurement (rad/s)");
ylabel("Angular Position Error (Measured)");
hold off;

%plot position error vs. truth data angular rate (second dataset)
figure();
hold on; grid on;
scatter(data_2(:,3),error_2,15,'blue','.');
title("Angular Position Error vs. Encoder Rate (0.5 Hz, 0.25 A)");
xlabel("Encoder Rate Measurement (rad/s)");
ylabel("Angular Position Error (Measured)");
hold off;

%plot position error vs. truth data angular rate (third dataset)
figure();
hold on; grid on;
scatter(data_3(:,3),error_3,15,'blue','.');
title("Angular Position Error vs. Encoder Rate (0.5 Hz, 0.75 A)");
xlabel("Encoder Rate Measurement (rad/s)");
ylabel("Angular Position Error (Measured)");
hold off;

%print all generated figures to png files

% for i = 1:23
%     print(i,("Lab_3_Figure_"+ string(i)),'-dpng','-r300');
% end

%% Question 3.3

Data10_5 = table2array(readtable("24_02_21_001_RWHEEL_10_5"));
Data10_5 = Data10_5(122:593,:); % cutoff values for when motor was producing torque
Data20_10 = table2array(readtable("24_02_21_001_RWHEEL_20_10"));
Data20_10 = Data20_10(130:392,:);
Data20_5 = table2array(readtable("24_02_21_001_RWHEEL_20_5"));
Data20_5 = Data20_5(124:391,:);
Data5_10 = table2array(readtable("24_02_21_001_RWHEEL_5_10"));
Data5_10 = Data5_10(125:1085,:);
Data5_5 = table2array(readtable("24_02_21_001_RWHEEL_5_5"));
Data5_5 = Data5_5(133:600,:);

Data10_5(:,3) = Data10_5(:,3) * 0.1047;  %rpm to rad/s
Data20_10(:,3) = Data20_10(:,3) * 0.1047;
Data20_5(:,3) = Data20_5(:,3) * 0.1047;
Data5_10(:,3) = Data5_10(:,3) * 0.1047;
Data5_5(:,3) = Data5_5(:,3) * 0.1047;

Data10_5(:,4) = Data10_5(:,4) * 0.0335;  %Amps to Nm
Data20_10(:,4) = Data20_10(:,4) * 0.0335;
Data20_5(:,4) = Data20_5(:,4) * 0.0335;
Data5_10(:,4) = Data5_10(:,4) * 0.0335;
Data5_5(:,4) = Data5_5(:,4) * 0.0335;

%Time10_5 = Data10_5(:,1);
%Time20_10 = Data20_10(:,1);
%Time20_5 = Data20_5(:,1);
%Time5_10 = Data5_10(:,1);
%Time5_5 = Data5_5(:,1);


%Plotting Angular velocity over time
figure(1)
subplot(3,2,1)
plot(Data10_5(:,1), Data10_5(:,3))
p10_5 = polyfit(Data10_5(:,1),Data10_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("10 Nm over 5 seconds")

subplot(3,2,2)
plot(Data20_10(:,1), Data20_10(:,3))
p20_10 = polyfit(Data20_10(:,1),Data20_10(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("20 Nm over 10 seconds")

subplot(3,2,3)
plot(Data20_5(:,1), Data20_5(:,3))
p20_5 = polyfit(Data20_5(:,1),Data20_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("20 Nm over 5 seconds")

subplot(3,2,4)
plot(Data5_10(:,1), Data5_10(:,3))
p5_10= polyfit(Data5_10(:,1),Data5_10(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("5 Nm over 10 seconds")

subplot(3,2,5)
plot(Data5_5(:,1), Data5_5(:,3))
p5_5 = polyfit(Data5_5(:,1),Data5_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("5 Nm over 5 seconds")

%Plotting torque over time , for getting cutoff values
%{
figure(1)
plot(Data10_5(:,1), Data10_5(:,4))
p10_5 = polyfit(Data10_5(:,1),Data10_5(:,3),1);

figure(2)
plot(Data20_10(:,1), Data20_10(:,4))
p20_10 = polyfit(Data20_10(:,1),Data20_10(:,3),1);

figure(3)
plot(Data20_5(:,1), Data20_5(:,4))
p20_5 = polyfit(Data20_5(:,1),Data20_5(:,3),1);

figure(4)
plot(Data5_10(:,1), Data5_10(:,4))
p5_10= polyfit(Data5_10(:,1),Data5_10(:,3),1);

figure(5)
plot(Data5_5(:,1), Data5_5(:,4))
p5_5 = polyfit(Data5_5(:,1),Data5_5(:,3),1);
%}

%MOI = Torque / Angular acceleration
I10_5 = mean(Data10_5(:,4))/p10_5(1);
I20_10 = mean(Data20_10(:,4))/p20_10(1);
I20_5 = mean(Data20_5(:,4))/p20_5(1);
I5_10 = mean(Data5_10(:,4))/p5_10(1);
I5_5 = mean(Data5_5(:,4))/p5_5(1);
I = [I10_5,I20_10,I20_5,I5_10,I5_5];

MeanI = mean(I); %kg m^2
StdDev = std(I);


% Angular Momentum Capacity
%H = I*Omega
H = MeanI*4000*0.1047;  %kg m^2 rad/s

%% Control Data

data = readmatrix("24_03_06_Control_76.47_26.76",'FileType','text',NumHeaderLines=2);
clean_data = data(:,1:4);
clean_data(:,1) = clean_data(:,1) - clean_data(1,1);
clean_data(:,1) = clean_data(:,1)./1000;

figure(); hold on; grid on;
plot(clean_data(:,1),clean_data(:,2));
plot(clean_data(:,1),clean_data(:,3));
legend("Reference Height","Measured Position","Location",'northoutside');
title("Controlled Output with K1 = 76.46, K2 = 26.76");
xlabel("Time in Seconds");
ylabel("Angular Position (radians)");
hold off;

subset = find(clean_data(:,1)>=10 & clean_data(:,1) <= 22);
clean_data_subset = clean_data(subset,:);
figure(); hold on; grid on;
plot(clean_data_subset(:,1),clean_data_subset(:,2));
plot(clean_data_subset(:,1),clean_data_subset(:,3));
legend("Reference Height","Measured Position","Location",'northoutside');
set(gca,'YDir','reverse');
title("Controlled Output with K1 = 76.46, K2 = 26.76");
xlabel("Time in Seconds");
ylabel("Angular Position (radians)");
hold off;






%% functions
function data = readin(filename)
%INPUTS     filename    string value of a filename to read
%
%OUTPUTS    data        cleaned data file read in from filename
%
%   readin reads the values from a file output from LabVIEW set up for the
%   rotating spacecraft lab experiment.  It converts the time array to
%   start at zero seconds and converts the encoder data from RPM to rad/s.
%   It then plots the time history of the measurements.

%read in data from file
data = readmatrix((filename),"FileType",'text',"NumHeaderLines",2);
data(:,1) = data(:,1) - data(1,1);  %start time vector at zero
data(:,3) = data(:,3).*(pi/30);     %convert encoder data from RPM to rad/s

%plot time history
figure();
hold on; grid on;
plot(data(:,1),data(:,2));
plot(data(:,1),data(:,3));
str_breakdown = split(filename,'_');    %parse out filename string for labeling
frequency = str_breakdown(6);
current = str_breakdown(7);
title("Time History for " + frequency + " Hz and " + current + " Amplitude");
xlabel("Time (s)");
ylabel("Angular Rate (rad/s)");
legend("Gyro Data","Encoder Data");
hold off;

%plot time history with gyro rate sign corrected
figure();
hold on; grid on;
plot(data(:,1),-data(:,2));
plot(data(:,1),data(:,3));
str_breakdown = split(filename,'_');
frequency = str_breakdown(6);
current = str_breakdown(7);
title("Time History for " + frequency + " Hz and " + current + " Amplitude");
xlabel("Time (s)");
ylabel("Angular Rate (rad/s)");
legend("Gyro Data","Encoder Data");
hold off;

end

function [K,b] = biasCalc(filename)
%INPUTS         filename    string value of filename (same as readin)
%
%OUTPUTS        K           scalar value of scaling factor
%               b           scalar value of bias
%
%METHODOLOGY    biasCalc reads in data from a file specified by the
%filename input (intended to be a LabVIEW file from the rotating spacecraft
% experiment, cleans the data, generates a linear fit between the gyro data
% and encoder data, uses that linear fit to generate bias and scaling
% factor calibration values for the gyro, and then plots the data and fit
% lines.

%read in data
data = readmatrix((filename),"FileType",'text',"NumHeaderLines",2);
data(:,1) = data(:,1) - data(1,1);  %start time vector at zero
data(:,3) = data(:,3).*(pi/30);     %convert encoder data from RPM to rad/s

p = polyfit(data(:,3),data(:,2),1); %generate linear fit line
K=1/p(1); b=-p(2)/p(1); %assign output values

%plot encoder vs. gyro data
figure();
hold on; grid on;
scatter(data(:,3),data(:,2),10,'blue','.');
x = linspace(min(data(:,3)),max(data(:,3)),100);
plot(x,polyval(p,x),"Color","r","LineStyle","--","LineWidth",2);
yline(b,"LineWidth",2,"LineStyle","--","Color",'k');
str_breakdown = split(filename,'_');
frequency = str_breakdown(6);
current = str_breakdown(7);
title("Truth vs. Measured for " + frequency + " Hz and " + current + " Amplitude");
xlabel("Encoder (Truth) Data (rad/s)");
ylabel("Gyro (Measured) Data (rad/s)");
%legend("Gyro Data","Encoder Data");
hold off;

end