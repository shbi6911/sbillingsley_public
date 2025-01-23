%By:        Shane Billingsley
%Class:     ASEN 2804 Vehicle Design Lab
%Date:      Spring 2023

% Housekeeping
clc; clear; close all;

% Array to help read in thrust files with a for loop
tFiles = ["T01" "T02" "T03" "T04" "T05" "T06" "T07" "T08" "T09" "T10" "T11" "T12" "T13" "T14" "T15"];

Tdata = load('ThrustModel.mat');

thrustCorrection = Tdata.thrustCorrection;
t = Tdata.t;
Thrust = Tdata.Thrust;

% Pre-allocating matrices for efficient use of memory
dataRaw = zeros(14000,15);
m = zeros(1,15);
index = zeros(1,15);
data = zeros(514,15);
dataSmooth = zeros(514,15);
sumRaw = zeros(514,1);
sumSmooth = zeros(514,1);

% Allocate time vector based on average elapsed time
T = (0:513);
Ts = T/1652;

T0 = linspace(0,0.239,396);
thrustMod = T0*40.902;
thrustMod = thrustMod';

thrustMod0 = (ones(118,1))*9.8220943;

correction = zeros(514,1);
correction(1:396) = thrustMod;
correction(397:514) = thrustMod0;



% Loop to extract thrust data from files, discarding first two columns of each and 
% storing the 3rd in the matrix 'testData'. After a file is read in, the maximum
% thrust value and its index are found and stored in 'm' and 'index'. The
% portion of the data starting at the max thrust value and all points through 0.3s 
% are then stored in 'data'. The thrust data is then smoothed and stored in
% 'dataSmooth'
for i = 1:15
    testData = readmatrix('LA_Test_'+tFiles(i)+'_W1000mL_B2000mL');
    dataRaw(:,i) = testData(:,3);
    [m(i),index(i)] = max(dataRaw(:,i));
    data(:,i) = dataRaw(index(i):index(i)+513,i);
    dataSmooth = smoothdata(data,"sgolay");
end

% Loop that averages all raw data and smoothed data vectors
for i = 1:length(T)
    sumRaw(i,:) = sum(data(i,:))/15;
    
    sumSmooth(i,:) = sum(dataSmooth(i,:))/15;
end

% Makes a series of 4 plots of the first raw thrust data set,
% the average of all raw thrust data sets, the first smoothed data set,
% and the averaged smoothed data set



figure(1)
hold on 
grid on
plot(t(1:245),thrustCorrection(1:245),'LineWidth',2)
plot(T/1652,correction)
xlabel("Time [s]")
ylabel("Corrected Thrust Data")
hold off

% figure()
% hold on
% grid on
% plot((1:14000)/1652,dataRaw(:,1:4)*4.448)
% xlabel("Time [s]")
% ylabel("Thrust Data [N]")
% hold off

% figure(2)
% hold on
% grid on
% plot(t, Thrust,"LineWidth",2)
% plot(T/1652,data*4.448)
% xlim([0 0.35])
% title("Raw Thrust vs Time")
% xlabel("Time [s]")
% ylabel("Thrust [N]")
% legend("Model Data","Experimental Data")
% hold off
% 
% figure(3)
% hold on
% grid on
% plot(t, Thrust,"LineWidth",2)
% plot(T/1652,sumRaw*4.448,"r")
% xlim([0 0.35])
% title("Raw Sum Thrust vs Time")
% xlabel("Time [s]")
% ylabel("Thrust [N]")
% legend("Model Data","Experimental Data")
% hold off

figure(4)
hold on
grid on
plot(t, Thrust,"LineWidth",2)
plot(T/1652,(dataSmooth*4.448)-correction)
xlim([0 0.35])
title("Smoothed Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [N]")
legend("Model Data","Experimental Data")
hold off
% 
% figure(5)
% hold on
% grid on
% plot(t, Thrust,"LineWidth",2)
% plot(T/1652,sumSmooth*4.448,"r")
% xlim([0 0.35])
% title("Smoothed Sum Thrust Thrust vs Time")
% xlabel("Time [s]")
% ylabel("Thrust [N]")
% legend("Model Data","Experimental Data")
% hold off

% Loop to plot each of the 15 raw and smoothed data sets for comparison.
% It's computationally expensive and produces 15 windows so it is commented
% out, just un-comment lines 70-85 and run to produce plots

% for i = 1:15
% figure(i)
% grid on
% hold on
% subplot(2,1,1)
% plot(T/1652,dataSmooth(:,i)*4.448,"r")
% title("Smoothed Thrust Data vs Time")
% xlabel("Time [s]")
% ylabel("Thrust [N]")
% subplot(2,1,2)
% plot(T/1652,data(:,i)*4.448)
% title("RawThrust Data vs Time")
% xlabel("Time [s]")
% ylabel("Thrust [N]")
% hold off
% end

% Housekeeping
% clear("i","testData","i","tFiles","dataRaw");