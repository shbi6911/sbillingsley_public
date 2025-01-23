%By:        Shane Billingsley
%Class:     ASEN 2804 Vehicle Design Lab
%Date:      Spring 2023

% Housekeeping
clc; clear; close all;

% Array to help read in thrust files with a for loop
tFiles = ["T01" "T02" "T03" "T04" "T05" "T06" "T07" "T08" "T09" "T10" "T11" "T12" "T13" "T14" "T15"];

Tdata08 = load('ThrustModel_CD08.mat');
Tdata075 = load('ThrustModel_CD075.mat');

thrustCorrection075 = Tdata075.thrustCorrection;
t075 = Tdata075.t(1:377);
Thrust075 = Tdata075.Thrust(1:377);
fix75 = Tdata075.thrustCorrection(1:514);


thrustCorrection08 = Tdata08.thrustCorrection;
t08 = Tdata08.t(1:365);
Thrust08 = Tdata08.Thrust(1:365);
fix08 = Tdata08.thrustCorrection(1:514);



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

TotalImpulseModel_CD08 = trapz(Thrust08);
TotalImpulseModel_CD075 = trapz(Thrust075);

for i = 1:15
TotalImpulseExperiment08(i) = (trapz(((dataSmooth(:,i))*4.448)-fix08));
TotalImpulseExperiment75(i) = (trapz(((dataSmooth(:,i))*4.448)-fix75));
PercentDifference08(i) = 100*abs(TotalImpulseModel_CD08-TotalImpulseExperiment08(i))/(TotalImpulseModel_CD08);
PercentDifference075(i) = 100*abs(TotalImpulseModel_CD075-TotalImpulseExperiment75(i))/(TotalImpulseModel_CD075);
end

totalDelta08 = sum(PercentDifference08)/15;
totalDelta75 = sum(PercentDifference075)/15;

figure(1)
hold on
grid on
plot(t08, Thrust08,"LineWidth",2)
plot(T/1652,(dataSmooth*4.448)-fix08)
title("Cd = 0.8: Smoothed Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [N]")
legend("2L Model Data","Experimental Data")
hold off

figure(2)
hold on
grid on
plot(t075, Thrust075,"LineWidth",2)
plot(T/1652,(dataSmooth*4.448)-fix75)
xlim([0 0.35])
title("Cd = 0.75: Smoothed Thrust vs Time")
xlabel("Time [s]")
ylabel("Thrust [N]")
legend("2L Model Data","Experimental Data")
hold off