%This portion of code will calculate all work with the large Windtunnel and
%chambered airfoils. 
clc
clear all
close all
Data = readmatrix("ASEN2802_InfiniteWing_FullRange.csv");
PortLoc = readmatrix("ClarkY14_PortLocations.xlsx");
PortLocReal = readmatrix("ClarkY14_PortLocations.xlsx");
%Set up variables to compute a local Cp for each point and angle of attack
pmean = mean(Data(:,3));

%use this new vector as an index to grab one of every different angle of
%attack row. 
DataGrab = 10:20:640;
%ModData = Data(DataGrab,:);
AvgGrab = 1:20:641;
AvgCp = mean(Data(1:20:end,:));
for i = 1:32
    AvgCp(i,:) = mean(Data([AvgGrab(i),AvgGrab(i+1)-1],:));
    i = i+1;
end
ModData = AvgCp;

%Cp matrix will now compute CP of each data point
Cp = ModData(:,15:30)./((1/2)*pmean*(ModData(:,4)).^2);

%Begin the process of graphing all figures and their proper labels. 

f1 = figure('units','normalized','outerposition',[0 0 1 1]);


%Create modified PortLoc matrix that normalizes the chord. 
PortLocMod = (PortLoc(1:17,:))/PortLoc(10,3);

%Add the fisrt colum of CP to the end so that it forms closed loop for
%intergration
Cpadd= Cp(:,1);
Cp(:,17) = Cpadd;
%Reorganize data so that TE point can be calculated and placed at the
%correct xlocation
Cp(:,18) = 1;
Cp = Cp(:,[1 2 3 4 5 6 7 8 9 18 10 11 12 13 14 15 16 17]);

%Computing the Trailing Edge point for each angle of attack. 
for i =1:length(DataGrab)
    %calculate line of best fit for top
    mtop=(Cp(i,8)-Cp(i,9))/(PortLocMod(8,3)-PortLocMod(9,3));
    %ADDING +1 to all rows of PortLocMod
    b = Cp(i,9)-(mtop*PortLocMod(9,3));
    Ytop = mtop+b;
    %line of best fit for bottom
    mbot = (Cp(i,12)-Cp(i,11))/(PortLocMod(12,3)-PortLocMod(11,3));
    b = Cp(i,11)-(mbot*PortLocMod(11,3));
    Ybot = mbot + b;
    %take average of both points. 
    Yavg = (Ytop+Ybot)/2;
    Cp(i,10) = Yavg;
    %increase index
    i = i+1;
end

PortLocMod(18,:) = [18 0 0 0];
PortLoc(18,:) = [18 0 0 0];
%Now plot with new averaged TE point
for i=1:32
    
    subplot(6,6,i);
    plot(PortLocMod(:,3),Cp(i,:))
    set(gca, 'YDir','reverse')
    ylim([-3.5 1])
    xlim([0 1.1])
    title("Airspeed = " + ModData(i,4) + "[m/s], AoA = "+ ModData(i,8) +"[deg]");
    xlabel 'Normalized Chord';
    ylabel 'Cp';
  
    i = i+1;
    
end


%Plot the base model of the Clark Y-14 airfoil
subplot(6,6,36);
plot(PortLoc(:,3),PortLoc(:,4),"-o","color","blue")
ylim([-0.3 0.5])
xlim([-0.3 3.7])
title("Clark Y-14 Airfoil");
xlabel 'y-axis (in)';
ylabel 'z-axis (in)';


Cl = zeros(32,1);
Cd = zeros(32,1);

for i=-15:16
    k=i+16;
    [Cl(k),Cd(k)] = ClCd_Calc(PortLoc,Cp,i);
end
% AoA = 9
% row = AoA +16;
% IntegralN = trapz(PortLoc(:,3)./(3.5031),Cp(row,:));
% IntegralA = trapz(PortLoc(:,4),Cp(row,:));
% Cn = -1*IntegralN;
% Ca = IntegralA;
% Cl = Cn*cos(AoA)-Ca*sin(AoA);
% Cd = Cn*sin(AoA)+Ca*cos(AoA);




X = [-15:1:16];
f2 = figure;

plot(X,Cl);
set(gca, 'YDir','normal')
xlabel("Coefficient of Lift [Cl]")
ylabel("Angle of Attack [AoA]")
title("Cl across AoA")
f3 = figure;
plot(X,Cd);
xlabel("Coefficient of Drag [Cd]")
ylabel("Angle of Attack [AoA]")
title("Cd across AoA")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cl and Cd functions 

function [Cl,Cd] = ClCd_Calc(Xdata,Ydata,AoA)%add "constants"?
row = AoA +16;
IntegralN = -(1/max(Xdata(:,3)))*trapz(Xdata(:,3),Ydata(row,:));
IntegralA = (1/max(Xdata(:,3)))*trapz(Xdata(:,4),Ydata(row,:));
Cn = 1*IntegralN;
Ca = 1*IntegralA;
Cl = Cn*cos(AoA*pi/180)-Ca*sin(AoA*pi/180);
Cd = Cn*sin(AoA*pi/180)+Ca*cos(AoA*pi/180);

end
