%By:        Shane Billingsley
%Class:     ASEN 2804 Vehicle Design Lab
%Date:      Spring 2023

clc; clear; close all;

model_250 = load('ThrustModel_250.mat');
thrust_250 = model_250.Thrust(1:353);
t_250 = model_250.t(1:353);
fix_250 = model_250.thrustCorrection(1:229);

model_313 = load('ThrustModel_313.mat');
thrust_313 = model_313.Thrust(1:375);
t_313 = model_313.t(1:375);
fix_313 = model_313.thrustCorrection(1:219);

model_375 = load('ThrustModel_375.mat');
thrust_375 = model_375.Thrust(1:355);
t_375 = model_375.t(1:355);
fix_375 = model_375.thrustCorrection(1:224);

model_438 = load('ThrustModel_438.mat');
thrust_438 = model_438.Thrust(1:386);
t_438 = model_438.t(1:386);
fix_438 = model_438.thrustCorrection(1:276);

model_500 = load('ThrustModel_500.mat');
thrust_500 = model_500.Thrust(1:330);
t_500 = model_500.t(1:330);
fix_500 = model_500.thrustCorrection(1:208);

model_563 = load('ThrustModel_563.mat');
thrust_563 = model_563.Thrust(1:335);
t_563 = model_563.t(1:335);
fix_563 = model_563.thrustCorrection(1:233);

model_625 = load('ThrustModel_625.mat');
thrust_625 = model_625.Thrust(1:257);
t_625 = model_625.t(1:257);
fix_625 = model_625.thrustCorrection(1:224);

load('PressureModel.mat');

dat_250 = ["250_2" "250_3" "250_4"];
dat_313 = ["313_1" "313_2" "313_3" "313_4" "313_5"];
dat_375 = ["375_1" "375_2" "375_3" "375_4" "375_5"];
dat_438 = ["438_1" "438_2" "438_3" "438_4" "438_5"];
dat_500 = ["500_1" "500_2" "500_3" "500_4" "500_5"];
dat_563 = ["563_1" "563_2" "563_3" "563_4" "563_5"];
dat_625 = ["625_1" "625_2" "625_3" "625_4" "625_5"];
data250 = zeros(14000,3);
data313 = zeros(14000,5);
data375 = zeros(14000,5);
data438 = zeros(14000,5);
data500 = zeros(14000,5);
data563 = zeros(14000,5);
data625 = zeros(14000,5);
max_250 = zeros(3,1);
index_250 = zeros(3,1);
max_313 = zeros(5,1);
index_313 = zeros(5,1);
max_375 = zeros(5,1);
index_375 = zeros(5,1);
max_438 = zeros(5,1);
index_438 = zeros(5,1);
max_500 = zeros(5,1);
index_500 = zeros(5,1);
max_563 = zeros(5,1);
index_563 = zeros(5,1);
max_625 = zeros(5,1);
index_625 = zeros(5,1);
test_250 = zeros(251,3);
test_313 = zeros(251,5);
test_375 = zeros(251,5);
test_438 = zeros(251,5);
test_500 = zeros(251,5);
test_563 = zeros(251,5);
test_625 = zeros(251,5);
smooth_250 = zeros(251,3);
smooth_313 = zeros(251,5);
smooth_375 = zeros(251,5);
smooth_438 = zeros(251,5);
smooth_500 = zeros(251,5);
smooth_563 = zeros(251,5);
smooth_625 = zeros(251,5);
t = 0:250;
for i = 1:3
data_250 = readmatrix("1250mL Bottle/Test_"+dat_250(i));
data250(:,i) = data_250(5:end,3);
[max_250(i) index_250(i)] = max(data250(:,i));
test_250(:,i) = data250(index_250(i):index_250(i)+250,i);
smooth_250(:,i) = smoothdata(test_250(:,i),"sgolay");
end
for i = 1:5
data_313 = readmatrix("1250mL Bottle/Test_"+dat_313(i));
data313(:,i) = data_313(5:end,3);
[max_313(i) index_313(i)] = max(data313(:,i));
test_313(:,i) = data313(index_313(i):index_313(i)+250,i);
smooth_313(:,i) = smoothdata(test_313(:,i),"sgolay");
data_375 = readmatrix("1250mL Bottle/Test_"+dat_375(i));
data375(:,i) = data_375(5:end,3);
[max_375(i) index_375(i)] = max(data375(:,i));
test_375(:,i) = data375(index_375(i):index_375(i)+250,i);
smooth_375(:,i) = smoothdata(test_375(:,i),"sgolay");
data_438 = readmatrix("1250mL Bottle/Test_"+dat_438(i));
data438(:,i) = data_438(5:end,3);
[max_438(i) index_438(i)] = max(data438(:,i));
test_438(:,i) = data438(index_438(i):index_438(i)+250,i);
smooth_438(:,i) = smoothdata(test_438(:,i),"sgolay");
data_500 = readmatrix("1250mL Bottle/Test_"+dat_500(i));
data500(:,i) = data_500(5:end,3);
[max_500(i) index_500(i)] = max(data500(:,i));
test_500(:,i) = data500(index_500(i):index_500(i)+250,i);
smooth_500(:,i) = smoothdata(test_500(:,i),"sgolay");
data_563 = readmatrix("1250mL Bottle/Test_"+dat_563(i));
data563(:,i) = data_563(5:end,3);
[max_563(i) index_563(i)] = max(data563(:,i));
test_563(:,i) = data563(index_563(i):index_563(i)+250,i);
smooth_563(:,i) = smoothdata(test_563(:,i),"sgolay");
data_625 = readmatrix("1250mL Bottle/Test_"+dat_625(i));
data625(:,i) = data_625(5:end,3);
[max_625(i) index_625(i)] = max(data625(:,i));
test_625(:,i) = data625(index_625(i):index_625(i)+250,i);
smooth_625(:,i) = smoothdata(test_625(:,i),"sgolay");
end
% figure(1)
% 
% grid on
% subplot(2,1,1)
% plot(t/1652,test_250*4.448)
% title('1.25L Bottle Thrust 20/80 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_250*4.448)
% plot(t_250,thrust_250,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(2)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_313*4.448)
% title('1.25L Bottle Thrust 25/75 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_313*4.448)
% plot(t_313,thrust_313,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(3)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_375*4.448)
% title('1.25L Bottle Thrust 30/70 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_375*4.448)
% plot(t_375,thrust_375,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(4)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_438*4.448)
% title('1.25L Bottle Thrust 35/65 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_438*4.448)
% plot(t_438,thrust_438,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(5)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_500*4.448)
% title('1.25L Bottle Thrust 40/60 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_500*4.448)
% plot(t_500,thrust_500,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(6)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_563*4.448)
% title('1.25L Bottle Thrust 45/55 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_563*4.448)
% plot(t_563,thrust_563,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
% 
% figure(7)
% grid on
% subplot(2,1,1)
% plot(t/1652,test_625*4.448)
% title('1.25L Bottle Thrust 50/50 Ratio')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% subplot(2,1,2)
% hold on
% plot(t/1652,smooth_625*4.448)
% plot(t_625,thrust_625,"LineWidth",2)
% title('Smoothed vs Model')
% xlabel('Time [s]')
% ylabel('Thrust [N]')
% hold off
clear("dat_250","dat_313","dat_375","dat_438","dat_500","dat_563","dat_625", ...
    "data250","data313","data375","data438","data500","data563","data625", ...
    "data_250","data_313","data_375","data_438","data_500","data_563","data_625");

tests = ["250","313","375","438","500","563","625"];
end_time_new = [0.1 0.019 0.025 0.04 0.05 0.065 0.075];
C_dis_global_1250ml=zeros(1,6);
time_exp = t/1652;
const.At = pi * 0.021^2 * (1/4);
const.Pa = 83403.57;
indices = zeros(1,length(end_time));
figure(8);
hold on; grid on;

for i = 1:length(end_time)
    index = find((time_exp > end_time_new(i)),1,'first');
    indices(i)=index;
    dstring = ("data"+tests(i));
    data_exp=(eval("smooth_"+tests(i))*4.448);
    data_exp=data_exp(1:index,:);
    fit_exp = polyval(P_all(i,:),time_exp(1:index));
    value_exp = (2*const.At*(fit_exp-const.Pa))';
    C_dis_all.(dstring) = data_exp./value_exp;
    %plot(time_exp(1:index),C_dis_all.(dstring));
end

C_dis_edited.data313 = C_dis_all.data313(:,1:4);
C_dis_edited.data375 = C_dis_all.data375(:,2:5);
C_dis_edited.data438 = C_dis_all.data438(:,2:5);
C_dis_edited.data500 = C_dis_all.data500;
C_dis_edited.data563 = C_dis_all.data563;
C_dis_edited.data625 = C_dis_all.data625(:,1:3);

for i = 1:length(end_time_new)-1
    dstring = ("data"+tests(i+1));
    plot(time_exp(1:indices(i+1)),C_dis_edited.(dstring));
    C_dis_global_1250ml(i) = mean(C_dis_edited.(dstring),'all');
    
end
title("Experimental Discharge Coefficient 1250 mL Bottle");
ylim([0.5 1.5]);
xlabel("Time [s]");
ylabel("Discharge Coefficient");
yline(mean(C_dis_global_1250ml),'LineWidth',2,'Color','red',label='Mean Value');