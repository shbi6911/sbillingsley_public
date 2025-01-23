clc; clear; close all;

model_250 = load('ThrustModel_250_CD84.mat');
thrust_250 = model_250.Thrust(1:363);
t_250 = model_250.t(1:363);
mod250 = model_250.thrustCorrection(1:251);

model_313 = load('ThrustModel_313_CD84.mat');
thrust_313 = model_313.Thrust(1:343);
t_313 = model_313.t(1:343);
mod313 = model_313.thrustCorrection(1:251);


model_375 = load('ThrustModel_375_CD84.mat');
thrust_375 = model_375.Thrust(1:358);
t_375 = model_375.t(1:358);
mod375 = model_375.thrustCorrection(1:251);


model_438 = load('ThrustModel_438_CD84.mat');
thrust_438 = model_438.Thrust(1:354);
t_438 = model_438.t(1:354);
mod438 = model_438.thrustCorrection(1:251);

model_500 = load('ThrustModel_500_CD84.mat');
thrust_500 = model_500.Thrust(1:338);
t_500 = model_500.t(1:338);
mod500 = model_500.thrustCorrection(1:251);


model_563 = load('ThrustModel_563_CD84.mat');
thrust_563 = model_563.Thrust(1:324);
t_563 = model_563.t(1:324);
mod563 = model_563.thrustCorrection(1:251);


model_625 = load('ThrustModel_625_CD84.mat');
thrust_625 = model_625.Thrust(1:274);
t_625 = model_625.t(1:274);
mod625 = model_625.thrustCorrection(1:251);

t = 0:250;


dat_250 = ["250_2" "250_3" "250_4"];
dat_313 = ["313_1" "313_2" "313_3" "313_4" "313_5"];
dat_375 = ["375_1" "375_2" "375_3" "375_4" "375_5"];
dat_438 = ["438_1" "438_2" "438_3" "438_4" "438_5"];
dat_500 = ["500_1" "500_2" "500_3" "500_4" "500_5"];
dat_563 = ["563_1" "563_2" "563_3" "563_4" "563_5"];
dat_625 = ["625_1" "625_2" "625_3" "625_4" "625_5"];


for i = 1:3
data_250 = readmatrix("1250mL Bottle/Test_"+dat_250(i));
data250(:,i) = data_250(5:end,3);
[max_250(i) index_250(i)] = max(data250(:,i));
test_250(:,i) = data250(index_250(i):index_250(i)+250,i);
smooth_250(:,i) = (smoothdata(test_250(:,i),"sgolay")*4.448)-mod250;
end
for i = 1:5
data_313 = readmatrix("1250mL Bottle/Test_"+dat_313(i));
data313(:,i) = data_313(5:end,3);
[max_313(i) index_313(i)] = max(data313(:,i));
test_313(:,i) = data313(index_313(i):index_313(i)+250,i);
smooth_313(:,i) = (smoothdata(test_313(:,i),"sgolay")*4.448)-mod313;
data_375 = readmatrix("1250mL Bottle/Test_"+dat_375(i));
data375(:,i) = data_375(5:end,3);
[max_375(i) index_375(i)] = max(data375(:,i));
test_375(:,i) = data375(index_375(i):index_375(i)+250,i);
smooth_375(:,i) = (smoothdata(test_375(:,i),"sgolay")*4.448)-mod375;
data_438 = readmatrix("1250mL Bottle/Test_"+dat_438(i));
data438(:,i) = data_438(5:end,3);
[max_438(i) index_438(i)] = max(data438(:,i));
test_438(:,i) = data438(index_438(i):index_438(i)+250,i);
smooth_438(:,i) = (smoothdata(test_438(:,i),"sgolay")*4.448)-mod438;
data_500 = readmatrix("1250mL Bottle/Test_"+dat_500(i));
data500(:,i) = data_500(5:end,3);
[max_500(i) index_500(i)] = max(data500(:,i));
test_500(:,i) = data500(index_500(i):index_500(i)+250,i);
smooth_500(:,i) = (smoothdata(test_500(:,i),"sgolay")*4.448)-mod500;
data_563 = readmatrix("1250mL Bottle/Test_"+dat_563(i));
data563(:,i) = data_563(5:end,3);
[max_563(i) index_563(i)] = max(data563(:,i));
test_563(:,i) = data563(index_563(i):index_563(i)+250,i);
smooth_563(:,i) = (smoothdata(test_563(:,i),"sgolay")*4.448)-mod563;
data_625 = readmatrix("1250mL Bottle/Test_"+dat_625(i));
data625(:,i) = data_625(5:end,3);
[max_625(i) index_625(i)] = max(data625(:,i));
test_625(:,i) = data625(index_625(i):index_625(i)+250,i);
smooth_625(:,i) = (smoothdata(test_625(:,i),"sgolay")*4.448)-mod625;
end
It_250 = zeros(1,4);
It_313 = zeros(1,6);
It_375 = zeros(1,6);
It_438 = zeros(1,6);
It_500 = zeros(1,6);
It_563 = zeros(1,6);
It_625 = zeros(1,4);

for i = 1:6
if i == 6
It_313(i) = trapz(thrust_313);
It_375(i) = trapz(thrust_375);
It_438(i) = trapz(thrust_438);
It_500(i) = trapz(thrust_500);
It_563(i) = trapz(thrust_563);
It_625(i) = trapz(thrust_625);
else
It_313(i) = trapz(smooth_313(:,i));
It_375(i) = trapz(smooth_375(:,i));
It_438(i) = trapz(smooth_438(:,i));
It_500(i) = trapz(smooth_500(:,i));
It_563(i) = trapz(smooth_563(:,i));
It_625(i) = trapz(smooth_625(:,i));

end

end

for j = 1:5
PercentDiff_313(j) = abs(100*(((It_313(6)-It_313(j))/It_313(6))));
PercentDiff_375(j) = abs(100*(((It_375(6)-It_375(j))/It_375(6))));
PercentDiff_438(j) = abs(100*(((It_438(6)-It_438(j))/It_438(6))));
PercentDiff_500(j) = abs(100*(((It_500(6)-It_500(j))/It_500(6))));
PercentDiff_563(j) = abs(100*(((It_563(6)-It_563(j))/It_563(6))));
end

for i = 1:4


if i == 4
It_250(i) = trapz(thrust_250);
It_625(i) = trapz(thrust_625);
else
It_250(i) = trapz(smooth_250(:,i));
It_625(i) = trapz(smooth_625(:,i));
end

end

for j = 1:3
PercentDiff_250(j) = abs(100*(((It_250(4)-It_250(j))/It_250(4))));
PercentDiff_625(j) = abs(100*(((It_625(4)-It_625(j))/It_625(4))));
end
figure(1)

grid on
subplot(2,1,1)
plot(t/1652,test_250*4.448)
title('1.25L Bottle Thrust 20/80 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_250)
plot(t_250,thrust_250,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(2)
grid on
subplot(2,1,1)
plot(t/1652,test_313(:,1:end)*4.448)
title('1.25L Bottle Thrust 25/75 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_313(:,1:4))
plot(t_313,thrust_313,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(3)
grid on
subplot(2,1,1)
plot(t/1652,test_375(:,2:end)*4.448)
title('1.25L Bottle Thrust 30/70 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_375(:,2:end))
plot(t_375,thrust_375,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(4)
grid on
subplot(2,1,1)
plot(t/1652,test_438(:,2:end)*4.448)
title('1.25L Bottle Thrust 35/65 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_438(:,2:end))
plot(t_438,thrust_438,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(5)
grid on
subplot(2,1,1)
plot(t/1652,test_500*4.448)
title('1.25L Bottle Thrust 40/60 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_500)
plot(t_500,thrust_500,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(6)
grid on
subplot(2,1,1)
plot(t/1652,test_563*4.448)
title('1.25L Bottle Thrust 45/55 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_563)
plot(t_563,thrust_563,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off

figure(7)
grid on
subplot(2,1,1)
plot(t/1652,test_625(:,1:3)*4.448)
title('1.25L Bottle Thrust 50/50 Ratio')
xlabel('Time [s]')
ylabel('Thrust [N]')
subplot(2,1,2)
hold on
plot(t/1652,smooth_625(:,1:3))
plot(t_625,thrust_625,"LineWidth",2)
title('Smoothed vs Model')
xlabel('Time [s]')
ylabel('Thrust [N]')
hold off
clear("dat_250","dat_313","dat_375","dat_438","dat_500","dat_563","dat_625", ...
    "data250","data313","data375","data438","data500","data563","data625", ...
    "data_250","data_313","data_375","data_438","data_500","data_563","data_625");