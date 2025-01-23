% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% LoadASPENData
% Created: 1/31/24

%housecleaning
clc; clear;  close all;  tic;
%load in data
[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = ...
    LoadASPENData("3801_Sec001_Test2.csv");

% scatter plot dot size
sz = 5;
%% Question 3: plot both object trajectories in frame N
f = figure();   hold on;    grid on;
plot3(av_pos_inert(1,:),-av_pos_inert(2,:),-av_pos_inert(3,:),'b');
plot3(tar_pos_inert(1,:),-tar_pos_inert(2,:),-tar_pos_inert(3,:),'--r');
title("Trajectory of Vehicle and Target (Frame N)");
xlabel ("X-Position (mm)");
ylabel ("Y-Position (mm)");
zlabel ("Z-Position (mm)");
legend("Vehicle","Target");
hold off;

saveas(f,'Q3','png')
%% Question 4: 
% Position Plots
f = figure();
subplot(3,1,1);
plot(t_vec,av_pos_inert(1,:),'Color',[0 0 1])
title("Position in E Frame")
hold on
plot(t_vec,tar_pos_inert(1,:),"Color",[1 0 0])
xlabel("Time [s]")
ylabel("x-position [mm]")
subplot(3,1,2);
plot(t_vec,av_pos_inert(2,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_pos_inert(2,:),"Color",[1 0 0])
xlabel("Time [s]")
ylabel("y-position [mm]")
subplot(3,1,3);
plot(t_vec,av_pos_inert(3,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_pos_inert(3,:),"Color",[1 0 0])
xlabel("Time [s]")
ylabel("z-position [mm]")

saveas(f,'Q4Position','png')

% Euler Angle Plots
f = figure();
subplot(3,1,1);
scatter(t_vec,av_att(1,:),sz,'filled','Color',[0 0 1])
title("3-2-1 Euler Angles")
hold on
scatter(t_vec,tar_att(1,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Alpha [deg]")
subplot(3,1,2);
scatter(t_vec,av_att(2,:),sz,'filled','Color',[0 0 1])
hold on
scatter(t_vec,tar_att(2,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Beta [deg]")
subplot(3,1,3);
scatter(t_vec,av_att(3,:),sz,'filled','Color',[0 0 1])
hold on
scatter(t_vec,tar_att(3,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Gamma [deg]")

saveas(f,'Q4Angles','png')
%% Question 5:    determine and plot 3-1-3 Euler angles

av_att_313 = zeros(size(av_att));   %preallocate new variables
tar_att_313 = zeros(size(tar_att));

for ii = 1:length(t_vec)            %loop through data
    %find the DCM for the current 321 Euler angles
    av_DCM_321 = RotationMatrix321(av_att(:,ii));
    tar_DCM_321 = RotationMatrix321(tar_att(:,ii));
    %use EulerAngles313 to back the 313 angles out of the DCM
    av_att_313(:,ii) = EulerAngles313(av_DCM_321);
    tar_att_313(:,ii) = EulerAngles313(tar_DCM_321);
end

% Plot
f = figure();
subplot(3,1,1);
scatter(t_vec,av_att_313(1,:),sz,'filled','Color',[0 0 1])
title("3-1-3 Euler Angles")
hold on
scatter(t_vec,tar_att_313(1,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Alpha [deg]")
subplot(3,1,2);
scatter(t_vec,av_att_313(2,:),sz,'filled','Color',[0 0 1])
hold on
scatter(t_vec,tar_att_313(2,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Beta [deg]")
subplot(3,1,3);
scatter(t_vec,av_att_313(3,:),sz,'filled','Color',[0 0 1])
hold on
scatter(t_vec,tar_att_313(3,:),sz,'filled',"Color",[1 0 0])
xlabel("Time [s]")
ylabel("Gamma [deg]")

saveas(f,'Q5','png')
%% Question 6:  position vector of target relative to vehicle in Frame E

%position of target relative to vehicle is target position - vehicle
%position
tar_rel_av = tar_pos_inert - av_pos_inert;

% Plot
f = figure();
subplot(3,1,1);
plot(t_vec,tar_rel_av(1,:),'Color',[1 0 1])
title("Realtive Position in E Frame")
xlabel("Time [s]")
ylabel("x-position [mm]")
subplot(3,1,2);
plot(t_vec,tar_rel_av(2,:),"Color",[1 0 1])
xlabel("Time [s]")
ylabel("y-position [mm]")
subplot(3,1,3);
plot(t_vec,tar_rel_av(3,:),"Color",[1 0 1])
xlabel("Time [s]")
ylabel("z-position [mm]")

saveas(f,'Q6','png')

%% Question 7:  position vector of target relative to vehicle in Frame B
tar_rel_av_B = zeros(size(tar_rel_av));  %preallocate new variable
for jj = 1:length(tar_rel_av)   %loop through relative position
    %rotate each position vector of target relative to vehicle into vehicle
    %frame using DCM of 321 Euler angles at the time point
    tar_rel_av_B(:,jj) = RotationMatrix321(av_att(:,jj))*tar_rel_av(:,jj);
end

% Plot
f = figure();
subplot(3,1,1);
plot(t_vec,tar_rel_av_B(1,:),'Color',[1 0 1])
title("Realtive Position in B Frame")
xlabel("Time [s]")
ylabel("x-position [mm]")
subplot(3,1,2);
plot(t_vec,tar_rel_av_B(2,:),"Color",[1 0 1])
xlabel("Time [s]")
ylabel("y-position [mm]")
subplot(3,1,3);
plot(t_vec,tar_rel_av_B(3,:),"Color",[1 0 1])
xlabel("Time [s]")
ylabel("z-position [mm]")

saveas(f,'Q7','png')
