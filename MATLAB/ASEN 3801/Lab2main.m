% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% LoadASPENData
% Created: 1/31/24

%housecleaning
clc; clear;  close all;  tic;
%load in data
[t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = ...
    LoadASPENData("3801_Sec001_Test2.csv");


%% Question 3: plot both object trajectories in frame N
figure();   hold on;    grid on;
plot3(av_pos_inert(1,:),-av_pos_inert(2,:),-av_pos_inert(3,:),'b');
plot3(tar_pos_inert(1,:),-tar_pos_inert(2,:),-tar_pos_inert(3,:),'--r');
title("Trajectory of Vehicle and Target (Frame N)");
xlabel ("X-Position (mm)");
ylabel ("Y-Position (mm)");
zlabel ("Z-Position (mm)");
legend("Vehicle","Target");
hold off;

%% Question 4: 
figure();
subplot(1,3,1);
plot(t_vec,av_pos_inert(1,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_pos_inert(1,:),"Color",[1 0 0])
subplot(1,3,2);
plot(t_vec,av_pos_inert(2,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_pos_inert(2,:),"Color",[1 0 0])
subplot(1,3,3);
plot(t_vec,av_pos_inert(3,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_pos_inert(3,:),"Color",[1 0 0])

figure();
subplot(1,3,1);
plot(t_vec,av_att(1,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_att(1,:),"Color",[1 0 0])
subplot(1,3,2);
plot(t_vec,av_att(2,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_att(2,:),"Color",[1 0 0])
subplot(1,3,3);
plot(t_vec,av_att(3,:),'Color',[0 0 1])
hold on
plot(t_vec,tar_att(3,:),"Color",[1 0 0])


% %% Question 5:    determine and plot 3-1-3 Euler angles
% 
% av_att_313 = zeros(size(av_att));   %preallocate new variables
% tar_att_313 = zeros(size(tar_att));
% 
% for ii = 1:length(t_vec)            %loop through data
%     %find the DCM for the current 321 Euler angles
%     av_DCM_321 = Rotation321(av_att(:,ii));
%     tar_DCM_321 = Rotation321(tar_att(:,ii));
%     %use EulerAngles313 to back the 313 angles out of the DCM
%     av_att_313(:,ii) = EulerAngles313(av_DCM_321);
%     tar_att_313(:,ii) = EulerAngles313(tar_DCM_321);
% end
% 
% %% Question 6:  position vector of target relative to vehicle in Frame E
% 
% %position of target relative to vehicle is target position - vehicle
% %position
% tar_rel_av = tar_pos_inert - av_pos_inert;
% 
% %% Question 7:  position vector of target relative to vehicle in Frame B
% tar_rel_av_B = zeros(size(tar_rel_av));  %preallocate new variable
% for jj = 1:length(tar_rel_av)   %loop through relative position
%     %rotate each position vector of target relative to vehicle into vehicle
%     %frame using DCM of 321 Euler angles at the time point
%     tar_rel_av_B(:,jj) = RotationMatrix321(av_att(:,jj))*tar_rel_av(:,jj);
% end
% 

