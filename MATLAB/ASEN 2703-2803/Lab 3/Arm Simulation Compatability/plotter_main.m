%Bobby Hodgkinson
%ASEN 2003 Lab 6
%Simulation data loader

close all;
%Time (ms)	Hubangle(Theta in rad)	Tip Deflection(m)	Hub Angular Velocity (rads/s)	Tip Velocity (m/s)	Position Reference (rad)	Output Voltage (V)	K1 (Hub Prop)	K2 (Tip Prop)	K3 (Hub Deriv)	K4 (Tip Deriv) = data.simout.Data;
Time = data.simout.Data(:,1);
Hub_angle = data.simout.Data(:,2);
Tip = data.simout.Data(:,3);
Hub_velocity = data.simout.Data(:,4);
Tip_velocity = data.simout.Data(:,5);
Pos_ref = data.simout.Data(:,6);
Voltage = data.simout.Data(:,7);
K1 = data.simout.Data(:,8);
K2 = data.simout.Data(:,9);
K3 = data.simout.Data(:,10);
K4 = data.simout.Data(:,11);

figure
plot(Time,Hub_angle,'k')
hold on
plot(Time,Pos_ref,'r')
title(['Hub angle and reference for K1 =',num2str(K1(1)),' K2 = ',num2str(K2(1))])
xlabel('Time (ms)')
ylabel('Angle (rads)')
