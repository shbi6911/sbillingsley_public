%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

labdata = readmatrix("labdata3");
time = labdata(:,1);
angle = labdata(:,2);
angvelocity = labdata(:,4);
reference = labdata(:,6);
voltage = labdata(:,7);

plot (time, voltage);
xlabel ("Time (ms)");
ylabel ("Motor Voltage (V)");
title ("Voltage Applied to the Motor Over Time");
