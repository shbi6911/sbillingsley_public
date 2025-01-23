%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function [] = makePlot(xData,yData,xFit,yFit)
%MAKEPLOT creates a plot with two datasets and labels it
%Input: two vectors with x and y values for experimental data, and two 
% vectors with x and y values for a least squares best fit line.
%Output: none (prints a plot)
plot (xData, yData, 'k.');
hold on;
plot (xFit, yFit, 'r-');
xlabel ('Angle of Attack (degrees)');
ylabel ('Coefficient of lift, C_l');
legend ('Experimental Data', 'Piecewise Least Fit');


end