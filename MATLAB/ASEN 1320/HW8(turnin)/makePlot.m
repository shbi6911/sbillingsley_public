function [] = makePlot(xData,yData,xFit,yFit)
%input: xData, yData, xFit, yFit
%no output
%makePlot plots both the least squares regression lines and the initial
%data
%   It plots the initial x and y data first using a dotted blue line, then
%   on the same plot, the fit lines are added with a solid red line. The x
%   is labeled Angle of attack (degrees) and the y is labeled Coefficient
%   of lift, C_l. a legend then labels both sets of data
plot(xData, yData, '.b'); %dotted blue line
hold on;
plot(xFit,yFit, '-r'); %solid red line
ylabel('Coefficient of lift, C_l');
xlabel('Angle of attack (degrees)');
legend('Experimental Data','Piecewise Linear Fit');
end