function [m,b] = LeastSquares(xData, yData)
%LeastSquares calculates the slope and y intercept of a least squares regression line for equal lengths of x and y data
%   Calculates values for the length of the X vector, the sum of both data,
%   the sum of the vectors multiplied together(not matrix math), and the
%   sum of the squared xdata. It then uses all these values to calculate
%   the slope and y intercept which are both returned
N = length(xData);
A = sum(xData);
B = sum(yData);
XtimesY = xData.*yData;
C = sum(XtimesY);
Xsquared = xData.^2;
D = sum(Xsquared);
m = (A*B - N*C)/(A^2 - N*D);
b = (A*C - B*D)/(A^2 - N*D);
end