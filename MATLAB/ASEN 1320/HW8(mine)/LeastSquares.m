%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function [M,B] = LeastSquares(X,Y)
%LEASTSQUARES computes the slope and y-intercept for the least squares
%best-fit line for a given set of X and Y values
%Input: X and Y vectors, which must be the same length
%Output: a value M, which is the slope of the best fit line, and a value B,
%which is the y-intercept of the best fit line.
A = sum(X);
B = sum(Y);
temp = X.*Y;
C = sum(temp);
square = X.^2;
D = sum(square);
N = length(X);
M = ((A*B) - (N*C))/((A^2)-(N*D));
B = ((A*C) - (B*D))/((A^2)-(N*D));
end