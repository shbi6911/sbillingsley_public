%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

function [fit, fitline] = leastSquares(x,y)
N = length(x);
H = [ones(N,1) x];
fit = inv(H'*H)*(H'*y);
fitline = fit(2)*x + fit(1);
