%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function [yData] = GenerateData(xData)
%GENERATEDATA This function generates a polynomial curve and adds random
%noise
%Input: a vector of x-values
%Output: a corresponding vector of y-values
% The function evaluates a polynomial on xData, whose coefficients are read from file
% Pcoef.dat, and generates yData.  It then iterates through yData adding
% random values between -0.3 and 0.3 to the results.
load Pcoef.dat;
yData = polyval (Pcoef, xData);
low = -0.3;
high = 0.3;
rng(uint64(now*1000));
for ii = 1:length(yData)
    temp = rand*(high-low)+low;
    yData(ii) = yData(ii)+(temp);
end