function [yData] = GenerateData(xData)
%GenerateData takes an input of xData(angles of attack) and calculates the 
% yData at each of those angles
%   loads the Pcoef.dat file containing the given coefficients needed to
%   find the yData and stores it into polycoeffvector. Uses the polyval
%   function with the coefficients and xData to find the ydata(coefficient
%   of lift). To add noise to the data, we create a noise vector of the
%   same size as the xData with values between -0.3 and 0.3. This is then
%   added to the ydata and the noisy ydata is returned.
polycoeffvector = load('Pcoef.dat', '-ascii'); %expecting ascii data
yDatainitial = polyval(polycoeffvector, xData);
xDatalength = length(xData);
rng(uint64(now*1000))
noise = -0.3 + (0.6)*rand(1,xDatalength);
yData = yDatainitial + noise;
end