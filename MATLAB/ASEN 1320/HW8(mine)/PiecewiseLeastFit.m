%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function [yFit] = PiecewiseLeastFit(M,B,xFit)
%PIECEWISELEASTFIT creates a vector of y-values to match a given vector of
%x-values, according to a piecewise linear function
%Input: a vector of four slope values, a vector of four y-intercepts, and a
%vector of x-values
%Output: a vector of y-values
yFit = zeros(1,5001);
for ii = 1:length(xFit)
    if xFit(ii)<10
        yFit(ii)=(M(1)*xFit(ii)) + B(1);
    elseif xFit(ii)>= 10 && xFit(ii) < 15
        yFit(ii)=(M(2)*xFit(ii)) + B(2);
    elseif xFit(ii)>= 15 && xFit(ii) < 20
        yFit(ii)=(M(3)*xFit(ii)) + B(3);
    elseif xFit(ii)>=20
        yFit(ii)=(M(4)*xFit(ii)) + B(4);
    end
end