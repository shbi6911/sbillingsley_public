function [yFit] = PiecewiseLeastFit(M,B,xFit)
%PiecewiseLeastFit calculates the y values of the regression line given M B and the xFit vector
%   uses a loop that iterates four times through each range of x values, it
%   creates a logical vector to only get the values of x within the current
%   range and stores it into vector x, these same values are written into a
%   temporary yfit vector but are overwritten in the nested loop where each
%   y value is calculated. the temporary vector is then added to yfit which
%   increases the size of the vector until it has every y value
yFit = zeros(0,0);
xranges = [-5, 10, 15, 20, 100];
for ii = 1:4
    xlog = (xFit >= xranges(ii)) & (xFit < xranges(ii+1));
    x = xFit(xlog);
    yFittemp = xFit(xlog);
    for jj = 1:length(x)
        yFittemp(jj) = M(ii)*x(jj) + B(ii);
    end
    yFit = [yFit,yFittemp];
end
end