function [slope] = slopefinder(data_x,data_y,window)
%Slopefinder:   This function finds the slope of a dataset using a forward
%difference scheme.
%   INPUTS:     data_x: a vector of independent variable values
%               data_y: a vector of dependent variable values
%               window: a scalar value of the size of the window over which
%               to find a forward difference
%   Similar to the diff function, the function creates a vector of forward
%   difference values one shorter than the original vector.

slope=zeros(1,length(data_y)-1); %initialize output vector
for i=1:length(data_y)-1    %loop over input vector
    
    y0=data_y(i);           %find y-values
    if (i+window) > length(data_y)
        y1=data_y(end);
    else
        y1=data_y(i+window);
    end

    x0=data_x(i);           %find x-values
    if (i+window) > length(data_x)
        x1=data_x(end);
    else
        x1=data_x(i+window);
    end
slope(i)=(y1-y0)/(x1-x0);   %change in y over change in x
end
end