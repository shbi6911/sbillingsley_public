function [tss,steadyvalue] = steadystate(data_x,data_y,window,threshold)
%Steadystate:   This function finds the time to steady state and the steady
%state value of input data
%   INPUTS:     data_x: a vector of independent variable values
%               data_y: a vector of dependent variable values
%               window: a scalar value of the size of the window over which
%               to smooth the data and find a forward difference
%               threshold: when the slope drops below this value this is
%               considered steady-state
%   The function smooths the input y-data, and then finds the slope of the
%   smoothed data.  The minimum value of slope is taken to be the steady
%   state value.

smoothed = smoothdata(data_y,'gaussian',window);    %smooth the data
slope=slopefinder(data_x,smoothed,window);      %find slope
index=find(slope < threshold);              %locate index where slope drops below threshold
tss=data_x(index(1));                   %find x-value at that index
steadyvalue=data_y(index(1));           %find y-value at that index
end