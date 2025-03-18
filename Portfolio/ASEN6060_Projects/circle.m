function h = circle(x,y,r,color,name)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025
%
% circle is a simple function to plot a circle of defined center and radius
%
%INPUTS:    x,y     coordinates of the circle center
%           r       radius of the circle
%           color   color value of the plotted circular line
%           name    string or char array of the DisplayName of the line obj
%
%OUTPUTS    h       a Line object with properties as above
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'LineWidth',2);
    h.Color = color;
    h.DisplayName = name;
    hold off
end