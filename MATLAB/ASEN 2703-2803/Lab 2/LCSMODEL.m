%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%{
Calculates the vertical speed of the sleeve

Inputs :
    r {float} - radius of the rotating disk [mm]
    d {float} - distance between the center of the disk and the wall [mm]
    l {float} - length of rod connecting disk to collar [mm]
    theta {float} - angle of the disk relative to the vertical [deg]
    w {float} - angular velocity of the disk (derivative of theta) [deg/s]

Outputs:
    V_slide {float} - vertical speed of the collar [mm/s]
%}
function [V_slide] = LCSMODEL(r, d, l, theta, w)

% convert theta and w to radians and rads/s repsectively
theta = deg2rad(theta);
w = deg2rad(w);
                                     
beta = asin((d - (r .* sin(theta))) ./ l); % [ rad] angle between rod and collar

V_slide = -(w .* r .* sin(theta)) - (w .* r .* cos(theta) .* tan(beta));
end