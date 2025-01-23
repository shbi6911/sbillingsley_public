%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

function rot = WindRot(phi,lambda)
%WindRot constructs a rotation matrix to go from Earth-fixed axes to a
%ground-fixed frame 
%   INPUTS       latitude and longitude expressed in RADIANS
%
C1 = [0,0,1;0,1,0;-1,0,0]; %rotate 90 deg about y-axis
%rotate about intermediate x-axis by longitude
C2 = [1,0,0;0,cos(lambda),-sin(lambda);0,sin(lambda),cos(lambda)];
%rotate about intermediate y-axis by latitude
C3 = [cos(-phi),0,sin(-phi);0,1,0;-sin(-phi),0,cos(-phi)];
rot = C3*C2*C1;
end