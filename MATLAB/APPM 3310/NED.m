%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

function [N_local,E_local,D_local] = NED(invec)
%NED    This function constructs local NED axes for a position
%Inputs:    invec is the COLUMN VECTOR position in Earth Centered
%Outputs:   N_local: local North COLUMN VECTOR
%           E_local: local East COLUMN VECTOR
%           D_local: local Down COLUMN VECTOR

omega = atan2(invec(2),invec(1));  %longitude [rad]
a=6378137;                       %Earth dimensions [m]
b=6356752.3142;
                                 %latitude [rad]
alpha = atan((a^2/b^2)*(invec(3)/sqrt(invec(1)^2+invec(2)^2)));
N0=[0;0;1];                     %NED UNIT vectors @ (0,0,0)
E0=[0;1;0];
D0=[-1;0;0];

E_local=rotate(omega,N0,E0);        %local East vector
N_local=rotate(alpha,-E_local,N0)   %local North vector
D_local=cross(N_local,E_local)      %local Down vector
end