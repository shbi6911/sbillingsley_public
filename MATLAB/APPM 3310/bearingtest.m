%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

R_ac=[-3.92;3.47;-3.63]*10^6;
R_bc=[4.03;0.30;4.92]*10^6;
R_ba=R_bc-R_ac;
% omega = atan2(R_ac(2),R_ac(1));
% a=6378137;
% b=6356752.3142;
% alpha = atan((a^2/b^2)*(R_ac(3)/sqrt(R_ac(1)^2+R_ac(2)^2)));
% N0=[0;0;1];
% E0=[0;1;0];
% D0=[-1;0;0];
% 
% E_local=rotate(omega,N0,E0)
% N_local=rotate(alpha,-E_local,N0)
% D_local=cross(N_local,E_local)

[N_local,E_local,D_local]=NED(R_ac);