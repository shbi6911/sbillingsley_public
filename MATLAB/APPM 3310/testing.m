%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

axis = [0;0;1];
vector = [1;1;1];
theta = pi/4;
newvector = rotate(theta, axis, vector);

figure (1); hold on; grid on;
quiver3(0,0,0,1,1,1);
quiver3(0,0,0,newvector(1),newvector(2),newvector(3));