%By:        Shane Billingsley
%Class:     ASEN 3728 Aircraft Dynamics
%Date:      Spring 2024

clear; close all;

d = 0.1; d_m = d/sqrt(2);   k_m = 0.003;
mom = [-5;0;0.2;0.01];

mat = [-1 -1 -1 -1; -d_m -d_m d_m d_m;...
    d_m -d_m -d_m d_m; k_m -k_m k_m -k_m];

%force = mat\mom;
force = inv(mat)*mom;
%disp(force);

A = [0 1; -3.24 -1.8];
[V,D] = eig(A);
disp(D);
disp(V);

sys= ss(A,eye(2),eye(2),0);
step(sys);
