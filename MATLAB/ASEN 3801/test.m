close; clear all;
ap = ttwistor();
rho = stdatmo(1800);
v_vec = [20.99; 0; 0.5837];
u_0 = v_vec(1);

C_zalpha = -ap.CD0-ap.CLalpha;  Zw = 0.5*rho*u_0*ap.S*C_zalpha;

Mw = 0.5*rho*u_0*ap.c*ap.S*ap.Cmalpha;
Mw_dot = 0.25*rho*ap.c^2*ap.S*ap.Cmalphadot;
Mq = 0.25*rho*u_0*ap.c^2*ap.S*ap.Cmq;

A_lon_sp = [Zw/ap.m,u_0;...
    (1/ap.Iy)*(Mw + (Mw_dot*Zw)/ap.m),(1/ap.Iy)*(Mq + Mw_dot*u_0)];

Xu = -rho*u_0*ap.S*ap.CD0; 
C_weight = ap.W/(0.5*rho*u_0^2*ap.S);
Zu = -rho*u_0*ap.S*C_weight;

A_lon_ph = [Xu/ap.m, -ap.g; -Zu/(ap.m*v_vec(1)), 0];

sp_values = eig(A_lon_sp);
ph_values = eig(A_lon_ph);
% [V1,D1] = eig(A_lon); %get eigenvalues/eigenvectors
% values = diag(D1);
%find indices for short period and phugoid poles
% sp_val_ind = find(real(values) == min(real(values)));
% ph_val_ind = find(real(values) == max(real(values)));
% %sort values and vectors into short period and phugoid
% sp_values = values(sp_val_ind); ph_values = values(ph_val_ind);
% sp_vectors = V1(:,sp_val_ind);  ph_vectors = V1(:,ph_val_ind);
%find zeta and omega_n for short period and phugoid
sp_omega_n = sqrt(real(sp_values(1))^2 + imag(sp_values(1))^2);
sp_zeta = abs(real(sp_values(1))/sp_omega_n);
ph_omega_n = sqrt(real(ph_values(1))^2 + imag(ph_values(1))^2);
ph_zeta = abs(real(ph_values(1))/ph_omega_n);
%output results for reference
disp("Short Period Frequency is " + string(sp_omega_n)...
   + " and Damping Coefficient is "+ string(sp_zeta));
disp("Phugoid Frequency is " + string(ph_omega_n)...
   + " and Damping Coefficient is "+ string(ph_zeta));