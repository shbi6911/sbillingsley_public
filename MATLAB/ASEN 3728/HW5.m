%By:        Shane Billingsley
%Class:     ASEN 3728 Aircraft Dynamics
%Date:      Spring 2024

%% Question 1
close; clear all;
%givens
A_lon = [-0.020 0.016 -0.65 -32.17; -0.13 -1.019 454.21 0;...
     0 -0.0050 -1.38 0; 0 0 1.0 0]; %longitudinal matrix
V_a = 502;  %trim velocity in ft/s
[V1,D1] = eig(A_lon); %get eigenvalues/eigenvectors
values = diag(D1);
%find indices for short period and phugoid poles
sp_val_ind = find(real(values) == min(real(values)));
ph_val_ind = find(real(values) == max(real(values)));
%sort values and vectors into short period and phugoid
sp_values = values(sp_val_ind); ph_values = values(ph_val_ind);
sp_vectors = V1(:,sp_val_ind);  ph_vectors = V1(:,ph_val_ind);
%find zeta and omega_n for short period and phugoid
sp_omega_n = sqrt(real(sp_values(1))^2 + imag(sp_values(1))^2);
sp_zeta = abs(real(sp_values(1))/sp_omega_n);
ph_omega_n = sqrt(real(ph_values(1))^2 + imag(ph_values(1))^2);
ph_zeta = abs(real(ph_values(1))/ph_omega_n);
%output results for reference
%disp("Short Period Frequency is " + string(sp_omega_n)...
%    + " and Damping Coefficient is "+ string(sp_zeta));
%disp("Phugoid Frequency is " + string(ph_omega_n)...
%    + " and Damping Coefficient is "+ string(ph_zeta));

%normalize vectors
sp_vectors_norm = sp_vectors;   ph_vectors_norm = ph_vectors;
sp_vectors_norm = sp_vectors_norm./sp_vectors_norm(4,:);
ph_vectors_norm = ph_vectors_norm./ph_vectors_norm(4,:);
%display for reference
%disp("Short Period")
%disp(sp_vectors_norm);
%disp("Phugoid");
%disp(ph_vectors_norm);

%generate phasor plots

%nondimensionalize eigenvectors for phasor plot
sp_vectors_norm(1:2,:) = sp_vectors_norm(1:2,:)./V_a;
ph_vectors_norm(1:2,:) = ph_vectors_norm(1:2,:)./V_a;

% figure();   hold on;    grid on;
% plot([0 real(sp_vectors_norm(4,1))],[0 imag(sp_vectors_norm(4,1))], 'LineWidth',1.5,'Color','b');
% plot([0 real(sp_vectors_norm(3,1))],[0 imag(sp_vectors_norm(3,1))], 'LineWidth',1.5,'Color','g');
% plot([0 real(sp_vectors_norm(2,1))],[0 imag(sp_vectors_norm(2,1))], 'LineWidth',1.5,'Color','y');
% plot([0 real(sp_vectors_norm(1,1))],[0 imag(sp_vectors_norm(1,1))], 'LineWidth',1.5,'Color','r');
% legend("\delta \theta","\delta q","\delta w/u_0","\delta u/u_0");
% title("Phasor Plot for Short Period Mode");
% xlabel("Real");
% ylabel("Imaginary");
% hold off;
% 
% figure();   hold on;    grid on;
% plot([0 real(ph_vectors_norm(4,1))],[0 imag(ph_vectors_norm(4,1))], 'LineWidth',1.5,'Color','b');
% plot([0 real(ph_vectors_norm(3,1))],[0 imag(ph_vectors_norm(3,1))], 'LineWidth',1.5,'Color','g');
% plot([0 real(ph_vectors_norm(2,1))],[0 imag(ph_vectors_norm(2,1))], 'LineWidth',1.5,'Color','y');
% plot([0 real(ph_vectors_norm(1,1))],[0 imag(ph_vectors_norm(1,1))], 'LineWidth',1.5,'Color','r');
% legend("\delta \theta","\delta q","\delta w/u_0","\delta u/u_0");
% title("Phasor Plot for Phugoid Mode");
% xlabel("Real");
% ylabel("Imaginary");
% hold off;

%% Question 2
%givens
B = [-0.244;-1.46;-0.2;0];  %controls matrix
k1 = -5.0;  k2 = -0.005;    %gain values

A_sp = A_lon(2:3,2:3); %slice off short period approx from full A matrix
%get eigenvalues/eigenvectors
[V2,D2]=eig(A_sp);
values2 = diag(D2);
%find zeta and omega_n for short period
sp_approx_omega_n = sqrt(real(values2(1))^2 + imag(values2(1))^2);
sp_approx_zeta = abs(real(values2(1))/sp_approx_omega_n);
%disp("Short Period Approximation Frequency is " + string(sp_approx_omega_n)...
    %+ " and Damping Coefficient is "+ string(sp_approx_zeta));
%find closed loop matrix for approximation
K = [-k2/V_a -k1];
A_cl_approx = A_sp + B(2:3)*K;
%disp(A_cl_approx);
values2_cl_approx = eig(A_cl_approx);
%disp(values2_cl);

%find closed loop matrix for full longitudinal
BK = B*[0 K 0];
%disp(BK);
A_cl = A_lon + BK;
values2_cl = eig(A_cl);
%disp(A_cl);
%disp(values2_cl);

%% Question 3
%givens
A3 = [-0.0025 -30; 0.0001 0];
B3 = [10;0];    k1 = linspace(0.0001,0.01,100);    k2 = zeros(1,length(k1));

%get open loop eigenvalues
values3 = eig(A3);

%find natural freq and damping from eigenvalues
omega_n_3 = sqrt(real(values3(1))^2 + imag(values3(1))^2);
zeta_3 = abs(real(values3(1))/omega_n_3);
%disp("Phugoid Approximation Frequency is " + string(omega_n_3)...
%    + " and Damping Coefficient is "+ string(zeta_3));
%loop through possible gains and find damping coefficient
zeta_3_vec = zeros(1,length(k2));
for i = 1:length(k2)
    %define closed loop matrix
    A3_cl = A3 - [10*k1(i) 10*k2(i); 0 0];
    values3_cl = eig(A3_cl);
    omega_n_3_cl = sqrt(real(values3_cl(1))^2 + imag(values3_cl(1))^2);
    zeta_3_cl = abs(real(values3_cl(1))/omega_n_3_cl);
    zeta_3_vec(i) = zeta_3_cl;
end
figure();
plot(k1,zeta_3_vec);