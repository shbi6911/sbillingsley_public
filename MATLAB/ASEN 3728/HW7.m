%By:        Shane Billingsley
%Class:     ASEN 3728 Aircraft Dynamics
%Date:      Spring 2024

%% Question 1
%geometric givens
W = 126000; %weight in pounds
S = 2000;   %planform area in ft^2
b = 120;    %wingspan in ft
c = 18.94;  %average chord in ft
Ix = 115000;    %slugs*ft^2
Iy = 2450000;   
Iz = 4070000;
%stability derivatives
C_y_beta = -0.877;  C_l_beta = -0.196;  C_n_beta = 0.139;
C_l_p = -0.381; C_n_p = -0.049; C_l_r = 0.198;  C_n_r = -0.185;
C_l_delta_a = -0.038;   C_n_delta_a = 0.017;
C_y_delta_r = 0.216;    C_l_delta_r = 0.0226;
%coefficients for Part 1
C_L = 0.68; C_D = 0.08;
%pull air density
[rho,~,~,~,~,~] = stdatmo(0,0,'US');
%thrust is equal to drag which is determined by velocity needed for L=W
T = W*(C_D/C_L);
%Assume one engine is at max thrust, determine remaining thrust needed
T_left = 6000; %remaining left engine thrust in pounds
T_right_1 = 6000; %right engines thrust in pounds (per engine)
T_right_2 = T - (T_left+T_right_1);
%find yawing moment
arm = b/6;  %moment arm for engine
%find total moment due to asymmetrical thrust
N_thrust = -(T_right_1*arm) - (T_right_2*2*arm) + (T_left*arm);
delta_r = deg2rad(3);   %rudder deflection
C_n_delta_r = (N_thrust*C_L)/(W*b*delta_r);%necessary stability derivative
%find aileron deflection and sideslip angle by solving a linear system
%involving the yaw moment (from previous) and the roll moment (assumed 0)
mat = [C_n_beta C_n_delta_a; C_l_beta C_l_delta_a];
soln = mat\[(N_thrust*C_L)/(W*b);0];
beta = rad2deg(soln(1));    aileron = rad2deg(soln(2));

%% Question 2
u_0 = 10;   %velocity in m/s
theta_0 = 0;    %initial theta angle
A_lat = [-0.2472 -0.0671 -9.7797 9.8100; -0.7966 -16.5375 1.8114 0;
        0.4607 -0.3451 -0.4586 0; 0 1 0 0];%given lateral matrix
[V2,D2] = eig(A_lat);   %find eigenvalues/vectors
values = diag(D2);      %extract eigenvalues
real_values = values(find(imag(values) == 0));  %extract real eigenvalues
%dutch roll eigenvalues are the complex conjugate pair
d_roll_values = values(find(imag(values) ~= 0));
%spiral mode value is the smaller absolute value
spiral_value = real_values(find(abs(real_values) == min(abs(real_values))));
%roll mode value is the larger absolute value
roll_value = real_values(find(abs(real_values) == max(abs(real_values))));
tau_roll = 1/abs(roll_value);   %time constant is 1/lambda
tau_spiral = 1/abs(spiral_value);
%natural frequency is sqrt(a^2 + b^2)
d_roll_omega_n = sqrt(real(d_roll_values(1))^2 + imag(d_roll_values(1))^2);
%damping coefficient is a/omega_n
d_roll_zeta = abs(real(d_roll_values(1)))/d_roll_omega_n;

%find tau for roll mode approximation
tau_roll_approx = 1/abs(A_lat(2,2));
%find percent error
roll_error = abs((tau_roll - tau_roll_approx)/tau_roll)*100;

%Dutch roll approximation
A_lat_d_roll = [A_lat(1,1) -u_0; A_lat(3,1) A_lat(3,3)];%revised matrix
d_roll_approx_values = eig(A_lat_d_roll);   %eigenvalues
d_roll_omega_n_approx = norm(d_roll_approx_values(1));       %natural frequency
%damping coefficient
d_roll_zeta_approx = abs(real(d_roll_approx_values(1)))/d_roll_omega_n_approx;
%percent errors
d_roll_omega_n_error = abs((d_roll_omega_n - d_roll_omega_n_approx)/d_roll_omega_n)*100;
d_roll_zeta_error = abs((d_roll_zeta - d_roll_zeta_approx)/d_roll_zeta)*100;

%roll and spiral mode approximation
Lv = A_lat(2,1);    Lp = A_lat(2,2);    Lr = A_lat(2,3);    g=A_lat(1,4);
Nv = A_lat(3,1);    Np = A_lat(3,2);    Nr = A_lat(3,3);

C = u_0*Nv; D= u_0*((Lv*Np) - (Lp*Nv)) - (g*Lv);
E = g*((Lv*Nr) - (Lr*Nv));
%solve characteristic eqn for lambdas
roll_spiral_values_approx = [(-D + sqrt(D^2 - (4*C*E)))/(2*C);
                             (-D - sqrt(D^2 - (4*C*E)))/(2*C)];
%separate roll and spiral
spiral_value_approx = roll_spiral_values_approx(find(abs(roll_spiral_values_approx) == min(abs(roll_spiral_values_approx))));
roll_value_approx = roll_spiral_values_approx(find(abs(roll_spiral_values_approx) == max(abs(roll_spiral_values_approx))));
%time constants
tau_roll_approx_2 = 1/abs(roll_value_approx);
tau_spiral_approx = 1/abs(spiral_value_approx);
%percent errors
roll_error_2 = abs((tau_roll - tau_roll_approx_2)/tau_roll)*100;
spiral_error = abs((tau_spiral - tau_spiral_approx)/tau_spiral)*100;

%% Question 3
A_dr = [-0.056 -730; 0.0012 -0.15]; B_dr = [0 5.5;0.004 -0.47];

d_roll_approx_values_3 = eig(A_dr);   %eigenvalues
d_roll_omega_n_approx_3 = sqrt(real(d_roll_approx_values_3(1))^2 + ...
    imag(d_roll_approx_values_3(1))^2);       %natural frequency
%damping coefficient
d_roll_zeta_approx_3 = abs(real(d_roll_approx_values_3(1)))/d_roll_omega_n_approx_3;
%percent errors

%find k for specified yaw damping
zeta = 0.7; omega_n = 0.9656;   lambda = -zeta*omega_n + 1i*omega_n*sqrt(1-zeta^2);
a = A_dr(1,1);  b= A_dr(1,2);   c= A_dr(2,1);   d=A_dr(2,2);
e = B_dr(1,2);  g= B_dr(2,2);   f= B_dr(2,1);
k = (lambda^2 - (lambda*(a+d))+((a*d)-(c*b)))/((a*g)-(g*lambda)-(c*e));

%find transfer function
num = [0 e -(d*e + b*g)];
den = [1 -(a+d) -b*c];
disp(num);
disp(den);

%find maximum yaw rate
delta_r_max = deg2rad(25);
r_max = abs(delta_r_max/k);