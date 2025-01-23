%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%% Parameters


K_g = 33.3; % Gear ratio (unitless)
K_m = 0.0401; % Motor constant (V/rad/s)
J = 0.0005 + (0.2 * .2794^2) + 0.0015;   % Load Inertia (kg/m^2)
R_m = 19.2; % Motor resistance (ohms)

%K_p = input("Proportional Gain ");  % Proportional gain
%K_d = input("Derivative Gain ");    % Derivative gain
K_p = 10.5;
K_d = 0.351;
%% Closed Loop System

num = (K_p * K_g * K_m)/(J * R_m);

% Coefficient of s^2
d2 = 1;
% Coefficient of s
d1 = ((K_g * K_m * K_d)/(J * R_m)) + (((K_g * K_m)^2)/(J * R_m));
% Constant term
d0 = num;

den = [d2 d1 d0];
sysTF = tf(num,den);


%% Step Response
opts = stepDataOptions('StepAmplitude',0.5);
[x,t] = step(sysTF,opts);

plot(t,x);
yline(0.5+(0.5*0.2));
yline(0.5-(0.5*0.2));
yline(0.5+(0.5*0.05));
yline(0.5-(0.5*0.05));
xline(1);
xlabel("Time (s)");
ylabel("Angle (rad)");
hold on;