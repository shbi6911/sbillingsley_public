%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 2
clear; clc; const = getConst();
%givens
a = 10000;
n = sqrt(const.mu.Earth/a^3);
P = 2*pi*(1/n);
start_pos = [2;2;0];    start_vel = [-0.03;0.01;0.05];

[Phi_rr,Phi_rv,Phi_vr,Phi_vv] = CWmatrix(n,P/2);

one_plus = inv(Phi_rv)*([-2;2;0] - Phi_rr*start_pos);
delta_v1 = one_plus - start_vel;

two_minus = Phi_vr*start_pos + Phi_vv*one_plus;

delta_v2 = [0;0;0] - two_minus;

delta_v2_tot = norm(delta_v1) + norm(delta_v2);

%% plotting

t = linspace(0,P/2,1000);
x_pos = zeros(1,length(t));
y_pos = zeros(1,length(t));
for i = 1:length(t)
    [PHI_rr,PHI_rv,~,~] = CWmatrix(n,t(i));
    current_pos = PHI_rr*start_pos + PHI_rv*one_plus;
    x_pos(i) = current_pos(1);
    y_pos(i) = current_pos(2);
end

start_vel_arrow = start_vel(1:2)./norm(start_vel(1:2));
maneuver1_vel_arrow = one_plus(1:2)./norm(one_plus(1:2));
end_vel_arrow = two_minus(1:2)./norm(two_minus(1:2));

figure(); hold on; grid on;
plot(y_pos,x_pos);
quiver(y_pos(1),x_pos(1),start_vel_arrow(2),start_vel_arrow(1),'MaxHeadSize',1);
quiver(y_pos(1),x_pos(1),maneuver1_vel_arrow(2),maneuver1_vel_arrow(1),'MaxHeadSize',1);
quiver(y_pos(end),x_pos(end),end_vel_arrow(2),end_vel_arrow(1),'MaxHeadSize',1);
plot(y_pos(1),x_pos(1),'.','MarkerSize',15);
plot(y_pos(end),x_pos(end),'.','MarkerSize',15);
plot(0,0,'.','MarkerSize',20);
arrow([y_pos(250),x_pos(250)],[y_pos(255),x_pos(255)]);
arrow([y_pos(750),x_pos(750)],[y_pos(755),x_pos(755)]);
xlim([-3 3]);
ylim([-3,3]);
set(gca, 'XDir','reverse')
title("CubeSat Trajectory Relative to Primary Satellite");
xlabel("Radial Position (m)");
ylabel("Along-Track Position (m)");
legend("Trajectory","Initial Velocity","Velocity after Maneuver 1",...
    "Velocity before Maneuver 2","Initial Position","Final Position","Primary Spacecraft",'','');
axis equal


%axis equal

