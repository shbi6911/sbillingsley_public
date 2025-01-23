%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

function verify()

load ("project2verification.mat");

figure(1);
hold on; grid on;
plot(verification.time,verification.height,"LineWidth",2);       % plot z-position
title ("Bottle Rocket Trajectory (Z)");
xlabel ("Time (s)");
ylabel ("X-Position (m)");
hold off;

figure(2);
hold on; grid on;
plot(verification.time,verification.thrust,"LineWidth",2);       % plot thrust
title ("Bottle Rocket Thrust");
xlabel ("Time (s)");
ylabel ("Thrust (N)");
hold off;

figure(3);
hold on; grid on;
plot(verification.time,verification.distance,"LineWidth",2);       % plot x-position
title ("Bottle Rocket Trajectory (X)");
xlabel ("Time (s)");
ylabel ("Z-Position (m)");
hold off;

figure(4);
hold on; grid on;
plot(verification.time,verification.xVelocity,"LineWidth",2);       % plot x-velocity
title ("Bottle Rocket Velocity (X)");
xlabel ("Time (s)");
ylabel ("X-Velocity (m/s)");
hold off;

figure(5);
hold on; grid on;
plot(verification.time,verification.zVelocity,"LineWidth",2);       % plot z-velocity
title ("Bottle Rocket Velocity (Z)");
xlabel ("Time (s)");
ylabel ("Z-Velocity (m/s)");
hold off;

figure(6);
hold on; grid on;
plot(verification.time,verification.airVolume,"LineWidth",2);       % plot air volume
title ("Bottle Rocket Air Volume");
xlabel ("Time (s)");
ylabel ("Air Volume (m^3)");
hold off;

end