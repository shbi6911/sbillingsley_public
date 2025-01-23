%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

load ("project2verification.mat");

figure(1);
hold on; grid on;
plot(verification.time,verification.height,"LineWidth",2);       % plot z-position
title ("Bottle Rocket Trajectory (Z)");
xlabel ("Time (s)");
ylabel ("X-Position (m)");
hold off;