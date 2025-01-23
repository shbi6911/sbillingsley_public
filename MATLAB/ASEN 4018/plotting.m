%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

%plot some data
data = readmatrix('Deorbit47kg');
figure();
plot(data(:,1)./365,data(:,2));
title("Deorbit of 47kg CubeSat with $2m^2$ drag area",'interpreter','latex');
xlabel("Elapsed Time (years)");
ylabel("Altitude (km)");