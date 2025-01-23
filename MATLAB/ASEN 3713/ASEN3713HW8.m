%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023

b=150/(8000*0.01*570);

v=linspace(0.005,0.06,1000);

T=(exp(-b.*(3./v)))*(18-950)+950;

figure(1);
hold on;
title("Temperature vs. Velocity");
xlabel("Velocity m/s");
ylabel("Temperature (deg C)");
plot(v,T,'r');
print("tempplot",'-dpng','-r300');