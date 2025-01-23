%By:        Shane Billingsley
%Class:     ASEN 3713 Thermodynamics and Heat Transfer
%Date:      Fall 2023


C_1=3.74177*10^8;
C_2=1.43878*10^4;

T=5780;
lambda = logspace(-2,3,1000);
e = exp(C_2./(lambda.*T));
denominator = ((lambda).^5).*(e-1);
emisivity = C_1./denominator;

loglog(lambda,emisivity);
title("Emissivity of a Blackbody at 5780 K (Solar Temp)");
ylabel("Blackbody Emissive Power");
ylim([10^-6,10^9]);
xlabel("Wavelength");

print('ASEN3713HW9','-dpng','-r300');