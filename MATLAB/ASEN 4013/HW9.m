%By:    Shane Billingsley
%Class: ASEN 4013 Foundations of Propulsion
%Date:  Fall 2024

load('CEARUNdataHW9');

pip = linspace(1,10,16)./0.101;     g0 = 9.80665;
%parse out data
rho_chamber = CEARUNdataHW9.rho(1);
rho_throat = CEARUNdataHW9.rho(2);
rho = CEARUNdataHW9.rho(3:18);

temp_chamber = CEARUNdataHW9.t(1);
temp_throat = CEARUNdataHW9.t(2);
temp = CEARUNdataHW9.t(3:18);

isp_chamber = CEARUNdataHW9.isp(1)/g0;
isp_throat = CEARUNdataHW9.isp(2)/g0;
isp = CEARUNdataHW9.isp(3:18)./g0;

cf_chamber = CEARUNdataHW9.cf(1);
cf_throat = CEARUNdataHW9.cf(2);
cf = CEARUNdataHW9.cf(3:18);

mw_chamber = CEARUNdataHW9.mw(1);
mw_throat = CEARUNdataHW9.mw(2);
mw = CEARUNdataHW9.mw(3:18);

gam_chamber = CEARUNdataHW9.gam(1);
gam_throat = CEARUNdataHW9.gam(2);
gam = CEARUNdataHW9.gam(3:18);

cstar = (isp.*g0)./cf;

%plotting
figure();
plot(pip,rho);
title("Density");
xlabel("$\frac{P_c}{P_e}$ for $P_e$ = 101 kPa",'interpreter','latex');
ylabel("$\frac{kg}{m^3}$",'interpreter','latex','fontsize',20);

figure();
plot(pip,temp);
title("Temperature");
xlabel("$\frac{P_c}{P_e}$ for $P_e$ = 101 kPa",'interpreter','latex');
ylabel("Kelvin",'interpreter','latex');

figure();
plot(pip,isp);
title("Specific Impulse");
xlabel("$\frac{P_c}{P_e}$ for $P_e$ = 101 kPa",'interpreter','latex');
ylabel("seconds",'interpreter','latex');

figure();
plot(pip,cstar);
title("Characteristic Velocity");
ylim([1300,1400]);
xlabel("$\frac{P_c}{P_e}$ for $P_e$ = 101 kPa",'interpreter','latex');
ylabel("$\frac{m}{s}$",'interpreter','latex','fontsize',20);

ratio = 1000/14.7;
rho_comp = interp1(pip,rho,ratio);
temp_comp = interp1(pip,temp,ratio);
isp_comp = interp1(pip,isp,ratio);
mw_comp = interp1(pip,mw,ratio);
gamma_comp = interp1(pip,gam,ratio);