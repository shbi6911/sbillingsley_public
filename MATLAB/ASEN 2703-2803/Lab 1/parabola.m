%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

%define parabola
h0 = 125;
g= 9.81;
h = 0;
initialv = sqrt(2*g*(h0-h));
launchangle = pi/4;
a = -g/(2*initialv^2*(cos(launchangle))^2);
b = tan(launchangle);
c = h;
quad = roots([a b 0]);
finalx = quad(2);
X = linspace(0,finalx,1000);
Z = (a*X.^2)+X*b+c;
Y = zeros(1,1000);

%plot
figure (1);
plot3(X,Y,Z);
grid on; hold on;
X_check = -50:100;
Y_check = -50:100;
Z_check1 = ones(151,151)+44;
surf(X_check,Y_check,Z_check1);
hold off;

%calculate Gs
Zprime = b+(2*a*X);
Z2prime = 2*a;
rho_parabola = ((1+Zprime.^2).^(3/2))./abs(Z2prime);
phi_parabola = atan(Zprime);
G1 = g*rho_parabola;
G2 = 2*g*(h0-Z);
G_parabola = cos(phi_parabola)-(G2./G1);

figure(2);
plot(X,G_parabola);
ylim([-1 6]);
