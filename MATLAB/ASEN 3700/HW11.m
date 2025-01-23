%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%% Question 1
clear;  const = getConstOrbitz();

T = 2*60*60;    n = (2*pi)/T;   t1 = 30*60;
R_0 = [0;6000;0];
V_0 = [0;-3;0];

[PHI_rr,PHI_rv,PHI_vr,PHI_vv] = CWmatrix(n,t1);

R_1 = PHI_rr*R_0 + PHI_rv*V_0;
V_1 = PHI_vr*R_0 + PHI_vv*V_0;
r1 = norm(R_1);     v1 = norm(V_1);

%% Question 2
clear;  const = getConstOrbitz();

r_orbit = 6600;     T = 2*pi*sqrt(r_orbit^3/const.Earth.mu);
n = (2*pi)/T;
R_0 = [1000;1000;1000];
v_0 = [0;0;5];

t1 = T/3;
[PHI_rr,PHI_rv,PHI_vr,PHI_vv] = CWmatrix(n,t1);

delta_v1 = -inv(PHI_rv)*PHI_rr*R_0;
delta_v2 = PHI_vr*R_0 +PHI_vv*delta_v1;
delta_v_total = norm(delta_v2) + norm(delta_v1-v_0);

%% Question 3
clear;  const = getConstOrbitz();

altA = 300; rA = 300+const.Earth.radius;
T = 2*pi*sqrt(rA^3/const.Earth.mu);         n = (2*pi)/T;
R0B = [-30;0;0];
V0B = [0;60*n;0];

t = linspace(0,T,1000);
X = -30*cos(n*t);
Y = 60*sin(n*t);
plot(X,Y);
set(gca, 'YDir','reverse');
grid on;

t1 = T/3;       Rf = [0;100;0];
[PHI_rr,PHI_rv,PHI_vr,PHI_vv] = CWmatrix(n,t1);
V1B = inv(PHI_rv)*(Rf - PHI_rr*R0B);
delta_v1 = V1B - V0B;
V2B = PHI_vr*R0B +PHI_vv*V1B;
delta_v2 = -V2B;