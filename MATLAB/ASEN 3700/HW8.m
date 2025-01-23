%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%% Question 1
clear;
const = getConstOrbitz();
v = 2.23; %km/s
r = 402000; %km
theta = -150; %degrees
energy = ((v^2)/2) - (const.Earth.mu/r);
p1 = const.Earth.mu^2/(2*energy);
p2 = -const.Earth.mu*r*cosd(theta);
p3 = -(p1 + (const.Earth.mu*r));
e_vec = roots([p1 p2 p3]);
e = e_vec(2);
h = sqrt(const.Earth.mu*r*(1+(e*cosd(theta))));
rp = (h^2/const.Earth.mu)*(1/(1+e));
vp = h/rp;
%disp(e);
%disp(rp);
%disp(vp);

%% Question 2
clear;
const = getConstOrbitz();
r_p = 200 + const.Earth.radius;
r_a = 600 + const.Earth.radius;
r = 400 +const.Earth.radius;
e = (r_a - r_p)/(r_a + r_p);
a = 0.5*(r_a + r_p);
thingy = sqrt((1-e)/(1+e));
T = ((2*pi)/sqrt(const.Earth.mu))*a^(3/2);
theta = acos((a*(1-e^2)-r)/(r*e));
E = 2*atan(thingy*tan(theta/2));
M = E - e*sin(E);
t = (M*T)/(2*pi);
t_tot = T-2*t; t_tot_min = t_tot/60;
%disp(t_tot_min);

%% Question 3
clear;
const = getConstOrbitz();
T = 14*3600; %s
r_p = 10000; %km
t = 10*3600; %s
a = ((sqrt(const.Earth.mu)*T)/(2*pi))^(2/3);
r_a = 2*a - r_p;
e = (r_a - r_p)/(r_a + r_p);
M = (2*pi*t)/T;
E = kepler_E(e,M);
thingy = sqrt((1-e)/(1+e));
theta = 2*atan((1/thingy)*tan(E/2));
r = (a*(1-e^2))/(1+(e*cos(theta)));
v = sqrt(2*((const.Earth.mu/r) - (const.Earth.mu/(2*a))));
h = sqrt(const.Earth.mu*a*(1-e^2));
v_r = (const.Earth.mu/h)*e*sin(theta);
%disp(r);
%disp(v);
%disp(v_r);

%% Question 4
clear;
const = getConstOrbitz();
r_a = 16000; %km
r_p = 7500; %km
time = 40*60; %s
theta_0 = deg2rad(80); %rad
e = (r_a - r_p)/(r_a + r_p);
a = 0.5*(r_a + r_p);
T = ((2*pi)/sqrt(const.Earth.mu))*a^(3/2);
thingy = sqrt((1-e)/(1+e));
E = 2*atan(thingy*tan(theta_0/2));
M = E - e*sin(E);
t = (M*T)/(2*pi);
t2 = t+time;
M2 = (2*pi*t2)/T;
E2 = kepler_E(e,M2);
theta2 = 2*atan((1/thingy)*tan(E2/2));
disp(rad2deg(theta2));

%% Question 5
clear;
const = getConstOrbitz();
h = 1.403888e9; %km^2/s
T = 92.057251*3600*24; %s
theta1 = deg2rad(137.961581); %rad
a = ((sqrt(const.Sun.mu)*T)/(2*pi))^(2/3);
e = sqrt(1-(h^2/(const.Sun.mu*a)));
thingy = sqrt((1-e)/(1+e));
E = 2*atan(thingy*tan(theta1/2));
M = E - e*sin(E);
t = (M*T)/(2*pi);
t_days = t/(3600*24);
r = a*((1-e^2)/(1+(e*cos(theta1))));
rp = a*(1-e);
ra = 2*a - rp;
t_a = T/2;
t_remain = t_a - t;
t2 = t + (50*3600*24);
M2 = (2*pi*t2)/T;
E2 = kepler_E(e,M2);
theta2 = 2*atan((1/thingy)*tan(E2/2));
r2 = a*((1-e^2)/(1+(e*cos(theta2))));
v_theta = h/r2;
v_r = (const.Sun.mu/h)*e*sin(theta2);
b = a*sqrt(1-e^2);
p = a*(1-e^2);
%plotting
t_vec = linspace(0,(2*pi),500);
X = a*cos(t_vec);
Y = b*sin(t_vec);
plot(X,Y); hold on; grid on;
Y2 = a*sin(t_vec);
plot(X,Y2,'r');
plot(0,0,'.','MarkerSize',10,'Color','k');
plot(a*e,0,'.','MarkerSize',10,'Color','k');
X3 = linspace(-a,a,1000);
Y3 = X3*tan(E);
Y4 = (X3-(a*e))*tan(theta1);
plot(X3,Y3,'k');
plot(X3,Y4,'k');
xlim([-(1.1*a),(1.1*a)]);
ylim([-(1.1*a),(1.1*a)]);
axis equal


%% functions
function E = kepler_E(e, M)
error = 1e-8;
if M < pi
    E = M+ e/2;
else
    E = M-e/2;
end
ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1-e*cos(E));
    E = E-ratio;
end
end

