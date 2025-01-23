%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

%% Question 1
clear;  const = getConstOrbitz();
r1 = 300 + const.Earth.radius;
r2 = 3000 + const.Earth.radius;

a_t = (r1 + r2)/2;
v_init = sqrt(const.Earth.mu/r1);
v_fin = sqrt(const.Earth.mu/r2);

energy_t = -(const.Earth.mu/(2*a_t));
v_p_t = sqrt(2*(energy_t + (const.Earth.mu/r1)));
v_a_t = sqrt(2*(energy_t + (const.Earth.mu/r2)));

delta_v_1 = v_p_t - v_init;
delta_v_2 = v_fin - v_a_t;

delta_v = delta_v_1 + delta_v_2;
TOF = 0.5*((2*pi)/sqrt(const.Earth.mu))*a_t^(3/2);

%% Question 2
clear;  const = getConstOrbitz();
r1 = 300 + const.Earth.radius;  %given radii
r2 = 3000 + const.Earth.radius;
v_init = sqrt(const.Earth.mu/r1);   %velocity of initial and final orbits
v_fin = sqrt(const.Earth.mu/r2);
e_1 = 0.3;                          %given eccentricity
a_t1 = r1/(1-e_1);          %semi-major axis of transfer orbit 1
r_a_t1 = (2*a_t1) - r1;     %radius of transfer orbit 1 at apogee
energy_t1 = -(const.Earth.mu/(2*a_t1)); %specific energy of transfer orbit 1
%velocity of transfer 1 at perigee
v_p_t1 = sqrt(2*(energy_t1 + (const.Earth.mu/r1)));
%velocity of transfer 1 at apogee
v_a_t1 = sqrt(2*(energy_t1 + (const.Earth.mu/r_a_t1)));
delta_v_1 = v_p_t1 - v_init;    %delta-v from orbit 1 to transfer orbit 1
a_t2 = (r_a_t1 + r2)/2;     %semi-major axis of transfer orbit 2
energy_t2 = -(const.Earth.mu/(2*a_t2)); %specific energy of transfer orbit 2
%velocity of transfer 2 at apogee
v_a_t2 = sqrt(2*(energy_t2 + (const.Earth.mu/r_a_t1)));
delta_v_2 = v_a_t2 - v_a_t1;    %delta-v from transfer 1 to transfer 2
%velocity of transfer 2 at perigee
v_p_t2 = sqrt(2*(energy_t2 + (const.Earth.mu/r2)));
delta_v_3 = v_p_t2 - v_fin;     %delta-v from transfer 2 to orbit 2
delta_v_transfer = delta_v_1 + delta_v_2 + delta_v_3;   %total delta-v
%time of flight for transfer orbit 1
TOF1 = 0.5*((2*pi)/sqrt(const.Earth.mu))*a_t1^(3/2);
%time of flight for transfer orbit 2
TOF2 = 0.5*((2*pi)/sqrt(const.Earth.mu))*a_t2^(3/2);
TOF = TOF1 + TOF2;  %total time of flight

%% Question 4
clear;  const = getConstOrbitz();
r1 = 15000;     %radius of circular orbit 1
gamma1 = 0;     %flight path angle of circular orbit is zero
rp2 = 500+const.Earth.radius;     %perigee radius of orbit 2
ra2 = 22000;    %apogee radius of orbit 2
a2 = (ra2 + rp2)/2;     %semi-major axis of orbit 2
e2 = (ra2 - rp2)/(ra2 + rp2);   %eccentricity of orbit 2
v1 = sqrt(const.Earth.mu/r1);   %velocity of a circular orbit
h1 = r1*v1;     %all velocity in circular orbit is perpendicular
h2 = sqrt(const.Earth.mu*a2*(1-e2^2));  %angular momentum of orbit 2
v2_perp = h2/r1;    %perpendicular component of orbit 2 velocity
theta = acos((((v2_perp*h2)/const.Earth.mu)-1)/e2);  %true anomaly
v2_r = (const.Earth.mu/h2)*e2*sin(theta);   %radial component of orbit 2
delta_gamma = atan2(v2_r,v2_perp);     %change in flight path angle
v2 = sqrt(v2_r^2 + v2_perp^2);    %speed of orbit 2
delta_v = sqrt(v1^2 + v2^2 - (2*v1*v2*cos(delta_gamma)));   %law of cosines

%% Question 5
clear;  const = getConstOrbitz();
e = 0.3;    %eccentricity
h = 60000;  %specific angular momentum
delta_i = deg2rad(90);  %change in inclination
r_a = (h^2/const.Earth.mu)*(1/(1-e));   %radius at apogee
v_a = h/r_a;    %speed at apogee
delta_v = 2*v_a*sin(delta_i/2);     %required delta-v

%% Question 6
clear;  const = getConstOrbitz();
v1_act = [4;11;0];                  %accidental delta-v
v_e = sqrt(const.Sun.mu/const.Earth.a); %velocity in Earth orbit
v1_f = v1_act + [0;v_e;0];          %velocity vector after accidental burn
%find specific energy of transfer orbit
energy1 = ((norm(v1_f)^2/2) - (const.Sun.mu/const.Earth.a));
a_t = -(const.Sun.mu/(2*energy1));  %semi-major axis of transfer orbit
h_t = const.Earth.a*v1_f(2);        %angular momentum of transfer orbit
e_t = sqrt(-(1/((a_t*const.Sun.mu)/h_t^2)-1));  %eccentricity of xfer orbit
theta = acos((((1-e_t^2)/(const.Earth.a/a_t))-1)/e_t);  %true anomoly
r_a_t = (h_t^2/const.Sun.mu)*(1/(1-e_t));   %radius of xfer orbit at apogee

