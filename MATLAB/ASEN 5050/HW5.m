%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc; const = getConst();
%givens
a1 = 7045;      e1 = 0.23;  theta_star1 = deg2rad(-142);
mi = 1224;  %initial mass in kg
Isp = 212;  %specific impulse
g0 = 9.81;   %acceleration of gravity in m/s
%find specific angular momentum
p1 = a1*(1-e1^2);
h1 = sqrt(p1*const.mu.Moon);
%find flight path angle using h, r, and v
r1 = p1/(1+e1*cos(theta_star1));
v1 = sqrt(const.mu.Moon*((2/r1)-(1/a1)));
gamma = -acos(h1/(r1*v1));
%find components of velocity and express as a vector
v_r1 = v1*sin(gamma);   v_theta1 = v1*cos(gamma);
V1_roh = [v_r1;v_theta1;0];
%find vector v2 after maneuver given delta-v
delta_V1 = [0.30;-0.10;0];
V2_roh = V1_roh + delta_V1;
%calculate a and e of new orbit
v2 = norm(V2_roh);                      %velocity
a2 = const.mu.Moon/((2*const.mu.Moon/r1) - v2^2);   %semimajor axis
H2_roh = cross([r1;0;0],V2_roh);    h2 = norm(H2_roh);  %angular mom
e2 = sqrt(1- (h2^2)/(const.mu.Moon*a2));    %eccentricity

p2 = a2*(1-e2^2);       %semi-latus rectum
theta_star2 = acos((p2 - r1)/(r1*e2));  %true anomaly

%find change in argument of periapsis
evec1 = cross(V1_roh, cross([r1;0;0],V1_roh))./const.mu.Moon - [1;0;0];
evec2 = cross(V2_roh, H2_roh)./const.mu.Moon - [1;0;0];
delta_omega = acos(dot(evec1,evec2)/(e1*e2));
if norm(cross(evec1,evec2)) < 0
    delta_omega = -delta_omega;
end

%find change in argument of periapsis using different method
delta_omega_check = theta_star1 - theta_star2;

%find propellant used
exponent = (-norm(delta_V1)*1000)/(Isp*g0);
prop_mass = mi - (mi*exp(exponent));

%% Problem 2
clear; clc; const = getConst();
%givens
r1 = 6500;      E1 = pi/2;      rp1 = 5915;
rp2 = 5712;     ra2 = 7888;
%find semimajor axis and eccentricity of orbit 1
a1 = r1;        e1 = 1 - (rp1/a1);

%find specific angular momentum of orbit 1
p1 = a1*(1-e1^2);
h1 = sqrt(p1*const.mu.Mars);
%find flight path angle 1 using h, r, and v
v1 = sqrt(const.mu.Mars*((2/r1)-(1/a1)));
gamma_1 = acos(h1/(r1*v1));
%find components of velocity 1 and express as a vector
v_r1 = v1*sin(gamma_1);   v_theta1 = v1*cos(gamma_1);
V1_roh = [v_r1;v_theta1;0];

%find semimajor axis and eccentricity of orbit 2
e2 = (ra2-rp2)/(ra2+rp2);
a2 = ra2/(1+e2);            r2 = r1;
%find specific angular momentum of orbit 2
p2 = a2*(1-e2^2);
h2 = sqrt(p2*const.mu.Mars);
%find flight path angle 2 using h, r, and v
v2 = sqrt(const.mu.Mars*((2/r2)-(1/a2)));
gamma_2 = acos(h2/(r2*v2));
%use Law of Cosines to find magnitude of delta_v
delta_v = sqrt(v1^2 + v2^2 - (2*v1*v2*cos(gamma_2 - gamma_1)));

%% Problem 3
clear; clc; const = getConst();
%givens
a_ear = const.a.Earth*const.AU;   a_sat = const.a.Saturn*const.AU;
rb = 11*const.AU;
%find a and e of Hohmann transfer orbit
a_t1 = (a_ear+a_sat)/2;
e_t1 = (a_sat-a_ear)/(a_sat+a_ear);
%find apoapsis and periapsis velocities of transfer orbit
vp_t1 = sqrt(((2*const.mu.Sun)/a_ear) - (const.mu.Sun/a_t1));
va_t1 = sqrt(((2*const.mu.Sun)/a_sat) - (const.mu.Sun/a_t1));
%find total delta_v of Hohmann transfer
delta_v_t1 = (vp_t1 - sqrt(const.mu.Sun/a_ear)) + (sqrt(const.mu.Sun/a_sat) - va_t1);
%find TOF of Hohmann transfer
TOF_t1 = pi*sqrt(a_t1^3/const.mu.Sun);
TOF_t1_years = TOF_t1/(3600*24*365);
%find initial phase angle of Hohmann transfer
halfP_sat = pi*sqrt(a_sat^3/const.mu.Sun);
phi = pi*(1-(TOF_t1/halfP_sat));

%find a and e of bielliptic transfer orbits
a_t2 = (a_ear+rb)/2;    e_t2 = (rb-a_ear)/(rb+a_ear);
a_t3 = (rb + a_sat)/2;    e_t2 = (rb-a_sat)/(rb+a_sat);
%find apoapsis and periapsis velocities of transfer orbits
vp_t2 = sqrt(((2*const.mu.Sun)/a_ear) - (const.mu.Sun/a_t2));
va_t2 = sqrt(((2*const.mu.Sun)/rb) - (const.mu.Sun/a_t2));
vp_t3 = sqrt(((2*const.mu.Sun)/a_sat) - (const.mu.Sun/a_t3));
va_t3 = sqrt(((2*const.mu.Sun)/rb) - (const.mu.Sun/a_t3));
%find total delta_v of bielliptic transfer
delta_v1 = vp_t2 - sqrt(const.mu.Sun/a_ear);    %earth circ to xfer 2 @peri
delta_v2 = va_t3 - va_t2;                       %xfer 2 @ apo to xfer 3 @apo
delta_v3 = vp_t3 - sqrt(const.mu.Sun/a_sat);    %xfer 3 @ peri to saturn circ
delta_v_t2 = delta_v1 + delta_v2 + delta_v3;
%find TOF of bielliptic transfer
TOF_t2 = pi*sqrt(a_t2^3/const.mu.Sun) + pi*sqrt(a_t3^3/const.mu.Sun);
TOF_t2_years = TOF_t2/(3600*24*365);
%find relevant ratios to inform comparison
rb_ra = rb/a_ear;
rc_ra = a_sat/a_ear;

%% Problem 4
const.mu.Mars = 42828.4;        %mu used by STK in km^3/s^2
a4 = 6600;                      %km
e4 = 0.46;      i4 = deg2rad(72.4);     RAAN4 = deg2rad(207);
omega4 = deg2rad(10);       theta_star4 = deg2rad(180);

%read in initial state values from STK
h4_i = 14928.32436827371;       %km^2/s
energy4_i = -3.244573701074832; %km^2/s^2
E4_i = deg2rad(180);
date4_i = 2457113.5;            %UTC Julian date

alt4_f = 499.9999999999993179;  %km
E4_f = -((2*pi) -deg2rad(333.1699757431793));
theta_star4_f = -((2*pi) -deg2rad(317.1716183649411));
date4_f = 2457113.58639139;     %UTC Julian date

%calculate time of flight
TOF4_julian = (date4_f - date4_i)*(3600*24);

%below 500 km
date4_f2 = 2457113.60194506;
E4_f2 = deg2rad(26.57417547688323);
%TOF below 500 km
TOF4_2_julian = (date4_f2 - date4_f)*(3600*24);
TOF4_2_eccent = sqrt(a4^3/const.mu.Mars)*((E4_f2 - e4*sin(E4_f2)) - (E4_f - e4*sin(E4_f)));
disp("TOF using STK is " + string(TOF4_2_julian) + " sec and using MATLAB is " + string(TOF4_2_eccent))

%% Define constants
function const = getConst()
    const.mu.Mars = 4.305e4;    %km^3/s^2
    const.mu.Moon = 4902.799;
    const.mu.Sun = 1.32712428e11;
    const.mu.Saturn = 3.794e7;
    const.radius.Mars = 3397.2; %km
    const.radius.Moon = 1738;
    const.a.Earth = 1.0000010178;   %AU
    const.a.Saturn = 9.554909595;
    const.AU = 149597870.7;         %km
end

