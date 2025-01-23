%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc; const = getConst();

%characterize initial orbit
r_p1 = 1850;    r_a1 = 5400;    r1 = 3500;      %km
e1 = (r_a1 - r_p1)/(r_a1 + r_p1);   %eccentricity
a1 = (r_a1+r_p1)/2;                 %semimajor axis
p1 = a1*(1 - e1^2);                 %semilatus rectum
h1 = sqrt(const.mu.Moon*p1);        %specific angular momentum
v1 = sqrt(((2*const.mu.Moon)/r1) - (const.mu.Moon/a1)); %velocity
gamma1 = -acos(h1/(r1*v1));         %flight path angle
V1 = [v1*sin(gamma1);v1*cos(gamma1);0];     %velocity vector
theta_star1 = -acos((p1 - r1)/(r1*e1));  %true anomaly

%characterize final orbit
r_p2 = 3000;    r_a2 = 5400;    r2 = 3500;      %km
e2 = (r_a2 - r_p2)/(r_a2 + r_p2);   %eccentricity
a2 = (r_a2+r_p2)/2;                 %semimajor axis
p2 = a2*(1 - e2^2);                 %semilatus rectum
h2 = sqrt(const.mu.Moon*p2);        %specific angular momentum
v2 = sqrt(((2*const.mu.Moon)/r2) - (const.mu.Moon/a2)); %velocity
gamma2 = -acos(h2/(r2*v2));         %flight path angle
V2 = [v2*sin(gamma2);v2*cos(gamma2);0];     %velocity vector
theta_star2 = -acos((p2 - r2)/(r2*e2));  %true anomaly

%find delta_v
delta_v = sqrt(v1^2 + v2^2 - (2*v1*v2*cos(gamma2 - gamma1)));
delta_v_check = norm(V2 - V1);
%find change in argument of periapsis
delta_w = theta_star1 - theta_star2;


%% Problem 2
clear; clc; const = getConst();

%characterize initial orbit
Vvec_in = [2.2789;5.8841;0];    %km/s
a1 = 1.8e6;                     %semimajor axis of initial orbit
gamma = atan2(Vvec_in(1),Vvec_in(2));   %flight path angle
denom = dot(Vvec_in,Vvec_in)/2 + const.mu.Jupiter/(2*a1);
r1 = const.mu.Jupiter/denom;        %orbital radius from vis-viva
h1 = norm(cross([r1;0;0],Vvec_in)); %specific angular momentum
e1 = sqrt(1 - (h1^2/(const.mu.Jupiter*a1)));    %eccentricity   
p1 = a1*(1 - e1^2);                 %semilatus rectum
theta_star1 = acos((p1 - r1)/(r1*e1));  %true anomaly

v_moon = sqrt(const.mu.Jupiter/r1);     %orbital speed of moon X
V_infty_in = Vvec_in - [0;v_moon;0];       %velocity relative to moon

%characterize new orbit
v_out = 7.6759;     gamma_out = deg2rad(20.91);
Vvec_out = [v_out*sin(gamma_out);v_out*cos(gamma_out);0];
h2 = norm(cross([r1;0;0],Vvec_out)); %specific angular momentum
denom2 = dot(Vvec_out,Vvec_out) - ((2*const.mu.Jupiter)/r1);
a2 = -const.mu.Jupiter/denom2;      %semimajor axis
e2 = sqrt(1 - (h2^2/(const.mu.Jupiter*a2)));    %eccentricity   
p2 = a2*(1 - e2^2);                 %semilatus rectum
theta_star2 = acos((p2 - r1)/(r1*e2));  %true anomaly

V_infty_out = Vvec_out - [0;v_moon;0]; %velocity relative to moon after flyby

%characterize hyperbola
mu_moon = 4000;
a_h = -mu_moon/norm(V_infty_in)^2;  %semimajor axis
%turning angle as angle between incoming and outgoing vectors
delta = acos(dot(V_infty_in,V_infty_out)/(norm(V_infty_in)*norm(V_infty_out)));
e_h = 1/(sin(delta/2));             %eccentricity
r_p_h = a_h*(1-e_h);                %periapsis radius

%characterize impulsive maneuver
delta_Vvec = Vvec_out - Vvec_in;
delta_v = norm(delta_Vvec);
m_i = 2500; %initial mass
Isp = 300;  %specific impulse
g0 = 9.81;  %acceleration of gravity on Earth
m_f = m_i*exp(-delta_v/((Isp*g0)/1000));

%% Problem 3
clear; clc; const = getConst();

Vvec_in = [3.2476;13.1175;0];
Vvec_out = [3.3426;13.7356;0];
V_Jupiter = [0;sqrt(const.mu.Sun/(const.a.Jupiter*const.AU));0];
V_infty_in = V_Jupiter - Vvec_in;
V_infty_out = V_Jupiter - Vvec_out;

%characterize hyperbola
a_h = -const.mu.Jupiter/norm(V_infty_in)^2;  %semimajor axis
%turning angle as angle between incoming and outgoing vectors
delta = acos(dot(V_infty_in,V_infty_out)/(norm(V_infty_in)*norm(V_infty_out)));
e_h = 1/(sin(delta/2));             %eccentricity
r_p_h = a_h*(1-e_h);                %periapsis radius