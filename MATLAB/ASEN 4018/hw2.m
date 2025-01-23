%By:    Adrian Bryant
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

clear
clc

Re = 6387.1;
altitude = 0:999;

for x = 1:1000
    area(x) = (4)*pi*(Re+x)^2;
    n = MSISatmosphere1000(altitude(x));
    mass(x) = n.mass * 1e12;
end

mass = trapz(altitude,(mass.*area));



Re = 6387.1;
altitude = 0:999;

for x = 1:1000
    area(x) = (4)*pi*(Re+x)^2 * 1e6;
    b = MSISatmosphere1000(altitude(x));
    n_tot(x) = b.total * 1e6;
    
end

T = 300;
R = 8.3145 / 6.022e23;
V = 10*10*100*1e9;

n_for_all = trapz(altitude*1e3,area .* n_tot);

P = (n_for_all * R * T) / V;

%----
avogadro = 6.022e23;
k_b = 1.38e-23; % J/k
T = 300;
g = 9.81;
scale_h = @(m,g,T)  k_b*T/(m*g);
% 0.215 O2 and 0.785 N2
m_air = 0.215*(32/avogadro*1e-3) + 0.785*(28/avogadro*1e-3);
m_o = (16 / avogadro*1e-3);
m_he = (4 / avogadro*1e-3);
m_h = (1 / avogadro*1e-3);

fprintf("Scale Heights at 300 K: \n-------- \n")
fprintf("Scale Height of N2/O2 Mix: %d m \n",scale_h(m_air,g,T))
fprintf("Scale Height of O : %d m \n",scale_h(m_o,g,T))
fprintf("Scale Height of He : %d m \n",scale_h(m_he,g,T))
fprintf("Scale Height of H : %d m \n",scale_h(m_h,g,T))

%---
% Scale Heights Quantitatively
scale_h_est = @(z0,z1,n0,n1) (z0 - z1) / (log(n1) - log(n0)) ;

for x=1:999
    H_rel(x) = scale_h_est(x,x+1,n_tot(x),n_tot(x+1));
end

figure(3)
plot(log(n_tot),1:1000)
xlabel("Ln(Total Number Density)")
ylabel("Height [km]")

figure(4)
plot(H_rel,1:999)
xlabel("Scale Height Estimation for 1 km intervals [km]")
ylabel("Height [km]")

%---
P=0;

year = 60 * 60 * 24 * 365; % minutes
P_400k = 92.42 * 60; % seconds
P_500k = 94.47 * 60; % seconds
area = [30*30; 30*30; 10*10] * 0.0001; % m
p_opt = [P_400k, P_500k, P_400k];
tspan = [0, 5*year];


options = odeset('RelTol',0.4,'Events',@done);
[t,p] = ode45(@(t,X) dp_dt(t,X,area(1)), tspan, p_opt(1), options);
t_day_1 = t./(60*60*24);
p_min_1 = p./(60);
alt_1=( ((p.^2) .* 5.97e24 .* 6.674e-11) ./ (4 * pi^2) ).^(1/3)./1000 - Re;

[t,p] = ode45(@(t,X) dp_dt(t,X,area(2)), tspan, p_opt(2), options);
t_day_2 = t./(60*60*24);
p_min_2 = p./(60);

p_min_2 = p_min_2((p_min_2 > 87) & (p_min_2 < 95));
t_day_2 = t_day_2((p_min_2 > 87) & (p_min_2 < 95));
alt_2=( ((p.^2) .* 5.97e24 .* 6.674e-11) ./ (4 * pi^2) ).^(1/3) ./1000 - Re;
alt_2 = alt_2((p_min_2 > 87) & (p_min_2 < 95));

[t,p] = ode45(@(t,X) dp_dt(t,X,area(3)), tspan, p_opt(3), options);
t_day_3 = t./(60*60*24);
p_min_3 = p./(60);

p_min_3 = p_min_3( ~any( isnan( p_min_3 ) | isinf( p_min_3 ), 2 ),: );
t_day_3 = t_day_3( ~any( isnan( p_min_3 ) | isinf( p_min_3 ), 2 ),: );
alt_3=( ((p.^2) .* 5.97e24 .* 6.674e-11) ./ (4 * pi^2) ).^(1/3) ./1000 - Re;

p_min_3 = p_min_3((p_min_3 > 87) & (p_min_3 < 95));
alt_3 = alt_3((p_min_3 > 87) & (p_min_3 < 95));
t_day_3 = t_day_3((p_min_3 > 87) & (p_min_3 < 95));
figure(1)
hold on
grid on
plot(t_day_1, p_min_1)
plot(t_day_2, p_min_2)
plot(t_day_3, p_min_3)
title("Orbital Period Decay Due to Drag")
xlabel("Time [Days]")
ylabel("Orbital Period [Minutes]")
ylim([87.5 95])
legend(["400km Starting Orbit, +X Orientation", ...
    "500km Starting Orbit, +X Orientation", "400km Starting Orbit," + ...
    " +Z Orientation"])



figure(5)
hold on
grid on
plot(t_day_1, alt_1)
plot(t_day_2, alt_2)
plot(t_day_3, alt_3)
title("Orbital Altitude Decay Due to Drag")
xlabel("Time [Days]")
ylabel("Altitude [km]")
ylim([140 500])
legend(["400km Starting Orbit, +X Orientation", ...
    "500km Starting Orbit, +X Orientation", "400km Starting Orbit," + ...
    " +Z Orientation"])

fprintf("400km Starting Orbit, +X Orientation Deorbits in %i days\n", ceil(min(t_day_1(alt_1 < 140))))
fprintf("500km Starting Orbit, +X Orientation Deorbits in %i days\n", ceil(min(t_day_2(alt_2 < 140))))
fprintf("400km Starting Orbit, +Z Orientation Deorbits in %i days\n", ceil(min(t_day_3(alt_3 < 140))))

%---
options = odeset('RelTol',0.5,'Events',@done);
tspannew = (1:1:200) * 24 * 60 * 60;
err = -0.1;
[t,p] = ode45(@(t,X) dp_dt_withERROR(t,X,area(1),err), tspannew, p_opt(1), options);
t_err_day_1 = t./(60*60*24);
p_err_min_1 = p./(60);
alt_err_1=( ((p.^2) .* 5.97e24 .* 6.674e-11) ./ (4 * pi^2) ).^(1/3) ./1000 - Re;

options = odeset('RelTol',0.8,'Events',@done);

err = 0.1;
[t,p] = ode45(@(t,X) dp_dt_withERROR(t,X,area(1),err), tspan, p_opt(1), options);
t_err_day_2 = t./(60*60*24);
p_err_min_2 = p./(60);
alt_err_2=( ((p.^2) .* 5.97e24 .* 6.674e-11) ./ (4 * pi^2) ).^(1/3) ./1000 - Re;

figure(2)
hold on
grid on
plot(t_day_1, p_min_1)
plot(t_err_day_1, p_err_min_1)
plot(t_err_day_2,p_err_min_2)

title("\pm 10% Error within Density, Orbital Decay due to Drag")
xlabel("Time (days)")
ylabel("Period (Minutes)")
legend("No Error","-10% Error in Density","+10% Error in Density")

figure(6)
hold on
grid on
plot(t_day_1, alt_1)
plot(t_err_day_1, alt_err_1)
plot(t_err_day_2,alt_err_2)

title("\pm 10% Error within Density, Orbital Decay due to Drag")
xlabel("Time (days)")
ylabel("Altitude (Km)")
legend("No Error","-10% Error in Density","+10% Error in Density")



function [dP] = dp_dt_withERROR(t,P, Area, err)
    C_d = 2.2; 
    m = 4; % kg -- kinda arbitrary
    a = ( (P^2 * 5.97e24 * 6.674e-11) / (4 * pi^2) )^(1/3);
    if a < 100
        dP = 0
        return
    end
    params = MSISatmosphere1000((a/1000) - 6371);
    density = (params.mass + params.mass*err)  * 1000;
    
    dP = -3 * pi * a * density * C_d * Area / m;
end

function [value, isterminal, direction] = done(T, P)
    value      = (P < 85 || P > 95);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end

function [dP,a] = dp_dt(t,P, Area)
    C_d = 2.2; 
    m = 4; % kg -- kinda arbitrary
    a = ( (P^2 * 5.97e24 * 6.674e-11) / (4 * pi^2) )^(1/3);
    if a < 100
        return
    end
    params = MSISatmosphere1000((a/1000) - 6371);
    density = params.mass * 1000;
    
    dP = -3 * pi * a * density * C_d * Area / m;
end