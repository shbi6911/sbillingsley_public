%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Test Functionalized Code
clear; clc; const = getConst();
state1 = [-2.686982e7;1.326980e8;5.752566e7;-29.781722;-5.101774;-2.210394];   %km (Sun-centered)
state2 = [-5.642901e7;-8.571048e7;-3.499466e7;29.658341;-16.091100;-9.116674]; %km, km/s
TOF_1 = (2457871.5 - 2457754.5)*(3600*24);
%angleflag is 1 for delta_theta_star < 180, 2 for delta_theta_star > 180
angleflag = 1;

[state1_t, state2_t] = transfer(state1,state2,TOF_1,angleflag);

delta_V_1 = state1_t(4:6) - state1(4:6);
delta_V_2 = state2_t(4:6) - state2(4:6);
delta_v_tot = norm(delta_V_1) + norm(delta_V_2);

%% Run Through Transfer Parameter Space
clear; clc; const = getConst();

Earth_data = readmatrix('HW6_Ephem_Earth.txt','NumHeaderLines',3,'Delimiter',',');
Earth_data(:,2) = [];
Mars_data = readmatrix('HW6_Ephem_Mars.txt','NumHeaderLines',3,'Delimiter',',');
Mars_data(:,2) = [];
angleflag = 1;
%preallocate results arrays
departure_data = zeros(length(Mars_data),length(Earth_data));
arrival_data = zeros(length(Mars_data),length(Earth_data));
for i = 1:length(Earth_data)
    for j = 1:length(Mars_data)
        TOF = (Mars_data(j,1) - Earth_data(i,1))*3600*24;   %find time of flight
        [dep_state,arr_state] = transfer(Earth_data(i,2:7)', Mars_data(j,2:7)',TOF,angleflag);
        departure_data(j,i) = norm(dep_state(4:6) - Earth_data(i,5:7)');
        arrival_data(j,i) = norm(arr_state(4:6) - Mars_data(j,5:7)');
    end
end

figure();
[X,Y] = meshgrid(Earth_data(:,1) - min(Earth_data(:,1)),...
    Mars_data(:,1) - min(Earth_data(:,1)));
[C1,h1] = contour(X,Y,departure_data,[0:0.25:5,5.5:0.5:10]);
title("Porkchop Plot for Earth Departure V_{\infty}");
xlabel("Earth departure Epoch (Days past JD 2453541.5)");
ylabel("Mars Arrival Epoch (Days past JD 2453541.5)");
clabel(C1, h1, 'FontSize', 8, 'Color', 'k', 'LabelSpacing',1000); 
colormap('turbo')
colorbar

figure();
[C2,h2] = contour(X,Y,arrival_data,[0:0.25:5,5.5:0.5:10]);
title("Porkchop Plot for Mars Arrival V_{\infty}");
xlabel("Earth departure Epoch (Days past JD 2453541.5)");
ylabel("Mars Arrival Epoch (Days past JD 2453541.5)");
clabel(C2, h2, 'FontSize', 8, 'Color', 'k', 'LabelSpacing',1000); 
colormap('turbo')
colorbar

function [state1_t,state2_t] = transfer(state1,state2,TOF,angleflag)
%TRANSFER calculates a transfer orbit between two points given orbital
%state vectors at those points, a time of flight, and an angle flag
%determining whether a transfer angle greater than or less than 180 degrees
%is desired.  Currently hardcoded to consider heliocentric elliptical transfers
%
%INPUTS     state1,2        orbital state vectors [Rx;Ry;Rz;Vx;Vy;Vz] at
%                           positions 1 and 2
%           TOF             time of flight in seconds
%           angleflag       1 = < 180 degrees, 2 = > 180 degrees
%
%OUTPUTS    state1,2_t      orbital state vectors as above, but at
%                           positions 1 and 2 on the TRANSFER orbit

%set constants
const = getConst();
%break out input variables
Rvec_1 = state1(1:3);   R1 = norm(Rvec_1);
Rvec_2 = state2(1:3);   R2 = norm(Rvec_2);
%calculate transfer angle
if angleflag == 1
    delta_theta_star = acos(dot(Rvec_1,Rvec_2)/(R1*R2));
elseif angleflag == 2
    delta_theta_star = (2*pi) - acos(dot(Rvec_1,Rvec_2)/(R1*R2));
else
    disp("Angleflag not set correctly");
    return
end

%calculate geometric quantities
c = norm(Rvec_1 - Rvec_2);
s = (R1 + R2 + c)/2;
a_min = s/2;
% find parabolic TOF
TOF_p = (1/3)*sqrt(2/const.mu.Sun)*(s^(3/2) - (s-c)^(3/2));

%find TOF of minimum energy transfer
n_min = sqrt(const.mu.Sun/a_min^3);
alpha_min = pi;
if angleflag == 1
    beta_min = 2*asin(sqrt((s-c)/s));
elseif angleflag == 2
    beta_min = -2*asin(sqrt((s-c)/s));
else
    disp("Angleflag not set correctly");
    return
end
TOF_min = (1/n_min)*((alpha_min - beta_min) - (sin(alpha_min) - sin(beta_min)));

%set initial guess and solve
guess = a_min + (0.5*a_min);
fun = @(a)lambert(a,TOF,TOF_min,const.mu.Sun,s,c,angleflag);
opts = optimoptions(@fsolve,'FunctionTolerance',1e-10,'Display','off');
a_soln = fsolve(fun,guess,opts);
a = real(a_soln);

%set alpha_0 and beta_0 for solved a
alpha_0 = 2*asin(sqrt(s/(2*a)));
beta_0 = 2*asin(sqrt((s-c)/(2*a)));
   
%determine appropriate alpha and beta values
if angleflag == 1 && TOF < TOF_min
    alpha = alpha_0;            beta = beta_0;
elseif angleflag == 2 && TOF < TOF_min
    alpha = alpha_0;            beta = -beta_0;
elseif angleflag == 1 && TOF > TOF_min
    alpha = (2*pi) - alpha_0;   beta = beta_0;
elseif angleflag == 2 && TOF > TOF_min
    alpha = (2*pi) - alpha_0;   beta = -beta_0;
else
    disp("Error (final a)");
end
%find p and e
p = ((4*a*(s - R1)*(s - R2))/c^2)*(sin((alpha+beta)/2))^2;
e = sqrt(1-(p/a));

%angle1 = acos((1/e)*((p/R1)-1));
%angle2 = acos((1/e)*((p/R2)-1));
%[theta_star_1,theta_star_2] = findangle(angle1,angle2,delta_theta_star,10^-2);

f = 1-(R2./p).*(1-cos(delta_theta_star));
g = (R2.*R1.*sin(delta_theta_star))./sqrt(const.mu.Sun*p);
coeff = ((1-cos(delta_theta_star))./p) - (1./R2) - (1./R1); 
f_dot = sqrt(const.mu.Sun/p).*tan(delta_theta_star./2).*coeff;
g_dot = 1 - (R1./p).*(1-cos(delta_theta_star));

Vvec1_t = (Rvec_2 - f.*Rvec_1)./g;
Vvec2_t = f_dot.*Rvec_1 + g_dot.*Vvec1_t;

% delta_V_1 = Vvec1_t - Vvec_1;
% delta_V_2 = Vvec2_t - Vvec_2;
% delta_v_tot = norm(delta_V_1) + norm(delta_V_2);

%set output states
state1_t = [Rvec_1;Vvec1_t];   state2_t = [Rvec_2;Vvec2_t];

end

function val = lambert(a,TOF,TOF_min,mu,s,c,angleflag)
%LAMBERT computes Lambert's equation to find the semi-major axis of an 
% elliptical transfer orbit between two points, defined by an input time of
% flight with given geometric properties.  The function is intended to be
% used with fsolve or similar rootfinding algorithm.  As such, it outputs a
% value which should converge to zero.
%
%INPUTS     a           guess for semi-major axis for this iteration
%           TOF         fixed time of flight for the transfer
%           TOF_min     fixed time of flight for minimum energy transfer
%           mu          fixed gravitational parameter for central body of
%                       the transfer orbit (km^3/s^2)
%           s           semiperimeter of the space triangle
%           c           chord of the space triangle
%           angleflag   flag to determine whether a transfer angle greater
%                       than (2) or less than (1) 180 degrees is desired
%OUTPUTS    val         a value of Lambert's equation with time of flight
%                       subtracted, such that at an appropriate value of a, 
%                       val=0

    %set alpha_0 and beta_0
    alpha_0 = 2*asin(sqrt(s/(2*a)));
    beta_0 = 2*asin(sqrt((s-c)/(2*a)));
   
    %determine appropriate alpha and beta values
    if angleflag == 1 && TOF < TOF_min
        alpha = alpha_0;            beta = beta_0;
    elseif angleflag == 2 && TOF < TOF_min
        alpha = alpha_0;            beta = -beta_0;
    elseif angleflag == 1 && TOF > TOF_min
        alpha = (2*pi) - alpha_0;   beta = beta_0;
    elseif angleflag == 2 && TOF > TOF_min
        alpha = (2*pi) - alpha_0;   beta = -beta_0;
    else
        disp("Error (lambert)");
        return;
    end
    %calculate transfer mean motion
     n = sqrt(mu/a^3);
    %find value of Lambert's equation
    val = (1/n)*((alpha - beta) - (sin(alpha) - sin(beta))) - TOF;
end

function [theta_star_1, theta_star_2] = findangle(angle1,angle2,delta_theta_star,tol)
%FINDANGLE determines the combination of two true anomalies fed to it which
%create an input transfer angle, within an input tolerance
%INPUTS     angle1      first true anomaly in radians
%           angle2      second true anomaly in radians
%                       these angles should be ordered (2=final 1=initial)
%           delta_theta_star    target transfer angle between the true
%                               anomalies
%           tol         tolerance defining how close the match should be
%OUTPUTS
%           theta_star_1,2  the combination of true anomalies for which 
%                           angle2 - angle1 produces the given transfer angle

    angle_2 = wrapTo2Pi([angle2, -angle2, angle2, -angle2]);%convert to 360
    angle_1 = wrapTo2Pi([angle1, -angle1, -angle1, angle1]);
    deltas = angle_2 - angle_1;                             %find differences
    index = find(abs(deltas - delta_theta_star) < tol); %detect match
    theta_star_1 = angle_1(index);      %set output angles
    theta_star_2 = angle_2(index);
end