clc
clear all
close all

tic

gm_sun = 1.32712428e11;
au = 49597870.7;

% Desired delta theta
angle = 0; %delta TA  under 180 enter 0, 1 if over 180

julian_date_1 = 2457754.5;
julian_date_2 = 2457871.5;

r1_bar = [-2.686982e7 1.326980e8 5.752566e7];
v1_bar = [-29.781722 -5.101774 -2.210394];
r1 = norm(r1_bar);
v1 = norm(v1_bar);

r2_bar = [-5.642901e7 -8.571048e7 -3.499466e7];
v2_bar = [29.658341 -16.091100 -9.116674];
r2 = norm(r2_bar);
v2 = norm(v2_bar);

tof = (julian_date_2 - julian_date_1) * 24 * 3600
%% step 1 Change in True Anamoly
dtheta_star = acosd(dot(r1_bar, r2_bar)/(norm(r1_bar) * norm(r2_bar)))

if angle == 0 && dtheta_star < 180
    dtheta_star_corrected = dtheta_star
elseif angle == 0 && dtheta_star > 180
    dtheta_star_corrected = 360 - dtheta_star
elseif angle == 1 && dtheta_star < 180
    dtheta_star_corrected = 360 - dtheta_star
elseif angle == 1 && dtheta_star > 180
    dtheta_star_corrected = dtheta_star
end


%% Step 2 Calculating Geometric Quantities
c = (r1^2 + r2^2 - 2 * r1 * r2 * cosd(dtheta_star))^(1/2)
s = (r1 + r2 + c)/2
%% step 3 Elliptical or Hyperbolic
tof_p = 1/3 * (2/gm_sun)^.5 * (s^(3/2) - (s - c)^(3/2))
tof_check = tof - tof_p  

%% Step 4  Shorter or Longer TOF

a_min = s/2;
n_min = (gm_sun/a_min^3)^(1/2);
alpha_min = pi;
beta_min = 2 * asin (((s - c)/s)^.5)

if dtheta_star_corrected < 180
    tof_min = 1/n_min * (alpha_min - beta_min - (sin(alpha_min) - sin(beta_min)))
elseif dtheta_star_corrected > 180
    tof_min = 1/n_min * (alpha_min - beta_min + (sin(alpha_min) - sin(beta_min)))
end


if tof < tof_min
    flightpath = 0%short
else
    flightpath = 1 %longer
end

%% Step 5
a_input = a_min + 10000

mu = gm_sun

%Solution to Lambert's Equation
options = optimoptions('fsolve', 'MaxIter', 1000, 'TolFun', 1e-8, 'StepTolerance', 1e-8, 'Display', 'iter')
a_solution = fsolve(@(a) Lambers_Equation_Solver_Function(dtheta_star_corrected, flightpath, a, s, c, mu, tof), a_input, options);
a = real(a_solution)

%% Step 6
alpha_0 = 2 * asin(sqrt(s / (2 * a)));
beta_0 = 2 * asin(sqrt((s - c) / (2 * a)));
n = sqrt(mu / a^3);

if dtheta_star_corrected < 180 && flightpath == 0
    alpha = alpha_0
    beta = beta_0
elseif dtheta_star_corrected < 180 && flightpath == 1
    alpha = 2*pi - alpha_0
    beta = beta_0
elseif dtheta_star_corrected > 180 && flightpath == 0
    alpha = alpha_0
    beta = -beta_0
elseif dtheta_star_corrected > 180 && flightpath == 1
    alpha = 2*pi - alpha_0
    beta = -beta_0
end

tof_test = 1/n*(alpha - beta - (sin(alpha) - sin(beta)))
tof_check = tof_test - tof

if tof_check > 1e-3
    fprintf('Break')
    return
end
%% Step 7

p = 4 * a * (s - r1) * (s - r2)/c^2 * (sin((alpha + beta)/2))^2;
e = (1 - p/a)^.5

%% Step 8

theta_star_1_i = acosd(1/e * (p/r1 - 1))
theta_star_2_i = acosd(1/e * (p/r2 - 1))

theta_check = theta_star_2_i - theta_star_1_i

tolerance = 1e-8;
found = false;
dtheta_star_loop = 0;
for sign1 = [-1, 1]
    for sign2 = [-1, 1]
        if dtheta_star_corrected > 180
            dtheta_star_loop = dtheta_star_corrected - 360;
        else
            dtheta_star_loop = dtheta_star_corrected
        end
        theta1 = sign1 * theta_star_1_i;
        theta2 = sign2 * theta_star_2_i;
        
        fprintf('theta1 = %.4f, theta2 = %.4f\n', theta1, theta2);
        if abs(theta2 - theta1 - dtheta_star_loop) < tolerance
            found = true;
            break; 
        end
    end
    if found
        break;
    end
end

theta_star_1 = theta1
theta_star_2 = theta2

%% Step 9

f = 1 - r2/p * (1 - cosd(dtheta_star_corrected))
g = r2 * r1/(mu * p)^.5 * sind(dtheta_star_corrected)

v1_barf = (r2_bar - f * r1_bar)/g

f_dot = (mu/p)^.5 * tand((dtheta_star_corrected)/2) * ((1 - cosd(dtheta_star_corrected))/p - 1/r1 - 1/r2)
g_dot = 1 - (r1 / p) * (1 - cosd(dtheta_star_corrected))

v2_bari = f_dot * r1_bar + g_dot * v1_barf

%% Step 10 

delta_v1 = (v1_barf - v1_bar)
delta_v2 = (v2_bar - v2_bari)

v_inf_earth = norm(delta_v1)
v_inf_mars = norm(delta_v2)

c3_earth = v_inf_earth^2;
c3_mars = v_inf_mars^2;





toc




