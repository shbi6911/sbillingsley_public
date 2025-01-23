clc
clear all
close all

tic
load('HW6EphemEarth')
load('HW6EphemMars')

num_rows_earth = 47; 
num_rows_mars = 41;  

c3_earth_matrix = zeros(num_rows_earth, num_rows_mars);
c3_mars_matrix = zeros(num_rows_earth, num_rows_mars);

v_inf_earth_matrix = zeros(num_rows_earth, num_rows_mars);
v_inf_mars_matrix = zeros(num_rows_earth, num_rows_mars);

%% Extract necessary columns for Earth and Mars
HW6EphemEarth_Numbers = HW6EphemEarth(:, {'VarName1', 'VarName5', 'VarName6', 'VarName7', 'VarName8', 'VarName9', 'VarName10'});
HW6EphemEarth_Matrix = table2array(HW6EphemEarth_Numbers);

HW6EphemMars_Numbers = HW6EphemMars(:, {'VarName1', 'VarName5', 'VarName6', 'VarName7', 'VarName8', 'VarName9', 'VarName10'});
HW6EphemMars_Matrix = table2array(HW6EphemMars_Numbers);

gm_sun = 1.32712428e11; 
au = 49597870.7; 

%% Nested Loops
for i = 1:num_rows_earth  
    for j = 1:num_rows_mars  

        %% Earth 
        r1_bar = [HW6EphemEarth_Matrix(i, 2), HW6EphemEarth_Matrix(i, 3), HW6EphemEarth_Matrix(i, 4)];
        v1_bar = [HW6EphemEarth_Matrix(i, 5), HW6EphemEarth_Matrix(i, 6), HW6EphemEarth_Matrix(i, 7)];
        julian_date_1 = HW6EphemEarth_Matrix(i, 1);

        %% Mars 
        r2_bar = [HW6EphemMars_Matrix(j, 2), HW6EphemMars_Matrix(j, 3), HW6EphemMars_Matrix(j, 4)];
        v2_bar = [HW6EphemMars_Matrix(j, 5), HW6EphemMars_Matrix(j, 6), HW6EphemMars_Matrix(j, 7)];
        julian_date_2 = HW6EphemMars_Matrix(j, 1);

        %% Problem Start
        angle = 0; %delta TA under 180 enter 0, 1 if over 180

        r1 = norm(r1_bar);
        v1 = norm(v1_bar);
        r2 = norm(r2_bar);
        v2 = norm(v2_bar);

        tof = (julian_date_2 - julian_date_1) * 24 * 3600;  % Time of flight in seconds
        
        if tof <= 0
            continue;  % Skip if TOF is negative or zero
        end

        %% Step 1: Change in True Anomaly
        dtheta_star = acosd(dot(r1_bar, r2_bar) / (norm(r1_bar) * norm(r2_bar)));

        if angle == 0 && dtheta_star < 180
            dtheta_star_corrected = dtheta_star;
        elseif angle == 0 && dtheta_star > 180
            dtheta_star_corrected = 360 - dtheta_star;
        elseif angle == 1 && dtheta_star < 180
            dtheta_star_corrected = 360 - dtheta_star;
        elseif angle == 1 && dtheta_star > 180
            dtheta_star_corrected = dtheta_star;
        end

        %% Step 2: Space Triangle
        c = sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cosd(dtheta_star));
        s = (r1 + r2 + c) / 2;

        %% Step 3: Elliptical or Hyperbolic Check
        tof_p = (2 / gm_sun) ^ 0.5 / 3 * (s^1.5 - (s - c)^1.5);
        
        %% Step 4
        a_min = s / 2;
        n_min = sqrt(gm_sun / a_min^3);
        alpha_min = pi;
        beta_min = 2 * asin(sqrt((s - c) / s));

        tof_min = 1 / n_min * (alpha_min - beta_min - (sin(alpha_min) - sin(beta_min)));
        
        if dtheta_star_corrected < 180
            tof_min = 1/n_min * (alpha_min - beta_min - (sin(alpha_min) - sin(beta_min)));
        elseif dtheta_star_corrected > 180
            tof_min = 1/n_min * (alpha_min - beta_min + (sin(alpha_min) - sin(beta_min)));
        end

        if tof < tof_min
            flightpath = 0; % Short flight path
        else
            flightpath = 1; % Longer flight path
        end
        

        %% Step 5: Solve Lambert's Equation
        a_input = a_min + 10000;  % Initial guess for semi-major axis
        options = optimoptions('fsolve', 'MaxIter', 1000, 'TolFun', 1e-30, 'TolX', 1e-30, 'Display', 'off');
        a_solution = fsolve(@(a) Lambers_Equation_Solver_Function(dtheta_star_corrected, flightpath, a, s, c, gm_sun, tof), a_input, options);
        a = real(a_solution);

        %% Step 6: Time Check

       
        alpha_0 = 2 * asin(sqrt(s / (2 * a)));
        beta_0 = 2 * asin(sqrt((s - c) / (2 * a)));
        n = sqrt(gm_sun / a^3);

        if dtheta_star_corrected < 180 && flightpath == 0
            alpha = alpha_0;
            beta = beta_0;
        elseif dtheta_star_corrected < 180 && flightpath == 1
            alpha = 2 * pi - alpha_0;
            beta = beta_0;
        elseif dtheta_star_corrected > 180 && flightpath == 0
            alpha = alpha_0;
            beta = -beta_0;
        elseif dtheta_star_corrected > 180 && flightpath == 1
            alpha = 2 * pi - alpha_0;
            beta = -beta_0;
        end

        tof_test = 1 / n * (alpha - beta - (sin(alpha) - sin(beta)));
        tof_check = tof_test - tof;
        
        if tof_check > 1e-3
            fprintf('Break')
            return
        end

        %% Step 7: Calculate p and e
        p = 4 * a * (s - r1) * (s - r2) / c^2 * (sin((alpha + beta) / 2))^2;
        e = sqrt(1 - p / a);

        %% Step 9: Calculate f and g funcations
        f = 1 - r2/p * (1 - cosd(dtheta_star_corrected));
        g = r2 * r1/(gm_sun * p)^.5 * sind(dtheta_star_corrected);
        
        v1_barf = (r2_bar - f * r1_bar)/g;
        
        f_dot = (gm_sun/p)^.5 * tand((dtheta_star_corrected)/2) * ((1 - cosd(dtheta_star_corrected))/p - 1/r1 - 1/r2);
        g_dot = 1 - (r1 / p) * (1 - cosd(dtheta_star_corrected));
        
        v2_bari = f_dot * r1_bar + g_dot * v1_barf;

        %% Step 10: Calculate delta Vs
        delta_v1 = norm(v1_barf - v1_bar);
        delta_v2 = norm(v2_bar - v2_bari);

        v_inf_earth = delta_v1;
        v_inf_mars = delta_v2;

        c3_earth = delta_v1^2;
        c3_mars = delta_v2^2;

        c3_earth_matrix(i, j) = c3_earth;
        c3_mars_matrix(i, j) = c3_mars;
        
%         if v_inf_earth >= 10;
%             v_inf_earth = NaN;
%         end
%         if v_inf_mars >= 7;
%             v_inf_mars = NaN;
%         end
        
        v_inf_earth_matrix(i, j) = v_inf_earth;
        v_inf_mars_matrix(i, j) = v_inf_mars;

    end
end

%% Plot Pork Chop Plot Stuff

% Julian days to dates
earth_departure_dates = datetime(HW6EphemEarth_Matrix(1:num_rows_earth, 1), 'convertfrom', 'juliandate');
mars_arrival_dates = datetime(HW6EphemMars_Matrix(1:num_rows_mars, 1), 'convertfrom', 'juliandate');

julian_earth_departure = HW6EphemEarth_Matrix(1:num_rows_earth, 1); 
julian_mars_arrival = HW6EphemMars_Matrix(1:num_rows_mars, 1);      

[X, Y] = meshgrid(julian_earth_departure, julian_mars_arrival);

%% Figures

% Choose every # of ticks
xticks = 1:10:num_rows_earth;  
yticks = 1:5:num_rows_mars; 

figure(1)
[C, h] = contour(X, Y, v_inf_earth_matrix', [0:.25:6]);
title('V_{\infty} (km/s) Earth Pork Chop Plot');
xlabel('Earth Departure Date');
ylabel('Mars Arrival Date');
colorbar
 

set(gca, 'XTick', julian_earth_departure(xticks));  
set(gca, 'YTick', julian_mars_arrival(yticks));     

xticklabels(datestr(earth_departure_dates(xticks), 'mmm dd, yyyy'));
yticklabels(datestr(mars_arrival_dates(yticks), 'mmm dd, yyyy'));
% [C, h] = contour(X, Y, v_inf_earth_matrix', [0:.25:6]);
clabel(C, h, 'FontSize', 8, 'Color', 'k', 'LabelSpacing',1000); 
colormap('winter')
colorbar
grid on;

figure(2)
[g, j] = contour(X, Y, v_inf_mars_matrix', [0:.25:5]);
title('V_{\infty} (km/s) Mars Pork Chop Plot');
xlabel('Earth Departure Date');
ylabel('Mars Arrival Date');


set(gca, 'XTick', julian_earth_departure(xticks));  
set(gca, 'YTick', julian_mars_arrival(yticks));     

xticklabels(datestr(earth_departure_dates(xticks), 'mmm dd, yyyy'));
yticklabels(datestr(mars_arrival_dates(yticks), 'mmm dd, yyyy'));

% [g, j] = contour(X, Y, v_inf_mars_matrix', [0:.25:5]);
clabel(g, j, 'FontSize', 8, 'Color', 'k', 'LabelSpacing',5000);
colormap('winter')
colorbar
grid on;

toc
