%By:        Shane Billingsley
%Class:     ASEN 3801 Aerospace Sciences Lab II
%Date:      Fall 2023

%this lab was divided and coded separately, thus this script was written
%entirely by Shane Billingsley

%this is a script to solve for the coefficients of Prandtl's Lifting Line
%theory and perform various calculations modeling a given wing geometry
%under this theory

%% initialize constants
clear; tic;
const = getConst();

%% compare wing efficiency with varying aspect ratio and taper ratio

%initialize vectors for aspect ratio comparison (see Fig 5.20 Anderson)
taper = linspace(0,1,100); %values of taper ratio (c_t/c_r)
ratio = [4 6 8 10];         %values of aspect ratio AR
delta = zeros(length(ratio),length(taper)); %matrix to store delta values

%copy const struct to make changes
const_ARtaper = getConst();
%set generic values for comparison of aspect and taper ratios
const_ARtaper.a0_r = 2*pi;  
const_ARtaper.a0_t = 2*pi;
const_ARtaper.geo_r = 3;
const_ARtaper.geo_t = 3;
const_ARtaper.aero_r = 0;
const_ARtaper.aero_t = 0;
const_ARtaper.c_r = 1;

for k = 1:length(ratio) %loop through values of aspect ratio
    for i = 1:length(taper) %loop through values of taper ratio
        %calculate tip chord using taper ratio
        const_ARtaper.c_t = const_ARtaper.c_r*taper(i); 
        %calculate span for this taper ratio
        const_ARtaper.b = (ratio(k)/2)*(const_ARtaper.c_r+const_ARtaper.c_t);
        %calculate planform area
        const_ARtaper.S = (const_ARtaper.b/2)*(const_ARtaper.c_r+const_ARtaper.c_t);
        [e,~,~] = PLLT(const_ARtaper,50); %find span efficiency factor
        delta(k,i) = (1/e)-1;   %back out delta value
    end
end

%plot comparison results (replicating Fig 5.20 Anderson)
figure (8);
color = ['r','b','k','m'];
hold on;
for j = 1:length(ratio)
    plot(taper,delta(j,:),'Color',color(j));    %plot each aspect ratio
end
xlabel ("Taper Ratio $\frac{c_{t}}{c_{r}}$",'Interpreter','latex');
ylabel ("\delta",'Rotation',0);
title("Wing Efficiency vs. Taper Ratio with Increasing Aspect Ratio");
legend("AR=4","AR=6","AR=8","AR=10");
hold off;

%% find N-values for convergence

%display convergence values for reference
disp("The span efficiency of this wing is "+string(const.e_abs));
disp("The lift coefficient is "+string(const.c_l_abs));
disp("The lift is "+string(const.lift_abs)+" lbs.");
disp("The coefficient of induced drag is "+string(const.c_di_abs));
disp("The drag is "+string(const.drag_i_abs)+" lbs.");

%use custom errorfinder function to determine convergence
lift_err_10 = errorfinder(const,[1 1000],10,2);    %converge to 10%
disp("Lift converges to 10% error at N= " + string(lift_err_10));
drag_err_10 = errorfinder(const,[1 1000],10,3);
disp("Drag converges to 10% error at N= " + string(drag_err_10));

lift_err_1 = errorfinder(const,[1 1000],1,2);      %converge to 1%
disp("Lift converges to 1% error at N= " + string(lift_err_1));
drag_err_1 = errorfinder(const,[1 1000],1,3);
disp("Drag converges to 1% error at N= " + string(drag_err_1));

lift_err_01 = errorfinder(const,[1 1000],0.1,2);   %converge to 0.1%
disp("Lift converges to 0.1% error at N= " + string(lift_err_01));
drag_err_01 = errorfinder(const,[1 1000],0.1,3);
disp("Drag converges to 0.1% error at N= " + string(drag_err_01));

%% plot convergence according to previously found values

%find size of plot (1.5 times the largest N-value needed for convergence)
size = round((max([lift_err_01 drag_err_01]))*1.5);    
%initialize vectors for PLLT values
e = zeros(1,size);
c_l = zeros(1,size);
c_di = zeros(1,size);

for i = 1:size      %run PLLT for each N-value
    [e(i),c_l(i),c_di(i)] = PLLT(const,i);
end

%turn coefficients into actual lift and drag values
lift_conv = c_l.*const.q.*const.S;  
drag_conv = c_di.*const.q.*const.S;

figure(9);  %plot convergence for lift
hold on;
%plot(lift,'o');
plot(lift_conv);
yline(const.lift_abs);
title("Lift Convergence");
ylabel("Lift (lb)");
xlabel("N-value used in Fourier series approximation");
legend("PLLT Lift Output","Convergence Value (N=1000)");
hold off;

figure(10);  %plot convergence for drag
hold on;
%plot(drag,'o');
plot(drag_conv);
yline(const.drag_i_abs);
title("Drag Convergence");
ylabel("Drag (lb)");
xlabel("N-value used in Fourier series approximation");
legend("PLLT Induced Drag Output","Convergence Value (N=1000)");
hold off;

%% compare aerodynamic efficiency over varying angle of attack

%copy const struct to make changes
const_LD = getConst();

%preallocate various vectors
%create vector of alpha (Angle of Attack) values
alpha = linspace(-5,12,50);
%create vector of lift coefficients for NACA0012
coeff_lift1 = linspace(-1.2,1.2,50);
%create vector of lift coefficients for NACA2412
coeff_lift2 = linspace(-0.6,1.5,50);
%create vector of induced drag coefficients at AoA values
wing_cdi = zeros(1,length(alpha));
%create vector of lift coefficients at AoA values
wing_cl = zeros(1,length(alpha));

load("NACAdata.mat");   %load experimental data

%generate polyfit line for NACA0012 lift vs AoA
NACA0012_a0 = polyfit(NACA0012_cl_alpha(:,1),NACA0012_cl_alpha(:,2),5);
figure (11); %plot polyfit line against data to check compliance
hold on;
scatter(NACA0012_cl_alpha(:,1),NACA0012_cl_alpha(:,2),'o','b');
plot(alpha,polyval(NACA0012_a0,alpha),'r');
title("Coefficient of Lift vs. Angle of Attack for a NACA 0012 Airfoil");
xlabel("Angle of Attack (degrees)");
ylabel("Coefficient of Lift");
legend("Experimental Data","Polyfit Line",'Location','northwest');
hold off;

%generate polyfit line for NACA2412 lift vs AoA
NACA2412_a0 = polyfit(NACA2412_cl_alpha(:,1),NACA2412_cl_alpha(:,2),4);
figure (12); %plot polyfit line against data to check compliance
hold on;
scatter(NACA2412_cl_alpha(:,1),NACA2412_cl_alpha(:,2),'o','b');
plot(alpha,polyval(NACA2412_a0,alpha),'r');
title("Coefficient of Lift vs. Angle of Attack for a NACA 2412 Airfoil");
xlabel("Angle of Attack (degrees)");
ylabel("Coefficient of Lift");
legend("Experimental Data","Polyfit Line",'Location','northwest');
hold off;

%generate polyfit line for c_d vs. c_l (drag polar) for NACA0012
NACA0012_polar = polyfit(NACA0012_cl_cd(:,1),NACA0012_cl_cd(:,2),4);
figure (13); %plot polyfit line against data to check compliance
hold on;
scatter(NACA0012_cl_cd(:,1),NACA0012_cl_cd(:,2),'o','b')
plot(coeff_lift1,polyval(NACA0012_polar,coeff_lift1),'r');
title("Drag Polar for a NACA 0012 Airfoil");
xlabel("Coefficient of Lift");
ylabel("Coefficient of Drag");
legend("Experimental Data","Polyfit Line",'Location','southwest');
hold off;

%generate polyfit line for c_d vs. c_l (drag polar) for NACA2412
NACA2412_polar = polyfit(NACA2412_cl_cd(:,1),NACA2412_cl_cd(:,2),5);
figure (14); %plot polyfit line against data to check compliance
hold on;
scatter(NACA2412_cl_cd(:,1),NACA2412_cl_cd(:,2),'o','b');
plot(coeff_lift2,polyval(NACA2412_polar,coeff_lift2),'r');
title("Drag Polar for a NACA 2412 Airfoil");
xlabel("Coefficient of Lift");
ylabel("Coefficient of Drag");
legend("Experimental Data","Polyfit Line",'Location','northwest');
hold off;

%generate c_d as a function of angle of attack by evaluating polyfit lines
%for lift slope and drag polar
NACA0012_cd = polyval(NACA0012_polar,polyval(NACA0012_a0,alpha));
%root airfoil has +1 degree geometric twist, hence alpha+1
NACA2412_cd = polyval(NACA2412_polar,polyval(NACA2412_a0,(alpha+1)));

%average root and tip c_d values to arrive at a sectional drag coefficient
%for the wing as a whole, as a function of angle of attack
wing_cd = mean([NACA0012_cd;NACA2412_cd]);

for i = 1:length(alpha) %loop through AoA values to run PLLT on each one
    %adjust geometric angles of attack to current values
    const_LD.geo_r = const_LD.wing_geo_r + alpha(i);
    const_LD.geo_t = const_LD.wing_geo_t + alpha(i);
    %find lift and induced drag coefficients for current AoA at N=50
    [~,wing_cl(i),wing_cdi(i)] = PLLT(const_LD,50);
end

%find L/D values (note c_l/c_d produces the same output)
%further note that C_d = c_d + c_di
aero_eff = wing_cl./(wing_cd+wing_cdi);

%display total drag coefficient for the given angle of attack 
disp("The total drag coefficient for the given wing is "+...
    string(wing_cd(27)+wing_cdi(27)));

%plot aerodynamic efficiency vs angle of attack
figure(15);
hold on;
plot(alpha,aero_eff);
title("Aerodynamic Efficiency L/D");
xlabel("Angle of Attack (degrees)");
ylabel("$\frac{L}{D}$",'Interpreter','latex','Rotation',0);
hold off;

disp("Elapsed Time = " + string(toc));

%% functions

function [err_index] = errorfinder(const,range,target,value)
% ERRORFINDER finds the number of N-values needed to converge to a given
% value of relative error for span efficiency, lift or drag, using a
% bifurcation search technique on the output of the PLLT function
%
%INPUTS:    const:  a STRUCT of constants for PLLT and other values (see
%                   getConst function
%           range:  a 2-element vector defining the guessed range for
%                   convergence
%           target: a targeted relative error value to converge to, as a
%                   percent
%           value:  a flag indicating which variable should be evaluated:
%           1=span efficiency, 2=lift, 3=drag
%OUTPUTS:   err_index:  the value of N needed to converge to within a
%                       relative error of *target*, evaluated against the
%                       value for N=1000
    if value == 1
        while diff(range) > 1   %looking to narrow the range down to two values
            index = floor((range(end)+range(1))/2); %value of range center
            [e,~,~] = PLLT(const,index); %PLLT values at this index
            e_error = abs((e - const.e_abs)/const.e_abs)*100;    %lift error at index
            if e_error < target
                range = [range(1) index];
            elseif e_error > target
                range = [index range(end)];
            elseif e_error == target
                range = [(index -1) index];
            end
        end
        err_index = range(2);
    elseif value == 2
        while diff(range) > 1   %looking to narrow the range down to two values
            index = floor((range(end)+range(1))/2); %value of range center
            [~,c_l_err,~] = PLLT(const,index); %PLLT values at this index
            lift_err = c_l_err*const.q*const.S;     %total lift at this index
            lift_error = abs((lift_err - const.lift_abs)/const.lift_abs)*100;    %lift error at index
            if lift_error < target
                range = [range(1) index];
            elseif lift_error > target
                range = [index range(end)];
            elseif lift_error == target
                range = [(index -1) index];
            end
        end
        err_index = range(2);

    elseif value == 3
        while diff(range) > 1   %looking to narrow the range down to two values
            index = floor((range(end)+range(1))/2); %value of range center
            [~,~,c_di_err] = PLLT(const,index); %PLLT values at this index
            drag_err = c_di_err*const.q*const.S;     %total lift at this index
            drag_error = abs((drag_err - const.drag_i_abs)/const.drag_i_abs)*100;    %lift error at index
            if drag_error < target
                range = [range(1) index];
            elseif drag_error > target
                range = [index range(end)];
            elseif drag_error == target
                range = [(index -1) index];
            end
        end
        err_index = range(2);
    else
        err_index = range(end);
    end
end

function const = getConst()
%getConst defines a struct of values describing the physical situation of a
%finite wing.  Geometry of the wing as well as some aspects of the flight
%condition are specified
%INPUTS:        None
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere 

    constPLLT = getConstPLLT();         %constants for PLLT
    %use PLLT to define constants for convergence, as well as some others
    constSYS = getConstSYS(constPLLT);
    %these lines concatenate the two preceding structs into one
    fn1 = fieldnames(constPLLT);
    fn2 = fieldnames(constSYS);
    fn = [fn1; fn2];
    c1 = struct2cell(constPLLT);
    c2 = struct2cell(constSYS);
    c = [c1;c2];
    const = cell2struct(c,fn,1);

function const = getConstPLLT()
%getConstPLLT is a subfunction that defines a struct of constants used by
%the PLLT function, in the context of Question 5 of the ASEN 3802 lab 3
%report.  Thus the wing geometry defined herein is that required by that
%section.
%INPUTS:        None
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere 

    % geometry of the wing
    const.b = 33 + (4/12);    %wingspan in feet
    const.c_r = 5 + (4/12);   %root chord in feet
    const.c_t = 3 +(8.5/12);  %tip chord in feet
    const.S = (const.b/2)*(const.c_r+const.c_t);    %planform area
    
    %root chord airfoil is a NACA 2412
    const.a0_r = 6.7789;    %lift slope in rad^-1
    const.aero_r = -2.2972;    %zero-lift angle of attack in degrees
    const.wing_geo_r = 1;      %geometric angle of attack in degrees
    %tip chord airfoil is a NACA 0012
    const.a0_t = 6.7786;    %lift slope in rad^-1
    const.aero_t = 0;    %zero-lift angle of attack in degrees
    const.wing_geo_t = 0;      %geometric angle of attack in degrees

    const.aot = 4;    %angle of attack the wing is flying in degrees
    %add this angle to the geo angles of the wing
    const.geo_r = const.wing_geo_r + const.aot;
    const.geo_t = const.wing_geo_t + const.aot;
end

function const = getConstSYS(constPLLT)
%getConstSYS is a subfunction that defines a struct of constants used by
%the errorfinder function, in the context of Question 5 of the ASEN 3802 lab 3
%report.  Thus the wing geometry defined herein is that required by that
%section.  This subfunction is separate from getConstPLLT because
%getConstSYS must call the PLLT function in order to define convergence
%values.
%INPUTS:        constPLLT:  A struct of constants needed by the PLLT
%                           function
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere
 
    % values of the physical system
    const.vel = 100 *1.68781;   %velocity of the plane in ft/sec
    const.altitude = 10000;       %altitude in feet
    
    %find atmospheric values
    %convert altitude to meters to use atmoscoesa function
    const.altitude_m = const.altitude * 0.3048;     %altitude in meters
    %find pressure and density at given altitude
    [~,~,const.P,const.rho] = atmoscoesa(const.altitude_m); 
    const.P = const.P* 0.0208854342; %convert Pascals to lb/ft^2
    const.rho = const.rho * 0.00194032; %convert kg/m^3 to slug/ft^3
   
    %find lift and drag coefficients with a very large N
    %assume these are the convergence values
    [const.e_abs,const.c_l_abs,const.c_di_abs] = PLLT(constPLLT,1000);
    
    %calculate convergence values for lift and induced drag using coefficients
    const.q = 0.5*const.rho*const.vel^2;  %dynamic pressure in Pa
    const.lift_abs = const.c_l_abs*const.q*constPLLT.S;     %total lift
    const.drag_i_abs = const.c_di_abs*const.q*constPLLT.S;  %total induced drag
end

end

function [e,c_l,c_di] = PLLT(const,N)
%PLLT applies the Fundamental Equation of Prandtl's Lifting Line Theory to
%a given airfoil, defined by input parameters
%   INPUTS: const:  A STRUCT containing the following fields:
%
%               const.b:      scalar value of wingspan (unit distance)
%               (degrees per unit distance)
%               const.a0_t:   scalar cross-sectional lift slope at wingtip
%               const.a0_r:   scalar cross-sectional lift slope at wing root
%               (unit distance)
%               const.c_t:    scalar value of chord at wingtip
%               const.c_r:    scalar value of chord at wing root
%               (degrees)
%               const.aero_t: scalar value of zero-lift angle of attack at
%               wingtip
%               const.aero_r: scalar value of zero-lift angle of attack at wing
%               root
%               (degrees)
%               const.geo_t:  scalar value of geometric angle of attack at the
%               wing tip
%               const.geo_r:  scalar value of geometric angle of attack at the
%               wing root
%           N:      scalar number of ODD terms to use in the Fourier
%               series expansion within PLLT
%
%   OUTPUTS:    e:      scalar value of span efficiency factor
%               c_l:    scalar value of coefficient of lift for the given
%               airfoil geometry
%               c_di:   scalar value of induced drag coefficient for the
%               given airfoil geometry

Nvec = linspace(1, (2*N)-1, N);     %linearly spaced vector of odd values of N
theta = (Nvec*pi)./(2*Nvec(end));       %corresponding angular locations along span
AR = (const.b^2)/const.S;      %calculate aspect ratio

%convert input values to radians
%a0_t = a0_t*(pi/180); a0_r = a0_r*(pi/180);
aero_t = const.aero_t*(pi/180); aero_r = const.aero_r*(pi/180);
geo_t = const.geo_t*(pi/180); geo_r = const.geo_r*(pi/180);

%get vectors of values at theta locations assuming linear variation along span
a0_vec = span(const.a0_t,const.a0_r,theta); 
c_vec = span(const.c_t,const.c_r,theta);
aero_vec = span(aero_t,aero_r,theta);
geo_vec = span(geo_t,geo_r,theta);

%generate matrix to solve for Fourier series coefficients (A_n)
%this matrix  is being generated as the sum of two matrices, one for each
%term in the Fundamental Equation of PLLT

%generate a matrix of size NxN of terms of the form n*theta
ntheta = ones(length(Nvec),length(theta));
ntheta = theta'.*ntheta;
ntheta = Nvec.*ntheta;

%generate a matrix of size NxN of terms of the form (4b/(a_0*c))sin(ntheta)
coeff_vec = (4*const.b)./(a0_vec.*c_vec); %leading term of Fund. Eq. of PLLT
mat1 = coeff_vec' .* sin(ntheta); %first coefficient matrix

%generate a matrix of size NxN of terms of the form
%n(sin(ntheta)/sin(theta))
mat2 = Nvec.*sin(ntheta);
mat2 = (1./sin(theta))' .* mat2;

%add these matrices and solve the system
finalmat = mat1 + mat2;
A_nvec = finalmat\(geo_vec - aero_vec)';

%calculate outputs using the Fourier series coefficients
e = 1/(1+sum(((A_nvec(2:end)./A_nvec(1)).^2).*Nvec(2:end)'));
c_l = A_nvec(1)*pi*AR;
c_di = (c_l^2)/(pi*e*AR);

function [y] = span(t, r, theta)
    % SPAN outputs an array of values at y locations, assuming linear
    % variation, decreasing from r at pi/2 to t at zero.
    %   INPUTS: t:      scalar value at theta=0
    %           r:      scalar value at theta=pi/2
    %           theta:  vector of angular locations between 0 and pi/2 to get values
    %                   assumed to be between zero and pi/2, increasing
    %   OUTPUTS:    span:   vector of values at y locations
    y = (t-r).*cos(theta) + r;  %calculate values
    
end

end