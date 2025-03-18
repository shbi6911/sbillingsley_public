% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     2-11-2025

%% Problem 1
%find Earth-Moon equilibrium points, with Jacobi constants, examine the
%effect of mass ratio on equilibrium point position
clear; clc; 
const = getConst(1);

%get equilibrium points for the Earth-Moon system
eq_pts_em = eqPts(const.mu.em);

%find Jacobi constants of equilibrium points in Earth-Moon
eq_pts_jc = jacobi(eq_pts_em,const.mu.em);

%examine the effect of mu on position of equilibrium points
mu_values = 0.5*logspace(-10,0); %vector of possible mu values
eq_pts_var = zeros(5,6,length(mu_values));  %preallocate 3d array of results
for i = 1:length(mu_values)
    eq_pts_var(:,:,i) = eqPts(mu_values(i));
end
L1_xcoords = squeeze(eq_pts_var(1,1,:));
L2_xcoords = squeeze(eq_pts_var(2,1,:));
L3_xcoords = squeeze(eq_pts_var(3,1,:));
L45_dist = sqrt((0.5 - mu_values).^2 + (3/4));
if const.vrb
    figure();
    semilogx(mu_values,L1_xcoords,'DisplayName',"L1");
    hold on;
    semilogx(mu_values,L2_xcoords,'DisplayName',"L2");
    semilogx(mu_values,L3_xcoords,'DisplayName',"L3");
    semilogx(mu_values,(0.5 - mu_values),'DisplayName',"L4&5");
    xlabel("$\mu$, Log Scale",'interpreter','latex');
    ylabel("Distance (Nondimensional Rotating Frame)");
    title("X-Distance in Rotating Frame")
    legend('Location','eastoutside');

    figure();
    semilogx(mu_values,L1_xcoords+mu_values','DisplayName',"L1");
    hold on;
    semilogx(mu_values,L2_xcoords+mu_values','DisplayName',"L2");
    semilogx(mu_values,-L3_xcoords-mu_values','DisplayName',"L3");
    semilogx(mu_values,L45_dist,'DisplayName','L4&5');
    xlabel("$\mu$, Log Scale",'interpreter','latex');
    ylabel("Distance (Nondimensional Rotating Frame)");
    title("Distance from P1");
    legend('Location','eastoutside');

    figure();
    semilogx(mu_values,(1-mu_values')-L1_xcoords,'DisplayName',"L1");
    hold on;
    semilogx(mu_values,L2_xcoords-(1-mu_values'),'DisplayName',"L2");
    semilogx(mu_values,-L3_xcoords+(1-mu_values'),'DisplayName',"L3");
    semilogx(mu_values,L45_dist,'DisplayName','L4&5');
    xlabel("$\mu$, Log Scale",'interpreter','latex');
    ylabel("X-position (Nondimensional Rotating Frame)");
    title("Distance from P2");
    ylim([-1 3]);
    legend('Location','eastoutside');
end

%% Problem 2
%find in-plane and out-of plane eigenvalues for Sun-Earth and Earth-Moon
%systems 
clear; clc; 
const = getConst(0);

%find Earth-Moon stability matrices
[inplane_em,outplane_em] = eqStab(const.mu.em,const.vrb);
%find Sun-Earth stability matrices
[inplane_se,outplane_se] = eqStab(const.mu.se,const.vrb);
eigval_se = zeros(5,6); eigval_em = zeros(5,6); %preallocate results
%find six eigenvalues for all five equilibrium points for both systems
for i = 1:5
    eigval_se(i,1:4) = eig(inplane_se(:,:,i));
    eigval_se(i,5:6) = eig(outplane_se(:,:,i));
    eigval_em(i,1:4) = eig(inplane_em(:,:,i));
    eigval_em(i,5:6) = eig(outplane_em(:,:,i));
end
%save certain eigenvalues/vectors for use in later sections
[V1_em,D1_em] = eig(inplane_em(:,:,1));
save('L1_eigvals.mat','V1_em','D1_em');

[V4_em,D4_em] = eig(inplane_em(:,:,4));
save('L4_eigvals.mat','V4_em','D4_em');


%% Problem 3
clear; clc; 
const = getConst(1);

mass_ratio = (25+3*sqrt(69))/2;
crit_val = 1/(1+mass_ratio);
mu_values = linspace(crit_val - crit_val*0.1,crit_val + crit_val*0.1); %vector of possible mu values
L4_eigs = zeros(6,length(mu_values));   %preallocate results
L5_eigs = zeros(6,length(mu_values));   % 6xn for 6 eigs and n mu values
%loop over mu values and extract eigenvalues for each one
for i = 1:length(mu_values)
    [inplane,outplane] = eqStab(mu_values(i));
    L4_eigs(1:4,i) = eig(inplane(:,:,4));
    L4_eigs(5:6,i) = eig(outplane(:,:,4));
    L5_eigs(1:4,i) = eig(inplane(:,:,5));
    L5_eigs(5:6,i) = eig(outplane(:,:,5));
end

if const.vrb
    figure();
    for i=1:6
        plot(mu_values, real(L4_eigs(i,:)));
        xline(crit_val)
        hold on;
    end
    
end
%% Problem 4
clear; clc; 
const = getConst(0);

%define inverse matrix for mathematical derivation
if const.vrb
    syms L1 L3 a1 a3
    syms xi eta xi_dot eta_dot
    B = [xi;xi_dot;eta;eta_dot];
    A = [1 1 1 1; L1 -L1 L3 -L3; a1 -a1 a3 -a3; a1*L1 a1*L1 a3*L3 a3*L3];
    disp(inv(A)*B);
end

%set initial conditions for position
xi_0 = -0.001;      eta_0 = 0;
%import values for Earth-Moon L1 from previous calculations
lambda_1 = 2.932055918598865;   lambda_3 = 2.334385875607174*1i;
L1_em_x = 0.836915131750382;    Uxx_L1_em = 11.295188987112699;
%calculate initial velocity based on given initial position
xi_dot_0 = (2*lambda_3^2*eta_0)/(lambda_3^2 - Uxx_L1_em);
eta_dot_0 = (xi_0/2)*(lambda_3^2 - Uxx_L1_em);
%calculate coefficients for ellipse
s = imag(lambda_3);     beta_3 = (s.^2 + Uxx_L1_em)/(2*s);
%generate sample time points to plot oscillatory solution
time_points = linspace(0,(2*pi)/s,100);
%calculate locations of one ellipse in terms of variations
xi = xi_0.*cos(s.*time_points) + (eta_0/beta_3).*sin(s.*time_points);
eta = eta_0.*cos(s.*time_points) - (beta_3*xi_0).*sin(s.*time_points);

%define initial state in the full rotating frame
initial_state = [L1_em_x + xi_0;0+eta_0;0;0+xi_dot_0;0+eta_dot_0;0];
tspan = [0,(2*pi)/s];   %nondimensional integration time
%integrate trajectory in the nonlinear system
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
[t, trajectory] = ode89(EOMfun, tspan, initial_state, opts);

const.vrb = 0;%plot oscillatory solution (linear & nonlinear)
if const.vrb
    figure(); hold on;
    plot(L1_em_x + xi,eta,'b','DisplayName','Linear Trajectory');
    plot(trajectory(:,1),trajectory(:,2),'r','DisplayName','NonLinear Trajectory');
    plot(L1_em_x,0,'Marker','.','MarkerSize',20,'Color','k','DisplayName','$L_1$');
    axis equal
    legend('interpreter','latex')
    title("Linear and Nonlinear Trajectory About Earth-Moon $L_1$",'interpreter','latex')
    xlabel("Nondimensional X-Position");
    ylabel("Nondimensional Y-Position");
end

%RUN SECTION 2 BEFORE THIS SECTION
load('L1_eigvals.mat');
eig_vector_3 = V1_em(:,3);      eig_vector_4 = V1_em(:,4);
basis1 = round(real(eig_vector_3),15);
basis2 = round(imag(eig_vector_3),15);


%% Problem 5
%note for future me; this code was written in a big hurry
clear; clc; 
const = getConst(1);

%RUN SECTION 2 BEFORE THIS SECTION
load('L4_eigvals.mat');
eig_values = diag(D4_em);       %extract eigenvalues
%parse out eigenvectors
sp_vector_1 = V4_em(:,1);   sp_vector_2 = V4_em(:,2);
lp_vector_1 = V4_em(:,3);   lp_vector_2 = V4_em(:,4);
%create basis vectors using real and imaginary parts
sp_basis_1 = real(sp_vector_1);     sp_basis_2 = imag(sp_vector_1);
lp_basis_1 = real(lp_vector_1);     lp_basis_2 = imag(lp_vector_2);
%generate scaling factor for a target position magnitude of 0.02
sp_factor_1 = 0.02/norm(sp_basis_1(1:2));
lp_factor_1 = 0.02/norm(lp_basis_1(1:2));
sp_factor_2 = 0.02/norm(sp_basis_2(1:2));
lp_factor_2 = 0.02/norm(lp_basis_2(1:2));
%generate initial conditions for propagation
sp_initial_state = round(sp_factor_1.*sp_basis_1,15);
lp_initial_state = round(lp_factor_1.*lp_basis_1,15);
%generate periods
sp_period = (2*pi)/imag(eig_values(1));
lp_period = (2*pi)/imag(eig_values(3));

%generate A2D matrix for Earth-Moon L4 (copy-pasting from eqStab function)
eq_pts = eqPts(const.mu.em);    mu = const.mu.em;
r1 = vecnorm([eq_pts(:,1)+mu,eq_pts(:,2),eq_pts(:,3)],2,2);
r2 = vecnorm([(eq_pts(:,1) - 1 + mu),eq_pts(:,2),eq_pts(:,3)],2,2);
Uxx = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*((eq_pts(:,1)+mu).^2))./r1.^5 +...
    (3.*mu.*(eq_pts(:,1) -1 + mu).^2)./r2.^5;
Uyy = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*(eq_pts(:,2).^2))./r1.^5 +...
    (3.*mu.*eq_pts(:,2).^2)./r2.^5;
Uxy = (3.*(1-mu).*(eq_pts(:,1)+mu).*eq_pts(:,2))./r1.^5 + ...
    (3.*mu.*(eq_pts(:,1) -1 + mu).*eq_pts(:,2))./r2.^5;
inplane = zeros(4,4);
inplane(1,3) = 1;    inplane(2,4) = 1;
inplane(3,4) = 2;    inplane(4,3) = -2;
inplane(3,1) = Uxx(4);  inplane(4,2) = Uyy(4);
inplane(3,2) = Uxy(4);  inplane(4,1) = Uxy(4);

%integrate linear trajectories
A2D_L4 = inplane;           %set A matrix
linearEOM = @(t,state) A2D_L4*state;    %define anonymous integration function
tspan_sp = [0,sp_period];
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
%integrate linear trajectory for short period
[~,linear_trajectory_sp] = ode89(linearEOM,tspan_sp,sp_initial_state,opts);
%integrate linear trajectory for long period
tspan_lp = [0,lp_period];
[~,linear_trajectory_lp] = ode89(linearEOM,tspan_lp,lp_initial_state,opts);

%generate nonlinear trajectories
non_sp_initial_state = [sp_initial_state(1:2); 0; sp_initial_state(3:4);0] + eq_pts(4,:)'; 
non_lp_initial_state = [lp_initial_state(1:2); 0; lp_initial_state(3:4);0] + eq_pts(4,:)';
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);    %anon function for nonlinear
%integrate nonlinear trajectory for short period
[~,nonlinear_trajectory_sp] = ode89(EOMfun,tspan_sp,non_sp_initial_state,opts);
%integrate nonlinear trajectory for long period
[~,nonlinear_trajectory_lp] = ode89(EOMfun,tspan_lp,non_lp_initial_state,opts);

if const.vrb
    figure();   hold on; grid on;
    plot(linear_trajectory_sp(:,1) + eq_pts(4,1),linear_trajectory_sp(:,2) + eq_pts(4,2),...
        'b','DisplayName','Linear');
    plot(nonlinear_trajectory_sp(:,1),nonlinear_trajectory_sp(:,2),...
        'r','DisplayName','Nonlinear');
    plot(eq_pts(4,1),eq_pts(4,2),'Marker','.','MarkerSize',20,'Color','k','DisplayName','$L_4$');
    plot(non_sp_initial_state(1),non_sp_initial_state(2),'Marker','.','MarkerSize',20,'Color','m','DisplayName','Initial Position');
    legend('interpreter','latex');
    title("Short Period Mode");
    axis equal

    figure();   hold on; grid on;
    plot(linear_trajectory_lp(:,1) + eq_pts(4,1),linear_trajectory_lp(:,2) + eq_pts(4,2),...
        'b','DisplayName','Linear');
    plot(nonlinear_trajectory_lp(:,1),nonlinear_trajectory_lp(:,2),...
        'r','DisplayName','Nonlinear');
    plot(eq_pts(4,1),eq_pts(4,2),'Marker','.','MarkerSize',20,'Color','k','DisplayName','$L_4$');
    plot(non_lp_initial_state(1),non_lp_initial_state(2),'Marker','.','MarkerSize',20,'Color','m','DisplayName','Initial Position');
    legend('interpreter','latex');
    title("Long Period Mode");
    axis equal
end
    

%% Functions

function eq_pts = eqPts(mu,vrb,guess)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     2-11-2025
%
% eqPts is a function to determine the equilibrium points of a system under
% the CR3BP, with the system defined by a value for mass ratio
%
%INPUTS:    mu      required input value of mass ratio
%           vrb     optional input for verbosity, controls output
%           guess   optional initial guesses for L1,L2,L3, without input
%                   these are determined programmatically
%OUTPUTS:   eq_pts  5x6 matrix of state vectors at each of the equilibrium
%                   points, in order from top to bottom (L1;L2;L3;L4;L5)
    arguments
        mu double
        vrb logical =0
        guess (1,3) =[0,0,0]
    end
    %define function for collinear equilibrium point location
    fun = @(x,mu) x - ((1 - mu)*(x + mu))/((abs(x+mu))^3) - ...
        (mu*(x - 1 + mu))/((abs(x - 1 + mu))^3);
    %preallocate output matrix
    eq_pts = zeros(5,6);
    %define initial guesses
    if norm(guess) == 0
        guess = [(1-(mu/3)^(1/3));(1+(mu/3)^(1/3));...
                 (-(1 +(5/12)*mu))];
    end
    %set options for rootfinding function
    if vrb
        opts = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1e-16);
    else
        opts = optimoptions(@fsolve,'Display','none','FunctionTolerance',1e-16);
    end
    
    %find all three collinear equilibrium points using fsolve
    for i = 1:length(guess)
        eq_pts(i,1) = fsolve(@(x)fun(x,mu),guess(i),opts);
    end
    % calculate the locations of L4 and L5 using analytical solutions
    eq_pts(4,1:2) = [0.5 - mu,sqrt(3)/2];
    eq_pts(5,1:2) = [0.5 - mu,-sqrt(3)/2];
end

function [inplane,outplane] = eqStab(mu,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     2-11-2025
%
% eqStab is a function to output the linearized stability matrices for the
% five equilibrium points of a CR3BP system
%DEPENDENCY:        this function calls the function eqPts
%
%INPUTS:    mu      required input value of mass ratio
%           vrb     optional input for verbosity, controls output
%OUTPUTS:   inplane     a THREE dimensional matrix of the size 4,4,5, where
%                   each of the five pages is a 4x4 linearized stability
%                   matrix for the in-plane stability modes of an
%                   equilibrium point, L1....L5
%           outplane    a THREE dimensional matrix of the size 2,2,5 where
%                   each page is a 2x2 stability matrix for the out of
%                   plane modes, for L1....L5
    arguments
        mu double
        vrb logical =0
    end
    eq_pts = eqPts(mu,vrb);     %locate equilibrium points
    inplane = zeros(4,4,5);     outplane = zeros(2,2,5);    %preallocate
    %calculate 5x1 vector of distances from equilibrium pts to P1 and P2
    r1 = vecnorm([eq_pts(:,1)+mu,eq_pts(:,2),eq_pts(:,3)],2,2);
    r2 = vecnorm([(eq_pts(:,1) - 1 + mu),eq_pts(:,2),eq_pts(:,3)],2,2);
    %calculate 5x1 vector of partial derivatives for each equilibrium pt
    Uzz = -(1-mu)./r1.^3 - mu./r2.^3 + (3.*(1 - mu).*eq_pts(:,3).^2)./r1.^5 +...
        (3.*mu.*eq_pts(:,3).^2)./r2.^5;
    Uxx = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*((eq_pts(:,1)+mu).^2))./r1.^5 +...
        (3.*mu.*(eq_pts(:,1) -1 + mu).^2)./r2.^5;
    Uyy = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*(eq_pts(:,2).^2))./r1.^5 +...
        (3.*mu.*eq_pts(:,2).^2)./r2.^5;
    Uxy = (3.*(1-mu).*(eq_pts(:,1)+mu).*eq_pts(:,2))./r1.^5 + ...
        (3.*mu.*(eq_pts(:,1) -1 + mu).*eq_pts(:,2))./r2.^5;
    %build out out-of-plane matrices
    for ii = 1:size(outplane,3)
        outplane(1,2,ii) = 1;
        outplane(2,1,ii) = Uzz(ii);
    end
    %build in-plane matrices
    for jj = 1:size(inplane,3)
        inplane(1,3,jj) = 1;    inplane(2,4,jj) = 1;
        inplane(3,4,jj) = 2;    inplane(4,3,jj) = -2;
        inplane(3,1,jj) = Uxx(jj);  inplane(4,2,jj) = Uyy(jj);
        inplane(3,2,jj) = Uxy(jj);  inplane(4,1,jj) = Uxy(jj);
    end
end