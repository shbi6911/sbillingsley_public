% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-04-2025

%% Problem 1
% write a script to numerically integrate a state vector and a state
% transition matrix along a trajectory in the CR3BP
%
clear; clc; const = getConst(1);
% set initial state and time span
initial_state = [0.5;0;0;0;0;0;reshape(eye(6),[],1)];
tspan = [0 2];
%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);

if const.vrb
    plot(state_matrix(:,1),state_matrix(:,2));
    circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r','');
end

%% Problem 2
%generate an initial guess for an L1 Lyapunov orbit, correct to a periodic
%orbit and plot (use modified constraint corrections algorithm)

clear; clc; const = getConst(1);
%generate an initial guess for an L1 Lyapunov orbit
%find Earth-Moon stability matrices
[inplane_em,~] = eqStab(const.mu.em);
%find oscillatory eigenvector for L1
[V1_em,D1_em] = eig(inplane_em(:,:,1));
% create basis for eigenspace of lambda3
basis1 = round(real(V1_em(:,3)),15);
basis2 = round(imag(V1_em(:,3)),15);
% generate initial state for integration @ L1
eq_pts = eqPts(const.mu.em);        %find equilibrium points
L1_initial = eq_pts(1,:)';
% add initial guess with initial STM (eye(6))
initial_state = [L1_initial;reshape(eye(6),[],1)];
initial_state(1:2) = initial_state(1:2) + basis1(1:2).*(1/250);
initial_state(4:5) = initial_state(4:5) + basis1(3:4).*(1/250);
tspan = [0 (2*pi)/imag(D1_em(3,3))];    %time span using period from eigenvalue
%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);

% plot resulting L1 Lyapunov orbit initial guess
if const.vrb
    figure(1);
    plot(state_matrix(:,1),state_matrix(:,2),'b','DisplayName',"Initial Guess"); hold on;
    plot(eq_pts(1,1),0,'.','MarkerSize',10,'Color','k','DisplayName',"$L_1$");
    axis equal
end

%define initial free variable vector V
V = [state_matrix(1,1:6)';t(end)];
%define constraint vector F
F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
%redefine F for modified-constraint                 MODIFIED CONSTRAINT
F(5) = [];  F(6) = state_matrix(1,2);
%define Jacobian
jacob = zeros(length(F),length(V)); %preallocate
% extract state matrix
phi = reshape(state_matrix(end,7:42),6,6);
jacob(1:6,1:6) = phi - eye(6);
jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
%redefine Jacobian for modified-constraint          MODIFIED CONSTRAINT
jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
% define rootfinding variables
Nmax = 100; N = 1;  tol = 1e-12;
%preallocate and store constraint evolution
constraint_evo = zeros(1,Nmax);     constraint_evo(N) = norm(F);
%implement Newton's method using minimum norm solution
while (norm(F)>tol && N <= Nmax)
    %V = V - lsqminnorm(jacob,F);    %define new free variable vector
    V = V - jacob'*((jacob*jacob')\F);
    %reintegrate trajectory with new free variables
    initial_state = [V(1:6);reshape(eye(6),[],1)];  %add initial STM
    tspan = [0 V(7)];
    opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    %define constraint vector F
    F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
    %redefine F for modified-constraint                 MODIFIED CONSTRAINT
    F(5) = [];  F(6) = state_matrix(1,2);
    %define Jacobian
    jacob = zeros(length(F),length(V)); %preallocate
    % extract state matrix
    phi = reshape(state_matrix(end,7:42),6,6);
    jacob(1:6,1:6) = phi - eye(6);
    jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
    %redefine Jacobian for modified-constraint          MODIFIED CONSTRAINT
    jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
    N = N+1;            %increment iteration counter
    constraint_evo(N) = norm(F);    %store constraint evolution
    %some debugging output
    if const.vrb
        figure(2);
        plot(state_matrix(:,1),state_matrix(:,2)); hold on; axis equal;
        if norm(F)<tol
            fprintf("Converged successfully to norm(F) = %d using %i iterations",norm(F),N);
        elseif  N >Nmax
            fprintf("Failed convergence, too many iterations");
        end
    end
end
%set final determination of Lyapunov initial state
constraint_evo(constraint_evo == 0) = [];   %truncate unused values
L1_IC = state_matrix(1,1:6);
tspan = [0 V(7)];
%reintegrate to verify periodic orbit
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
[t,state_matrix] = ode113(EOMfun, tspan, L1_IC, opts);

% plot final periodic Lyapunov orbit
if const.vrb
    figure(1); 
    plot(state_matrix(:,1),state_matrix(:,2),'r','DisplayName',"Final Orbit");
    %plot(eq_pts(1,1),0,'.','MarkerSize',20,'Color','k');
    %plot(eq_pts(2,1),0,'.','MarkerSize',20,'Color','r');
    %circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
    title("$L_1$ Lyapunov Orbit Correction",'interpreter','latex');
    xlabel("X-Position (Nondimensional Units, Earth-Moon Frame");
    ylabel("Y-Position");
    legend('interpreter','latex');
end

%plot evolution of the constraint vector
if const.vrb
    figure();
    semilogy(constraint_evo);
    xlabel("Iterations");
    ylabel("Norm of Constraint Vector")
    title("Evolution of Constraint Vector over Iteration")
end

%% Problem 3
% calcluate and plot the family of L1 Lyapunov orbits using natural
% parameter continuation in x-position

%set continuation step size
delta_x = 0.001;



%% Problem 4
%use a given initial guess to correct to an L1 halo orbit, using a modified
%constraint formulation

clear; clc; const = getConst(1);
% given initial guess for an L1 halo orbit
initial_state = [0.82340;0;-0.026755;0;0.13742;0;reshape(eye(6),[],1)];
tspan = [0 2.7477];

%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);

% plot resulting L1 halo orbit
if const.vrb
    figure(1);
    plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3),...
        'b','DisplayName',"Initial Guess"); hold on;
    plot3(const.eq_pts.em(1,1),const.eq_pts.em(1,2),const.eq_pts.em(1,3),...
        '.','MarkerSize',10,'Color','k','DisplayName',"$L_1$"); 
end

%define initial free variable vector V
V = [state_matrix(1,1:6)';t(end)];
%define constraint vector F
F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
%redefine F for modified-constraint                 MODIFIED CONSTRAINT
F(5) = [];  F(6) = state_matrix(1,2);
%define Jacobian
jacob = zeros(length(F),length(V)); %preallocate
% extract state matrix
phi = reshape(state_matrix(end,7:42),6,6);
jacob(1:6,1:6) = phi - eye(6);
jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
%redefine Jacobian for modified-constraint          MODIFIED CONSTRAINT
jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
% define rootfinding variables
Nmax = 100; N = 1;  tol = 1e-12;
%preallocate and store constraint evolution
constraint_evo = zeros(1,Nmax);     constraint_evo(N) = norm(F);
%implement Newton's method using minimum norm solution
while (norm(F)>tol && N <= Nmax)
    %V = V - lsqminnorm(jacob,F);    %define new free variable vector
    V = V - jacob'*((jacob*jacob')\F);
    %reintegrate trajectory with new free variables
    initial_state = [V(1:6);reshape(eye(6),[],1)];  %add initial STM
    tspan = [0 V(7)];
    opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    %define constraint vector F
    F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
    %redefine F for modified-constraint                 MODIFIED CONSTRAINT
    F(5) = [];  F(6) = state_matrix(1,2);
    %define Jacobian
    jacob = zeros(length(F),length(V)); %preallocate
    % extract state matrix
    phi = reshape(state_matrix(end,7:42),6,6);
    jacob(1:6,1:6) = phi - eye(6);
    jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
    %redefine Jacobian for modified-constraint          MODIFIED CONSTRAINT
    jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
    N = N+1;            %increment iteration counter
    constraint_evo(N) = norm(F);    %store constraint evolution
    %some debugging output
    if const.vrb
        figure(2);
        %plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3));
        if norm(F)<tol
            fprintf("Converged successfully to norm(F) = %d using %i iterations",norm(F),N);
        elseif  N >Nmax
            fprintf("Failed convergence, too many iterations");
        end
    end
end
%set final determination of halo orbit initial state
constraint_evo(constraint_evo == 0) = [];   %truncate unused values
L1_IC = state_matrix(1,1:6);
tspan = [0 V(7)];
%reintegrate to verify periodic orbit
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
[t,state_matrix] = ode113(EOMfun, tspan, L1_IC, opts);

% plot final periodic halo orbit
if const.vrb
    figure(1); 
    plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3),...
        'r','DisplayName',"Final Halo Orbit");
    %plot3(const.eq_pts.em(1,1),const.eq_pts.em(1,2),const.eq_pts.em(1,3),...
    %    '.','MarkerSize',10,'Color','k');
    xlabel("X");    ylabel("Y");    zlabel("Z");
    hold on; grid on;   
    %ylim([min(state_matrix(:,2))*1.25 max(state_matrix(:,2))*1.25]);
    %xlim([min(state_matrix(:,1))*0.9 max(state_matrix(:,1))*1.1]);
    axis equal;
    legend('interpreter','latex')
end

%plot evolution of the constraint vector
if const.vrb
    figure();
    semilogy(constraint_evo);
    xlabel("Iterations");
    ylabel("Norm of Constraint Vector")
    title("Evolution of Constraint Vector over Iteration")
end

