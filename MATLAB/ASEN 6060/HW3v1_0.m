% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-04-2025

%% Problem 1
% write a script to numerically integrate a state vector and a state
% transition matrix along a trajectory in the CR3BP
%
clear; clc; const = getConst(0);
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

clear; clc; const = getConst(1);
%generate an initial guess for an L1 Lyapunov orbit
%find Earth-Moon stability matrices
[inplane_em,~] = eqStab(const.mu.em);
%find oscillatory eigenvector for L1
[V1_em,D1_em] = eig(inplane_em(:,:,1));
% create basis for eigenspace of lambda3
basis1 = round(real(V1_em(:,3)),15);
basis2 = round(imag(V1_em(:,3)),15);
% generate initial state for integration
eq_pts = eqPts(const.mu.em);        %find equilibrium points
L1_initial = eq_pts(1,:)';
% add initial guess with initial STM (eye(6))
initial_state = [L1_initial;reshape(eye(6),[],1)];
initial_state(1:2) = initial_state(1:2) + basis1(1:2).*(1/200);
initial_state(4:5) = initial_state(4:5) + basis1(3:4).*(1/200);
tspan = [0 (2*pi)/imag(D1_em(3,3))];
%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);

% plot resulting L1 Lyapunov orbit initial guess
if const.vrb
    figure(1);
    plot(state_matrix(:,1),state_matrix(:,2)); hold on;
    plot(eq_pts(1,1),0,'.','MarkerSize',10,'Color','k');
    axis equal
end

%define initial free variable vector V
V = [state_matrix(1,1:6)';t(end)];
%define constraint vector F
F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
%define Jacobian
jacob = zeros(length(F),length(V)); %preallocate
% extract state matrix
phi = reshape(state_matrix(end,7:42),6,6);
jacob(1:6,1:6) = phi - eye(6);
jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
% define rootfinding variables
Nmax = 100; N = 0;  tol = 1e-12;
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
    %define Jacobian
    jacob = zeros(length(F),length(V)); %preallocate
    % extract state matrix
    phi = reshape(state_matrix(end,7:42),6,6);
    jacob(1:6,1:6) = phi - eye(6);
    jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
    N = N+1;            %increment iteration counter
    %some debugging output
    if const.vrb
        figure(1);
        plot(state_matrix(:,1),state_matrix(:,2));
        if norm(F)<tol
            fprintf("Converged successfully to norm(F) = %d\n",norm(F));
        elseif  N >Nmax
            fprintf("Too many iterations");
        end
    end
end
%set final determination of Lyapunov initial state
L1_IC = state_matrix(1,1:6);
tspan = [0 V(7)];
%reintegrate to verify periodic orbit
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
[t,state_matrix] = ode113(EOMfun, tspan, L1_IC, opts);

% plot final periodic Lyapunov orbit
if const.vrb
    figure(); 
    plot(state_matrix(:,1),state_matrix(:,2)); hold on; axis equal
    plot(eq_pts(1,1),0,'.','MarkerSize',20,'Color','k');
    %plot(eq_pts(2,1),0,'.','MarkerSize',20,'Color','r');
    %circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','');
end

%% Problem 3

%% Problem 4

clear; clc; const = getConst(1);
% given initial guess for an L4 halo orbit
initial_state = [0.82340;0;-0.026755;0;0.13742;0;reshape(eye(6),[],1)];
tspan = [0 2.7477];

%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);

% plot resulting L4 halo orbit
if const.vrb
    figure(1);
    plot(state_matrix(:,1),state_matrix(:,2)); hold on;
    plot(eq_pts(1,1),0,'.','MarkerSize',10,'Color','k');
    axis equal
end

%% Functions

function [value, isterminal, direction] = eventFun(t,y)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-4-2025
%
%eventFun is an event function for use with threeBP_relTraj and an ODE solver.
% It stops integration at the first point where y = 0 and the function y is
% increasing
%
%INPUTS:    t,y     time and state from the ode solver
%
%OUTPUTS    value       the state at which y=0
%           isterminal  a value of 1 indicates stop integration
%           direction   a value of 1 indicates the function is increasing
    value = y(2);
    isterminal = 1;
    direction = 1;
end