function [t,state_matrix,norms,ier] = corrArcLength(V,mu,n_star,vrb,constraint,step,Nmax,tol)
%correction implements a single-shooting correction scheme to converge an
%initial guess for a periodic orbit in the CR3BP toward a true periodic orbit.  
% By default, it implements the "general variable time formulation" for the 
% scheme.  An optional input allows the use of a "modified constraint" 
% formulation for each of the three velocity components.
%
% DEPENDENCY    This function calls threeBP_EOM and three_BPrefTraj
%               This function also calls correction
%
% INPUTS:   V          a 7X1 initial guess for a periodic orbit, with the
%                       first six elements being a state vector and the final 
%                       element being orbit period
%           mu          mass ratio of the system of interest
%           n_star      direction vector for pseudo-arclength continuation
%           vrb         optional verbosity input to control output
%           constraint  optional flag for a modified constraint
%                       0 = none, 1=x_dot, 2=y_dot, 3=z_dot
%           step        step size for pseudo-arclength continuation
%           Nmax        maximum number of iterations toward periodic orbit
%           tol         tolerance on the norm of constraint vector
%
% OUTPUTS:  t           integrated vector of time points for the final orbit
%           state_matrix    integrated state over time for the final orbit
%           norms       1xn vector of constraint vector norm over iteration
%           ier         success flag 0=convergence 1=failure
    arguments
        V(7,1) {mustBeNumeric}
        mu double
        n_star (7,1) {mustBeNumeric}
        vrb logical =0
        constraint {mustBeMember(constraint,[0,1,2,3])} =0
        step double =0
        Nmax double =100
        tol double =1e-12
    end
    %set initial conditions based on input
    initial_state = [V(1:6);reshape(eye(6),[],1)];    tspan = [0 V(7)];
    V_star = V;         %set fixed V-star vector
    %set tolerances and integrate
    opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
    EOMfun = @(t,r)threeBP_refTraj(t,r,mu);
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    %define constraint vector F
    F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
    %define Jacobian
    jacob = zeros(length(F),length(V)); %preallocate
    % extract state matrix
    phi = reshape(state_matrix(end,7:42),6,6);
    jacob(1:6,1:6) = phi - eye(6);
    jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',mu);
    % redefine constraint vector and Jacobian for modified-constraint
    % based on input flag
    if constraint == 1              %constrain x0 = 0
        F(4) = [];  F(6) = V(1);
        jacob(4,:) = [];    jacob(6,:) = [1,0,0,0,0,0,0];
    elseif constraint == 2          %constrain y0 = 0
        F(5) = [];  F(6) = V(2);
        jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
    elseif constraint == 3          %constrain z0 = 0
        F(6) = [];  F(6) = V(3);
        jacob(6,:) = [];    jacob(6,:) = [0,0,1,0,0,0,0];
    end
    %set final constraint for pseudo-arclength continuation
    F(7) = dot((V - V_star),n_star) - step;
    jacob(7,:) = n_star';

    N = 1;  %initialize iteration counter
    if vrb
        figure();   %create figure for function debugging output
    end
    %preallocate and store constraint evolution
    norms = zeros(length(F),Nmax);     norms(:,N) = F;
    %implement Newton's method using minimum norm solution
    while (norm(F)>tol && N <= Nmax)
        sz = size(jacob);
        if sz(1) == sz(2)       %update free variable vector
            V = V - jacob\F;        %if square, use simple linear solve
        else
            V = V - jacob'*((jacob*jacob')\F);  %use minimum norm solution
        end
        %reintegrate trajectory with new free variables
        initial_state = [V(1:6);reshape(eye(6),[],1)];  %add initial STM
        tspan = [0 V(7)];
        [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
        %redefine constraint vector F
        F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
        %define Jacobian
        phi = reshape(state_matrix(end,7:42),6,6);  % extract state matrix
        jacob(1:6,1:6) = phi - eye(6);
        jacob(:,7) = [threeBP_EOM(t(end),state_matrix(end,1:6)',mu);0];
        % redefine constraint vector and Jacobian for modified-constraint
        % based on input flag
        if constraint == 1              %constrain x0 = 0
            F(4) = [];  F(6) = V(1);
            jacob(4,:) = [];    jacob(6,:) = [1,0,0,0,0,0,0];
        elseif constraint == 2          %constrain y0 = 0
            F(5) = [];  F(6) = V(2);
            jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
        elseif constraint == 3          %constrain z0 = 0
            F(6) = [];  F(6) = V(3);
            jacob(6,:) = [];    jacob(6,:) = [0,0,1,0,0,0,0];
        end
        %set final constraint for pseudo-arclength continuation
        F(7) = dot((V - V_star),n_star) - step;
        jacob(7,:) = n_star';

        %increment counter and store constraint vector
        N = N+1;    norms(:,N) = F;
        %some debugging output
        if vrb
            if norm(state_matrix(:,3)) < 1e-10
                plot(state_matrix(:,1),state_matrix(:,2)); hold on; grid on; axis equal;
            else
                plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3)); hold on; grid on; axis equal;
            end

            if norm(F)<tol
                fprintf("Converged successfully to norm(F) = %d using %i iterations\n",norm(F),N);
            elseif  N >Nmax
                fprintf("Failed convergence, too many iterations\n");
            end %if norm(F)<tol or N>Nmax
        end %if vrb
    end %while (norm(F)>tol && N <= Nmax)
    %set output variables
    norms(:,vecnorm(norms) == 0) = [];     %truncate unused values
    if norm(F)<tol
        ier = 0;
    elseif  N >Nmax
        ier = 1;
    end
end %function