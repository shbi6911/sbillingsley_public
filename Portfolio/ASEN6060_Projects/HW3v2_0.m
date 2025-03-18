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
V_init = [eq_pts(1,:)';(2*pi)/imag(D1_em(3,3))];
% add initial guess
scale_factor = 1/250;   %set scaling factor for initial guess
V_init(1:2) = V_init(1:2) + basis1(1:2).*scale_factor;
V_init(4:5) = V_init(4:5) + basis1(3:4).*scale_factor;
%V_init = [ 0.836796532867765;0;0;0;0.000992942246581;0;2.691579559675719];
tspan = [0 V_init(7)];    %time span using period from eigenvalue
%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, tspan, [V_init(1:6);reshape(eye(6),[],1)], opts);

% plot resulting L1 Lyapunov orbit initial guess
if const.vrb
    figure(1); 
    plot(state_matrix(:,1),state_matrix(:,2),'b','DisplayName',"Initial Guess"); hold on; grid on;
    plot(const.eq_pts.em(1,1),0,'.','MarkerSize',10,'Color','k','DisplayName',"$L_1$");
    axis equal
end
vrb = 0;        %set output flag for this segment
%implement modified constraint correction to recover periodic orbit
%modify constraint vector in y_dot, constrain y0 to x-axis
[t,state_matrix,norms,~] = correction(V_init,const.mu.em,vrb,2);

% plot final periodic Lyapunov orbit
if const.vrb
    figure(1); 
    plot(state_matrix(:,1),state_matrix(:,2),'r','DisplayName',"Final Orbit");
    title("$L_1$ Lyapunov Orbit Correction",'interpreter','latex');
    xlabel("X-Position (Nondimensional Units, Earth-Moon Frame");
    ylabel("Y-Position");
    legend('Location','southwest','interpreter','latex');
end

%plot evolution of the constraint vector
if const.vrb
    figure(); 
    semilogy(norms);    grid on;
    xlabel("Iterations");
    ylabel("Norm of Constraint Vector")
    title("Evolution of Constraint Vector over Iteration")
end

%% Problem 3
% calcluate and plot the family of L1 Lyapunov orbits using natural
% parameter continuation in x-position
tic;
clear; clc; const = getConst(1);
%set continuation step size
step = -0.001;
%set maximum number of steps
step_max = 500;
%initialize step counter
n_orb = 1;

%get initial orbit
%find Earth-Moon stability matrices
[inplane_em,~] = eqStab(const.mu.em);
%find oscillatory eigenvector for L1
[V1_em,D1_em] = eig(inplane_em(:,:,1));
% create basis for eigenspace of lambda3
basis1 = round(real(V1_em(:,3)),15);
basis2 = round(imag(V1_em(:,3)),15);
% generate initial state for integration @ L1
V_init = [const.eq_pts.em(1,:)';(2*pi)/imag(D1_em(3,3))];
% add initial guess
scale_factor = 1/10000;   %set scaling factor for initial guess
V_init(1:2) = V_init(1:2) + basis1(1:2).*scale_factor;
V_init(4:5) = V_init(4:5) + basis1(3:4).*scale_factor;
%correct initial guess using modified constraint in y0
[t,state_matrix,~,ier] = correction(V_init,const.mu.em,0,2);
%store results in struct
lyap_L1.("orbit" + string(n_orb)).time = t;
lyap_L1.("orbit" + string(n_orb)).state = state_matrix;

%plot initial orbit
if const.vrb
    figure(1);   hold on;    grid on;
    plot(state_matrix(:,1),state_matrix(:,2),'r');
    h1 = plot(const.eq_pts.em(1,1),0,'.','MarkerSize',10,'Color','r','DisplayName',"$L_1$");
    %plot(const.eq_pts.em(2,1),0,'.','MarkerSize',10,'Color','r','DisplayName',"$L_2$");
    h2 = circle(1-const.mu.em,0,const.radius.Moon/const.l_star.em,'k','Moon'); hold on;
    h3 = circle(-const.mu.em,0,const.radius.Earth/const.l_star.em,'r','Earth'); hold on;
end

%begin loop for natural parameter continuation
while (n_orb < step_max +1 && ier ~= 1)
    V_new = [state_matrix(1,1:6)';t(end)];     %new guess is previous initial state
    %reintegrate with added step
    [t,state_matrix,~,ier] = corrNatParam(V_new,const.mu.em,0,2,1,step);
    n_orb = n_orb + 1;  %update step counter
    % if mod(n_orb,100) == 0
    %     step = step*2;          %variable step size
    % end
    %store results in struct
    lyap_L1.("orbit" + string(n_orb)).time = t;
    lyap_L1.("orbit" + string(n_orb)).state = state_matrix;
end
%load("L1_lyapunov.mat")    
%prepare requested output
periods = zeros(numel(fieldnames(lyap_L1)),1);  %preallocate
states = zeros(numel(fieldnames(lyap_L1)),6);
%extract period and initial state vector from each trajectory
for i = 1:numel(fieldnames(lyap_L1))            
    periods(i) = lyap_L1.("orbit" + string(i)).time(end);
    states(i,:) = lyap_L1.("orbit" + string(i)).state(1,1:6);
end
jacobi_consts = jacobi(states,const.mu.em);     %calculate jacobi constants
periods_days = (periods.*const.t_star.em)./(3600*24);   %dimensionalize period

%plot period vs. jacobi constant
if const.vrb
    figure();
    plot(periods_days,jacobi_consts,'LineWidth',2);
    xlabel("Orbit Period in Days");
    ylabel("Jacobi Constant");
    title("Orbit Period vs. Jacobi Constant Along the $L_1$ Lyapunov Family",...
        'interpreter','latex');
end %if const.vrb

%plot output trajectories
if const.vrb
    figure(1); 
    for i = 1:10:numel(fieldnames(lyap_L1))
        X = lyap_L1.("orbit" + string(i)).state(:,1);
        Y = lyap_L1.("orbit" + string(i)).state(:,2);
        plot(X,Y,'b');
    end %plotting for loop
    axis equal;
    xlabel("X-Position, Nondimensional Length Units");
    ylabel("Y-Position");
    title("Orbits in the $L_1$ Lyapunov Family",'interpreter','latex');
    legend([h1,h2,h3],'Location','southwest','interpreter','latex');
end %if const.vrb

save("L1_Lyapunov.mat",'lyap_L1');
fprintf("Elapsed Time is %d minutes \n",toc/60);

%% Problem 4
%use a given initial guess to correct to an L1 halo orbit, using a modified
%constraint formulation

clear; clc; const = getConst(1);
% given initial guess for an L1 halo orbit
V_init = [0.82340;0;-0.026755;0;0.13742;0;2.7477];

%set tolerances and integrate
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
EOMfun = @(t,r)threeBP_refTraj(t,r,const.mu.em);
[t,state_matrix] = ode113(EOMfun, [0 V_init(7)], [V_init(1:6);reshape(eye(6),[],1)], opts);

% plot resulting L1 halo orbit
if const.vrb
    figure(1);
    plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3),...
        'b','DisplayName',"Initial Guess"); hold on; grid on;
    plot3(const.eq_pts.em(1,1),const.eq_pts.em(1,2),const.eq_pts.em(1,3),...
        '.','MarkerSize',10,'Color','k','DisplayName',"$L_1$"); 
end

% use modified constraint in y0 to correct to final orbit
vrb = 0;        %set independent output flag for this segment
[t,state_matrix,norms,~] = correction(V_init,const.mu.em,vrb,2);

% plot final periodic halo orbit
if const.vrb
    figure(1); 
    plot3(state_matrix(:,1),state_matrix(:,2),state_matrix(:,3),...
        'r','DisplayName',"Final Halo Orbit");
    xlabel("X-Position, Nondimensional Length Units");    ylabel("Y-Position");    
    zlabel("Z-Position");
    hold on; grid on;   
    axis equal;
    legend('interpreter','latex')
end

%plot evolution of the constraint vector
if const.vrb
    figure(); 
    semilogy(norms);        grid on;
    xlabel("Iterations");
    ylabel("Norm of Constraint Vector")
    title("Evolution of Constraint Vector over Iteration")
end

%% Problem 5
% using pseudo-arclength continuation, continue the L1 halo orbit family
% from the final orbit calculated in problem 4, for as long as possible
tic;
clear; clc; const = getConst(1);
%set continuation step size
step = 0.001;
%set maximum number of steps
step_max = 1800;
%initialize step counter
n_orb = 1;
%set modified-constraint flag
constraint = 2;

%recover corrected orbit for Problem 4, starting orbit for continuation
% given initial guess for an L1 halo orbit
V_init = [0.82340;0;-0.026755;0;0.13742;0;2.7477];
%correct initial guess using modified constraint in y0
[t,state_matrix,~,ier] = correction(V_init,const.mu.em,0,2);
%store results in struct
halo_L1.("orbit" + string(n_orb)).time = t;
halo_L1.("orbit" + string(n_orb)).state = state_matrix;

%find n_star gradient vector for pseudo-arclength direction
V_new = [state_matrix(1,1:6)';t(end)];
 %define constraint vector F
F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
%define Jacobian
jacob = zeros(length(F),length(V_new)); %preallocate
% extract state matrix
phi = reshape(state_matrix(end,7:42),6,6);
jacob(1:6,1:6) = phi - eye(6);
jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
% redefine constraint vector and Jacobian for modified-constraint &
% continuation, based on input flag
if constraint == 1              %constrain x0 = 0
    F(4) = [];  F(6) = V_new(1);
    jacob(4,:) = [];    jacob(6,:) = [1,0,0,0,0,0,0];
elseif constraint == 2          %constrain y0 = 0
    F(5) = [];  F(6) = V_new(2);
    jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
elseif constraint == 3          %constrain z0 = 0
    F(6) = [];  F(6) = V_new(3);
    jacob(6,:) = [];    jacob(6,:) = [0,0,1,0,0,0,0];
end
n_star = null(jacob);       %n_star is in the nullspace of DF
if n_star(7) > 0
    n_star = -n_star;   %go in direction of decreasing period
end

%find distance to the Moon for stopping condition
index = find(abs(halo_L1.("orbit" + string(n_orb)).state(:,3)) == ...
    min(abs(halo_L1.("orbit" + string(n_orb)).state(:,3))));
r2 = sqrt((halo_L1.("orbit" + string(n_orb)).state(index,1) - 1 + const.mu.em)^2 + ...
    halo_L1.("orbit" + string(n_orb)).state(index,2)^2 + halo_L1.("orbit" + string(n_orb)).state(index,3)^2);

%begin loop for pseudo-arclength continuation stop at step_max, error or
%too close to the Moon
while (n_orb < step_max +1 && ier ~= 1 && r2 > (const.radius.Moon + 100)/const.l_star.em)
    V_new = [state_matrix(1,1:6)';t(end)];     %new guess is previous initial state
    %reintegrate with added step
    [t,state_matrix,~,ier] = corrArcLength(V_new,const.mu.em,n_star,0,2,step);

    %find n_star gradient vector for pseudo-arclength direction
     %define constraint vector F
    F = state_matrix(end,1:6)' - state_matrix(1,1:6)';
    %define Jacobian
    jacob = zeros(length(F),length(V_new)); %preallocate
    % extract state matrix
    phi = reshape(state_matrix(end,7:42),6,6);
    jacob(1:6,1:6) = phi - eye(6);
    jacob(:,7) = threeBP_EOM(t(end),state_matrix(end,1:6)',const.mu.em);
    % redefine constraint vector and Jacobian for modified-constraint &
    % continuation, based on input flag
    if constraint == 1              %constrain x0 = 0
        F(4) = [];  F(6) = V_new(1);
        jacob(4,:) = [];    jacob(6,:) = [1,0,0,0,0,0,0];
    elseif constraint == 2          %constrain y0 = 0
        F(5) = [];  F(6) = V_new(2);
        jacob(5,:) = [];    jacob(6,:) = [0,1,0,0,0,0,0];
    elseif constraint == 3          %constrain z0 = 0
        F(6) = [];  F(6) = V_new(3);
        jacob(6,:) = [];    jacob(6,:) = [0,0,1,0,0,0,0];
    end
    n_star_old = n_star;        %store previous direction vector
    n_star = null(jacob);       %n_star is in the nullspace of DF
    if dot(n_star,n_star_old) < 0   %go in same direction as previous n_star
        n_star = -n_star;
    end
   
    %store results in struct
    halo_L1.("orbit" + string(n_orb)).time = t;
    halo_L1.("orbit" + string(n_orb)).state = state_matrix;
    %find distance to the Moon for stopping condition
    index = find(abs(halo_L1.("orbit" + string(n_orb)).state(:,3)) == ...
        min(abs(halo_L1.("orbit" + string(n_orb)).state(:,3))));
    r2 = sqrt((halo_L1.("orbit" + string(n_orb)).state(index,1) - 1 + const.mu.em)^2 + ...
        halo_L1.("orbit" + string(n_orb)).state(index,2)^2 + ...
        halo_L1.("orbit" + string(n_orb)).state(index,3)^2);
    n_orb = n_orb + 1;  %update step counter
end

%plot output trajectories
if const.vrb
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    figure(); hold on; grid on;
    surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35);
    plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:10:numel(fieldnames(halo_L1))
        X = halo_L1.("orbit" + string(i)).state(:,1);
        Y = halo_L1.("orbit" + string(i)).state(:,2);
        Z = halo_L1.("orbit" + string(i)).state(:,3);
        plot3(X,Y,Z,'b');
    end %plotting for loop
    axis equal;
end %if const.vrb

%% Problem 5 - Halo Orbit Post-Processing
clear; clc; const = getConst(1);

%load in previous data
%load("L1_halo_northern.mat");

%find the index of the trajectory with minimum Z-norm, this should be the
%dividing orbit between the northern and southern family
% norms = zeros(numel(fieldnames(halo_L1)),1);
% for i = 1:numel(fieldnames(halo_L1))
%     Z = halo_L1.("orbit" + string(i)).state(:,3);
%     norms(i) = norm(Z);
% end
% index = find(norms == min(norms));
% count = 1;  %start counter for south_halo_L1 indexing
% %move the southern orbits to their own struct
% for i = index:-1:1
%     south_halo_L1.("orbit" + string(count)).state = halo_L1.("orbit" + string(i)).state;
%     south_halo_L1.("orbit" + string(count)).time = halo_L1.("orbit" + string(i)).time;
%     count = count +1;
% end
% %move the northern orbits to their own struct
% count = 1;  %reset counter for north_halo_L1 indexing
% for i = index+1:numel(fieldnames(halo_L1))
%     north_halo_L1.("orbit" + string(count)).state = halo_L1.("orbit" + string(i)).state;
%     north_halo_L1.("orbit" + string(count)).time = halo_L1.("orbit" + string(i)).time;
%     count = count +1;
% end
% %load in other data file, which replaces halo_L1
% load("L1_halo_southern.mat");
% %add these orbits to previously created southern struct
% for i = 1:numel(fieldnames(halo_L1))
%     south_halo_L1.("orbit" + string(index+i)).state = halo_L1.("orbit" + string(i)).state;
%     south_halo_L1.("orbit" + string(index+i)).time = halo_L1.("orbit" + string(i)).time;
%     count = count +1;
% end
load('L1_halo_orbits.mat');
figure(); hold on; grid on;

%define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    h1 = surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35,'DisplayName',"Moon");
    h2 = plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10,'DisplayName',"$L_1$");
%plot northern orbits/
for i = 1:20:numel(fieldnames(north_halo_L1))
    X = north_halo_L1.("orbit" + string(i)).state(:,1);
    Y = north_halo_L1.("orbit" + string(i)).state(:,2);
    Z = north_halo_L1.("orbit" + string(i)).state(:,3);
    if i == 1
        h3 = plot3(X,Y,Z,'r','DisplayName','Planar Orbit');
    else
        h4 = plot3(X,Y,Z,'b','DisplayName',"Northern Family");
    end
end
xlabel("X-Position (Nondimensional Length Units");
ylabel("Y-Position");
zlabel("Z-Position");
title("Earth-Moon $L_1$ Halo Orbits (Northern Family)",'interpreter','latex');
legend([h1 h2 h3 h4],'interpreter','latex')
axis equal;

%figure(); hold on; grid on;
%define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    h1 = surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35,'DisplayName',"Moon");
    h2 = plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10,'DisplayName',"$L_1$");
%plot southern orbits
for i = 1:20:numel(fieldnames(south_halo_L1))
    X = south_halo_L1.("orbit" + string(i)).state(:,1);
    Y = south_halo_L1.("orbit" + string(i)).state(:,2);
    Z = south_halo_L1.("orbit" + string(i)).state(:,3);
    if i == 1
        h3 = plot3(X,Y,Z,'r','DisplayName','Planar Orbit');
    else
        h5 = plot3(X,Y,Z,'k','DisplayName',"Southern Family");
    end
end
xlabel("X-Position (Nondimensional Length Units");
ylabel("Y-Position");
zlabel("Z-Position");
title("Earth-Moon $L_1$ Halo Orbits (Both Families)",'interpreter','latex');
legend([h1 h2 h3 h4 h5],'interpreter','latex')
axis equal;
% %save reorganized data
% save("L1_halo_orbits.mat",'north_halo_L1','south_halo_L1');
% delete("L1_halo_southern.mat","L1_halo_northern.mat");

%prepare requested output
periods_north = zeros(numel(fieldnames(north_halo_L1)),1);  %preallocate
states_north = zeros(numel(fieldnames(north_halo_L1)),6);
%extract period and initial state vector from each trajectory
for i = 1:numel(fieldnames(north_halo_L1))            
    periods_north(i) = north_halo_L1.("orbit" + string(i)).time(end);
    states_north(i,:) = north_halo_L1.("orbit" + string(i)).state(1,1:6);
end
%reverse order to start at NRHO, go to planar, then back to south NRHO

periods_south = zeros(numel(fieldnames(south_halo_L1)),1);  %preallocate
states_south = zeros(numel(fieldnames(south_halo_L1)),6);
for i = 1:numel(fieldnames(south_halo_L1))
    periods_south(i) = south_halo_L1.("orbit" + string(i)).time(end);
    states_south(i,:) = south_halo_L1.("orbit" + string(i)).state(1,1:6);
end
periods_days_north = (periods_north.*const.t_star.em)./(3600*24);
periods_days_south = (periods_south.*const.t_star.em)./(3600*24);
jacobi_consts_north = jacobi(states_north,const.mu.em);
jacobi_consts_south = jacobi(states_south,const.mu.em);
periods = [flip(periods_north);periods_south];  %merge vectors
states = [flip(states_north);states_south];
jacobi_consts = jacobi(states,const.mu.em);     %calculate jacobi constants
periods_days = (periods.*const.t_star.em)./(3600*24);   %dimensionalize period

%plot period vs. jacobi constant
if const.vrb
    figure();
    plot(periods_days,jacobi_consts,'LineWidth',2);
    xlabel("Orbit Period in Days");
    ylabel("Jacobi Constant");
    title("Orbit Period vs. Jacobi Constant Along the $L_1$ Halo Family",...
        'interpreter','latex');
end %if const.vrb

%plot period vs. jacobi constant for southern only
if const.vrb
    figure(); hold on;
    plot(periods_days_south(1:20:end),jacobi_consts_south(1:20:end),'Marker','o','Color','blue','DisplayName','Southern Family');
    plot(periods_days_north(10:20:end),jacobi_consts_north(10:20:end),'Marker','o','Color','red','DisplayName','Northern Family');
    xlabel("Orbit Period in Days");
    ylabel("Jacobi Constant");
    title("Orbit Period vs. Jacobi Constant Along the $L_1$ Halo Family",...
        'interpreter','latex');
    legend('Location','northwest');
end %if const.vrb