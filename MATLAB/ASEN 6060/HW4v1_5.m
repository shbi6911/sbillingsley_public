% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-18-2025

%% Problem 1a
%calculate Earth_Moon L2 Lyapunov family using pseudo-arclength continuation.
%Investigate stability using the eigenvalues and eigenvectors of the
%monodromy matrix
tic;
clear; clc; const = getConst(1);

%generate an initial guess for an L2 Lyapunov orbit
%construct linearized A matrix at L2 location
A = stab(const.mu.em,const.eq_pts.em(2,:));
%get eigenvalues/vectors from linearized A matrix
[V1,D1] = eig(A);
% create basis for eigenspace of lambda3
basis1 = round(real(V1(:,3)),15);
basis2 = round(imag(V1(:,3)),15);
% generate initial state for integration @ L2
V_init = [const.eq_pts.em(2,:)';(2*pi)/imag(D1(3,3))];
% add initial perturbation
scale_factor = 1/10000;   %set scaling factor for initial guess
V_init(1:6) = V_init(1:6) + basis1.*scale_factor;

%set continuation step size
step = 0.001;    
%set maximum number of steps
step_max = 5000;
%initialize step counter
n_orb = 1;
%set modified-constraint flag
constraint = 2; %constrain y0
%correct initial guess using modified constraint in y0
[t,state_matrix,~,ier] = correction(V_init,const.mu.em,0,2);
%store results in struct
lyap_L2.("orbit" + string(n_orb)).time = t;
lyap_L2.("orbit" + string(n_orb)).state = state_matrix;
%store eigenvalues and eigenvectors of monodromy matrix
[lyap_L2.("orbit" + string(n_orb)).vec,lyap_L2.("orbit"+string(n_orb)).val] = ...
    eig(reshape(state_matrix(end,7:42),6,6));

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
if n_star(7) < 0
    n_star = -n_star;   %go in direction of increasing period
end
 %find distance to the Moon
index = find(diff(sign(state_matrix(2:end-1,2))))+2;
r2 = sqrt((state_matrix(index,1) - 1 + const.mu.em)^2 + ...
state_matrix(index,2)^2 + state_matrix(index,3)^2);

%begin loop for pseudo-arclength continuation stop at step_max or error
while (n_orb < step_max && ier ~= 1 && r2 > const.radius.Moon/const.l_star.em)
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
    lyap_L2.("orbit" + string(n_orb)).time = t;
    lyap_L2.("orbit" + string(n_orb)).state = state_matrix;
    %store eigenvalues and eigenvectors of monodromy matrix
    [lyap_L2.("orbit" + string(n_orb)).vec,lyap_L2.("orbit"+string(n_orb)).val] = ...
        eig(reshape(state_matrix(end,7:42),6,6));
    n_orb = n_orb + 1;  %update step counter
    V_new = [state_matrix(1,1:6)';t(end)];  %update free variable vector
    %find distance to the Moon
    index = find(diff(sign(state_matrix(2:end-1,2))))+2;
    r2 = sqrt((state_matrix(index,1) - 1 + const.mu.em)^2 + ...
    state_matrix(index,2)^2 + state_matrix(index,3)^2);
end %end while (n_orb < step_max +1 && ier ~= 1 && r2 > const.radius.Moon/const.l_star.em)

%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10);
    circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k',''); hold on;
    for i = 1:50:numel(fieldnames(lyap_L2))
        X = lyap_L2.("orbit" + string(i)).state(:,1);
        Y = lyap_L2.("orbit" + string(i)).state(:,2);
        plot(X,Y,'b');
    end %plotting for loop
    axis equal;
end %if const.vrb
save("L2_lyapunov.mat",'lyap_L2');
fprintf("Elapsed Time is %d minutes.\n",toc/60);

%% Problem 1a Post-Processing
%this section loads trajectory data file created in the previous section
%sort eigenvalues and plot stability indices for L2 Lyapunov orbits
clear; clc; const = getConst(1);

load("L2_lyapunov.mat");
% extract eigenvalues and transform to stability indices
eig_values = zeros(numel(fieldnames(lyap_L2)),6);  %preallocate
%extract eigenvalues
for i = 1:numel(fieldnames(lyap_L2))
    eig_values(i,:) = diag(lyap_L2.("orbit" + string(i)).val);
end
%remove trivial values and sort eigenvalues into reciprocal pairs
[eig_trim,eig_sorted,~] = sortEigs(eig_values);
%generate vector of orbital periods for the family
periods = zeros(numel(fieldnames(lyap_L2)),1);  %preallocate
for i = 1:numel(fieldnames(lyap_L2))
    periods(i) = lyap_L2.("orbit" + string(i)).time(end);
end
periods_days = (periods.*const.t_star.em)./(3600*24);   %dimensionalize
%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    h1 = plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10,'DisplayName',"$L_2$");
    h2 = circle((1-const.mu.em),0,const.radius.Moon/const.l_star.em,'k','Moon'); hold on;
    for i = 1:50:numel(fieldnames(lyap_L2))
        X = lyap_L2.("orbit" + string(i)).state(:,1);
        Y = lyap_L2.("orbit" + string(i)).state(:,2);
        h3 = plot(X,Y,'b','DisplayName','Lyapunov Family');
    end %plotting for loop
    axis equal;
    xlabel("X-Position (Nondimensional Length Units");
    ylabel("Y-Position");
    title("Earth-Moon $L_2$ Lyapunov Orbit Family",'interpreter','latex');
    legend([h1 h2 h3],'Location','northwest','interpreter','latex');
end %if const.vrb

%plot stability indices
if const.vrb
    figure();   hold on;    grid on;
    plot(periods_days,eig_sorted(:,1) + eig_sorted(:,2),'LineWidth',2);
    yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_1 + \lambda_2$)",'interpreter','latex');
    title("Stability of $L_2$ Lyapunov Orbit Family, Inplane",'interpreter','latex');
    figure();   hold on;    grid on;
    plot(periods_days,eig_sorted(:,3) + eig_sorted(:,4),'LineWidth',2);
     yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_3 + \lambda_4$)",'interpreter','latex');
    title("Stability of $L_2$ Lyapunov Orbit Family, Outplane",'interpreter','latex');
end

%% Problem 1b
%%calculate Earth_Moon L2 halo family using pseudo-arclength continuation.
%Investigate stability using the eigenvalues and eigenvectors of the
%monodromy matrix

tic;
clear; clc; const = getConst(1);
%input given intial guess
V_init = [1.180462,0,-0.0209998,0,-0.158363,0,3.411921]';

%set continuation step size
step = 0.001;    
%set maximum number of steps
step_max = 5000;
%initialize step counter
n_orb = 1;
%set modified-constraint flag
constraint = 2; %constrain y0
%correct initial guess using modified constraint in y0
[t,state_matrix,~,ier] = correction(V_init,const.mu.em,0,2);
%store results in struct
halo_L2.("orbit" + string(n_orb)).time = t;
halo_L2.("orbit" + string(n_orb)).state = state_matrix;
%store eigenvalues and eigenvectors of monodromy matrix
[halo_L2.("orbit" + string(n_orb)).vec,halo_L2.("orbit"+string(n_orb)).val] = ...
    eig(reshape(state_matrix(end,7:42),6,6));

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
if n_star(7) < 0
    n_star = -n_star;   %go in direction of increasing period
end

%find distance to the Moon for stopping condition
index = find(abs(state_matrix(:,3)) == min(abs(state_matrix(:,3))));
r2 = sqrt((state_matrix(index,1) - 1 + const.mu.em)^2 + ...
    state_matrix(index,2)^2 + state_matrix(index,3)^2);

%begin loop for pseudo-arclength continuation stop at step_max or error
while (n_orb < step_max && ier ~= 1 && r2 > const.radius.Moon/const.l_star.em)
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
    halo_L2.("orbit" + string(n_orb)).time = t;
    halo_L2.("orbit" + string(n_orb)).state = state_matrix;
    %store eigenvalues and eigenvectors of monodromy matrix
    [halo_L2.("orbit" + string(n_orb)).vec,halo_L2.("orbit"+string(n_orb)).val] = ...
        eig(reshape(state_matrix(end,7:42),6,6));
    n_orb = n_orb + 1;  %update step counter
    %find distance to the Moon for stopping condition
    index = find(abs(state_matrix(:,3)) == min(abs(state_matrix(:,3))));
    r2 = sqrt((state_matrix(index,1) - 1 + const.mu.em)^2 + ...
    state_matrix(index,2)^2 + state_matrix(index,3)^2);

    V_new = [state_matrix(1,1:6)';t(end)];  %update free variable vector
end %end while (n_orb < step_max +1 && ier ~= 1)

%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35);
    plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:20:numel(fieldnames(halo_L2))
        X = halo_L2.("orbit" + string(i)).state(:,1);
        Y = halo_L2.("orbit" + string(i)).state(:,2);
        Z = halo_L2.("orbit" + string(i)).state(:,3);
        if i == 1
            plot3(X,Y,Z,'r');
        else
            plot3(X,Y,Z,'b');
        end
    end %plotting for loop
    axis equal;
end %if const.vrb
save("L2_halo.mat",'halo_L2');
fprintf("Elapsed Time is %d minutes.\n",toc/60);

%% Problem 1b Post-Processing
%this section loads trajectory data file created in the previous section
%sort eigenvalues and plot stability indices for L2 halo orbits of the
%northern family
clear; clc; const = getConst(1);

load("L2_halo.mat");
% extract eigenvalues and transform to stability indices
eig_values = zeros(numel(fieldnames(halo_L2)),6);  %preallocate
%extract eigenvalues
for i = 1:numel(fieldnames(halo_L2))
    eig_values(i,:) = diag(halo_L2.("orbit" + string(i)).val);
end
%remove trivial values and sort eigenvalues into reciprocal pairs
[eig_trim,eig_sorted,~] = sortEigs(eig_values);
%generate vector of orbital periods for the family
periods = zeros(numel(fieldnames(halo_L2)),1);  %preallocate
for i = 1:numel(fieldnames(halo_L2))
    periods(i) = halo_L2.("orbit" + string(i)).time(end);
end
periods_days = (periods.*const.t_star.em)./(3600*24);   %dimensionalize
%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    h1 = surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35,'DisplayName',"Moon");
    h2 = plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10,'DisplayName',"$L_2$");
    for i = 1:30:numel(fieldnames(halo_L2))
        X = halo_L2.("orbit" + string(i)).state(:,1);
        Y = halo_L2.("orbit" + string(i)).state(:,2);
        Z = halo_L2.("orbit" + string(i)).state(:,3);
        if i == 1
            h3 = plot3(X,Y,Z,'r','DisplayName',"Initial Orbit");
        else
            h4 = plot3(X,Y,Z,'b','DisplayName',"Northern Halo Family");
        end
    end %plotting for loop
    axis equal;
    xlabel("X-Position (Nondimensional Length Units");
    ylabel("Y-Position");
    zlabel("Z-Position");
    title("Earth-Moon $L_2$ Halo Orbit Family",'interpreter','latex');
    legend([h1 h2 h3 h4],'Location','northwest','interpreter','latex');
end %if const.vrb

%plot stability indices
if const.vrb
    figure();   grid on;
    plot(periods_days,eig_sorted(:,1) + eig_sorted(:,2),'LineWidth',2);
    yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_1 + \lambda_2$)",'interpreter','latex');
    title("Stability of $L_2$ Halo Orbit Family",'interpreter','latex');
    figure();   hold on;    grid on;
    plot(periods_days,eig_sorted(:,3) + eig_sorted(:,4),'LineWidth',2);
     yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_3 + \lambda_4$)",'interpreter','latex');
    title("Stability of $L_2$ Halo Orbit Family",'interpreter','latex');
end
%% Problem 1c
%%calculate Earth_Moon L2 axial family using pseudo-arclength continuation.
%Investigate stability using the eigenvalues and eigenvectors of the
%monodromy matrix

tic;
clear; clc; const = getConst(1);
%input given initial guess
%V_init = [1.0301513;0;0;0;0.7030025;0.1552945;4.312367];
%input initial guess taken from the far side of initial orbit
V_init = [1.218838;0;0;0;-0.425673;-0.0409125;4.312367]';

%set continuation step size
step = 0.0005;
%set maximum number of steps
step_max = 4000;
%initialize step counter & planar checker & planar flag
n_orb = 1;      points = zeros(step_max,5);     planar_flag = 0;
%set modified-constraint flag
constraint = 2;
%correct initial guess using modified constraint in y0
[t,state_matrix,~,ier] = correction(V_init,const.mu.em,0,2);

%store results in struct
axial_L2.("orbit" + string(n_orb)).time = t;
axial_L2.("orbit" + string(n_orb)).state = state_matrix;
%store eigenvalues and eigenvectors of monodromy matrix
[axial_L2.("orbit" + string(n_orb)).vec,axial_L2.("orbit"+string(n_orb)).val] = ...
    eig(reshape(state_matrix(end,7:42),6,6));
%store difference between Z-values at max and min y-value
points(n_orb,1) = find(state_matrix(:,2) == min(state_matrix(:,2)));
points(n_orb,2) = find(state_matrix(:,2) == max(state_matrix(:,2)));
points(n_orb,3) = state_matrix(points(n_orb,1),3);
points(n_orb,4) = state_matrix(points(n_orb,2),3);
points(n_orb,5) = abs(points(n_orb,4) - points(n_orb,3));


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

%begin loop for pseudo-arclength continuation stop at step_max or error
while (n_orb <= step_max && ier ~= 1 && planar_flag ~= 1)
    %reintegrate with added step
    [t,state_matrix,norms,ier] = corrArcLength(V_new,const.mu.em,n_star,0,2,step,100,1e-10);
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
    axial_L2.("orbit" + string(n_orb)).time = t;
    axial_L2.("orbit" + string(n_orb)).state = state_matrix;
    %store eigenvalues and eigenvectors of monodromy matrix
    [axial_L2.("orbit" + string(n_orb)).vec,axial_L2.("orbit"+string(n_orb)).val] = ...
        eig(reshape(state_matrix(end,7:42),6,6));
    n_orb = n_orb + 1;  %update step counter
    V_new = [state_matrix(1,1:6)';t(end)];  %update free variable vector
    %store difference between Z-values at max and min y-value
    points(n_orb,1) = find(state_matrix(:,2) == min(state_matrix(:,2)));
    points(n_orb,2) = find(state_matrix(:,2) == max(state_matrix(:,2)));
    points(n_orb,3) = state_matrix(points(n_orb,1),3);
    points(n_orb,4) = state_matrix(points(n_orb,2),3);
    points(n_orb,5) = abs(points(n_orb,4) - points(n_orb,3));
    %check for stopping condition
    %after the first 10 orbits, check to see if the difference between the
    %Z-coord @ min y-coord for this orbit and the first orbit is small
    if n_orb > 10
        planar_diff = abs(points(n_orb,3) - points(1,3));
        if planar_diff < 5e-4
            planar_flag = 1;
        end % if planar_diff < 5e-4
    end % if n_orb > 10
end %end while (n_orb < step_max +1 && ier ~= 1 && planar_flag ~= 1)

%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35);
    plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:30:numel(fieldnames(axial_L2))
        X = axial_L2.("orbit" + string(i)).state(:,1);
        Y = axial_L2.("orbit" + string(i)).state(:,2);
        Z = axial_L2.("orbit" + string(i)).state(:,3);
        if i == 1
            plot3(X,Y,Z,'r');
        else
            plot3(X,Y,Z,'b');
        end
    end %plotting for loop
    axis equal;
end %if const.vrb
save("L2_axial.mat",'axial_L2');
fprintf("Elapsed Time is %d minutes.\n",toc/60);

%% Problem 1c Post-Processing
%this section loads trajectory data file created in the previous section
%sort eigenvalues and plot stability indices for L2 axial orbits
clear; clc; const = getConst(1);

load("L2_axial.mat");
% extract eigenvalues and transform to stability indices
eig_values = zeros(numel(fieldnames(axial_L2)),6);  %preallocate
%extract eigenvalues
for i = 1:numel(fieldnames(axial_L2))
    eig_values(i,:) = diag(axial_L2.("orbit" + string(i)).val);
end
%remove trivial values and sort eigenvalues into reciprocal pairs
[eig_trim,eig_sorted,~] = sortEigs(eig_values);
%generate vector of orbital periods for the family
periods = zeros(numel(fieldnames(axial_L2)),1);  %preallocate
for i = 1:numel(fieldnames(axial_L2))
    periods(i) = axial_L2.("orbit" + string(i)).time(end);
end
periods_days = (periods.*const.t_star.em)./(3600*24);   %dimensionalize
%plot output trajectories
if const.vrb
    figure(); hold on; grid on;
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    h1 = surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','FaceAlpha',0.35,'DisplayName',"Moon");
    h2 = plot3(const.eq_pts.em(2,1),0,0,'.','Color','k','MarkerSize',10,'DisplayName',"$L_2$");
    for i = 1:50:numel(fieldnames(axial_L2))
        X = axial_L2.("orbit" + string(i)).state(:,1);
        Y = axial_L2.("orbit" + string(i)).state(:,2);
        Z = axial_L2.("orbit" + string(i)).state(:,3);
        if i == 1
            h3 = plot3(X,Y,Z,'r','DisplayName',"Initial Orbit");
        else
            h4 = plot3(X,Y,Z,'b','DisplayName',"Axial Family");
        end
    end %plotting for loop
    axis equal;
    xlabel("X-Position (Nondimensional Length Units");
    ylabel("Y-Position");
    zlabel("Z-Position");
    title("Earth-Moon $L_2$ Axial Orbit Family",'interpreter','latex');
    legend([h1 h2 h3 h4],'Location','northwest','interpreter','latex');
end %if const.vrb

%plot stability indices
if const.vrb
    figure();   grid on;
    plot(periods_days,eig_sorted(:,1) + eig_sorted(:,2),'LineWidth',2);
    yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_1 + \lambda_2$)",'interpreter','latex');
    title("Stability of $L_2$ Axial Orbit Family",'interpreter','latex');
    figure();   hold on;    grid on;
    plot(periods_days,eig_sorted(:,3) + eig_sorted(:,4),'LineWidth',2);
     yline(2,'r');   yline(-2,'r');
    xlabel("Orbit Period (days)");
    ylabel("Stability Index ($\lambda_3 + \lambda_4$)",'interpreter','latex');
    title("Stability of $L_2$ Axial Orbit Family",'interpreter','latex');
end

%% Problem 2a
%compute stable and unstable manifolds for an L1 Lyapunov orbit.  Plot
%until the x-coordinate of the trajectories reaches the Moon and the
%Earth's x-coordinate (for positive and negative trajectories respectively)
%Then plot for 6 nondimensional time units

clear; clc; const = getConst(1);
load("L1_lyapunov.mat");    %load previous data
periods = zeros(numel(fieldnames(lyap_L1)),1);  %preallocate
%locate an L1 Lyapunov orbit with a period ~2.72
for i = 1:numel(fieldnames(lyap_L1))
    periods(i) = lyap_L1.("orbit" + string(i)).time(end);
end
index = find(abs(periods-2.72) == min(abs(periods-2.72)));
%extract matrix of states for this periodic orbit
state_matrix = lyap_L1.("orbit" + string(index)).state;
N = 50;                     %set number of points
d = 50/const.l_star.em;     %set step distance
%find manifolds
[pos_stab_mani,neg_stab_mani,pos_un_mani,neg_un_mani] = compManifold(state_matrix,N,d);
%set tolerances and initial conditions
EventFun = @(t,y)orbitEvents(t,y,const.mu.em);
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14,'Events',EventFun);
tspan = [0 6];
EOMfun = @(t,r)threeBP_EOM(t,r,const.mu.em);
%integrate trajectories until event function triggers
for i = 1:length(pos_stab_mani)
    name = "point" + string(i);
    initial_state = pos_stab_mani(i,:)';
    % integrate stable positive direction
    [t,state_matrix] = ode113(EOMfun, flip(tspan), initial_state, opts);
    L1_mani.(name).short.stab.pos.time = t;           %store data in struct
    L1_mani.(name).short.stab.pos.state = state_matrix;
    % integrate stable negative direction
    initial_state = neg_stab_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, flip(tspan), initial_state, opts);
    L1_mani.(name).short.stab.neg.time = t;           %store data in struct
    L1_mani.(name).short.stab.neg.state = state_matrix;
    % integrate unstable pos direction
    initial_state = pos_un_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    L1_mani.(name).short.un.pos.time = t;           %store data in struct
    L1_mani.(name).short.un.pos.state = state_matrix;
    % integrate unstable neg direction
    initial_state = neg_un_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    L1_mani.(name).short.un.neg.time = t;           %store data in struct
    L1_mani.(name).short.un.neg.state = state_matrix;
end %for i = 1:length(pos_stab_mani)

%reset opts without event function
opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
%integrate trajectories for a full 6 time units
for i = 1:length(pos_stab_mani)
    name = "point" + string(i);
    initial_state = pos_stab_mani(i,:)';
    % integrate stable positive direction
    [t,state_matrix] = ode113(EOMfun, flip(tspan), initial_state, opts);
    L1_mani.(name).long.stab.pos.time = t;           %store data in struct
    L1_mani.(name).long.stab.pos.state = state_matrix;
    % integrate stable negative direction
    initial_state = neg_stab_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, flip(tspan), initial_state, opts);
    L1_mani.(name).long.stab.neg.time = t;           %store data in struct
    L1_mani.(name).long.stab.neg.state = state_matrix;
    % integrate unstable pos direction
    initial_state = pos_un_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    L1_mani.(name).long.un.pos.time = t;           %store data in struct
    L1_mani.(name).long.un.pos.state = state_matrix;
    % integrate unstable neg direction
    initial_state = neg_un_mani(i,:)';
    [t,state_matrix] = ode113(EOMfun, tspan, initial_state, opts);
    L1_mani.(name).long.un.neg.time = t;           %store data in struct
    L1_mani.(name).long.un.neg.state = state_matrix;
end %for i = 1:length(pos_stab_mani)

%plotting positive manifolds (Moon direction) until event function
if const.vrb
    figure();   hold on; grid on;
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.35);
    plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:numel(fieldnames(L1_mani))
        name = "point" + string(i);
        %plot positive stable manifold
        X = L1_mani.(name).short.stab.pos.state(:,1);
        Y = L1_mani.(name).short.stab.pos.state(:,2);
        Z = L1_mani.(name).short.stab.pos.state(:,3);
        plot3(X,Y,Z,'b');
        %plot positive unstable manifold
        X = L1_mani.(name).short.un.pos.state(:,1);
        Y = L1_mani.(name).short.un.pos.state(:,2);
        Z = L1_mani.(name).short.un.pos.state(:,3);
        plot3(X,Y,Z,'r');
    end % for i = 1:numel(fieldnames(L1_mani))
    axis equal;
end %if const.vrb

%plotting negative manifolds (Earth direction) until event
if const.vrb
    figure();   hold on; grid on;
    %define a sphere for the Earth
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Earth/const.l_star.em) + (- const.mu.em); 
    Y = Y.*(const.radius.Earth/const.l_star.em); 
    Z = Z.*(const.radius.Earth/const.l_star.em);
    surf(X,Y,Z,'FaceColor','b','EdgeColor','b','FaceAlpha',0.35);
    plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:numel(fieldnames(L1_mani))
        name = "point" + string(i);
        %plot negative stable manifold
        X = L1_mani.(name).short.stab.neg.state(:,1);
        Y = L1_mani.(name).short.stab.neg.state(:,2);
        Z = L1_mani.(name).short.stab.neg.state(:,3);
        plot3(X,Y,Z,'b');
        %plot negative unstable manifold
        X = L1_mani.(name).short.un.neg.state(:,1);
        Y = L1_mani.(name).short.un.neg.state(:,2);
        Z = L1_mani.(name).short.un.neg.state(:,3);
        plot3(X,Y,Z,'r');
    end % for i = 1:numel(fieldnames(L1_mani))
    axis equal;
end %if const.vrb

%plotting positive manifolds (Moon direction) for full 6 time units
if const.vrb
    figure();   hold on; grid on;
     %define a sphere for the Earth
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Earth/const.l_star.em) + (- const.mu.em); 
    Y = Y.*(const.radius.Earth/const.l_star.em); 
    Z = Z.*(const.radius.Earth/const.l_star.em);
    surf(X,Y,Z,'FaceColor','b','EdgeColor','b','FaceAlpha',0.35);
    %define a sphere for the Moon
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Moon/const.l_star.em) + (1 - const.mu.em); 
    Y = Y.*(const.radius.Moon/const.l_star.em); 
    Z = Z.*(const.radius.Moon/const.l_star.em);
    surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.35);
    plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:numel(fieldnames(L1_mani))
        name = "point" + string(i);
        %plot positive stable manifold
        X = L1_mani.(name).long.stab.pos.state(:,1);
        Y = L1_mani.(name).long.stab.pos.state(:,2);
        Z = L1_mani.(name).long.stab.pos.state(:,3);
        plot3(X,Y,Z,'b');
        %plot positive unstable manifold
        X = L1_mani.(name).long.un.pos.state(:,1);
        Y = L1_mani.(name).long.un.pos.state(:,2);
        Z = L1_mani.(name).long.un.pos.state(:,3);
        plot3(X,Y,Z,'r');
    end % for i = 1:numel(fieldnames(L1_mani))
    axis equal;
end %if const.vrb

%plotting negative manifolds (Earth direction) for full 6 time units
if const.vrb
    figure();   hold on; grid on;
    %define a sphere for the Earth
    [X,Y,Z] = sphere(100);
    X = X.*(const.radius.Earth/const.l_star.em) + (- const.mu.em); 
    Y = Y.*(const.radius.Earth/const.l_star.em); 
    Z = Z.*(const.radius.Earth/const.l_star.em);
    surf(X,Y,Z,'FaceColor','b','EdgeColor','b','FaceAlpha',0.35);
    plot3(const.eq_pts.em(1,1),0,0,'.','Color','k','MarkerSize',10);
    for i = 1:numel(fieldnames(L1_mani))
        name = "point" + string(i);
        %plot negative stable manifold
        X = L1_mani.(name).long.stab.neg.state(:,1);
        Y = L1_mani.(name).long.stab.neg.state(:,2);
        Z = L1_mani.(name).long.stab.neg.state(:,3);
        plot3(X,Y,Z,'b');
        %plot negative unstable manifold
        X = L1_mani.(name).long.un.neg.state(:,1);
        Y = L1_mani.(name).long.un.neg.state(:,2);
        Z = L1_mani.(name).long.un.neg.state(:,3);
        plot3(X,Y,Z,'r');
    end % for i = 1:numel(fieldnames(L1_mani))
    axis equal;
end %if const.vrb

%% Event Function

function [value,isterminal,direction] = orbitEvents(t,y,mu)
%this event function is designed for use with ode solvers in the CR3BP.  It
%stops integration when a trajectory reaches the yz plane defined at the
%location of P1 or P2.  It requires an additional input in the form of the
%mass ratio (mu) for the system.

    val1 = y(1) - (1-mu);   %at P2's x-coord
    val2 = y(1) - (-mu);    %at P1's x-coord
    value = [val1;val2];    %when either of the above equal zero
    isterminal = [1;1];     %stop integration
    direction = [];         %from any direction
end
