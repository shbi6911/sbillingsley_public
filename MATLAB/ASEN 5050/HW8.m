%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Problem 1
clear; clc; const = getConst();
%givens
r_p = 7500;     r_a = 8500;     R_planet = 6500;   %km
i = deg2rad(105);   P_sat = 110*60;
a_planet = 2.25*const.AU;           %planet's orbit radius in km
a_sat = (r_p + r_a)/2;              %satellite's semi-major axis
e_sat = (r_a - r_p)/(r_a + r_p);

%determine planetary mass
mu_planet = (a_sat^3*4*pi^2)/P_sat^2;
M_planet = mu_planet/const.G;

%find planetary orbit period
P_planet = 2*pi*sqrt(a_planet^3/const.mu.Sun);
Omega_dot = (2*pi)/P_planet;

%determine J2 coefficient
J2 = ((-2/3)*(Omega_dot/cos(i))*((1-e_sat^2)^2)*(a_sat^(7/2)))/(sqrt(mu_planet)*R_planet^2);

%% Problem 2
clear; clc; const = getConst();
%givens
R0_vec = [2489.63813;-3916.07418;-5679.05524];
V0_vec = [9.13452;-1.91212;2.57306];
R0 = norm(R0_vec);      V0 = norm(V0_vec);
X0 = [R0_vec;V0_vec];

%characterize initial orbit
[a0,e0,i0,RAAN_0,omega_0,theta_star_0] = CartesianToKepler(R0_vec,V0_vec,const.mu.Earth);
energy_0 = -const.mu.Earth/(2*a0);
h0 = sqrt(const.mu.Earth*a0*(1 - e0^2));
R1 = a0*(1+e0);
p0 = a0*(1-e0^2);
delta_theta = pi - theta_star_0;

%determine reference vector at apoapsis using f and g functions
f = 1-(R1./p0).*(1-cos(delta_theta));
g = (R1.*R0.*sin(delta_theta))./sqrt(const.mu.Earth*p0);
coeff = ((1-cos(delta_theta))./p0) - (1./R1) - (1./R0); 
f_dot = sqrt(const.mu.Earth/p0).*tan(delta_theta./2).*coeff;
g_dot = 1 - (R0./p0).*(1-cos(delta_theta));

R1_vec = R0_vec*f + V0_vec*g;                    %calculate arrays of vectors
V1_vec = R0_vec*f_dot + V0_vec*g_dot;
X_ref = [R1_vec;V1_vec];

%determine elapsed time from t0 to t1
thing = sqrt((1-e0)/(1+e0));
E0 = 2*atan(thing*tan(theta_star_0/2));
P0 = (2*pi*sqrt(a0^3/const.mu.Earth));
t1_0 = P0/2 - ((P0/(2*pi))*(E0 - e0*sin(E0)));

%numerically integrate trajectory
tspan = [0 t1_0];
opts = odeset('RelTol',10e-8,'AbsTol',10e-8);
EOMfun = @(t,state)OrbitEOM(t,state,const.mu.Earth);
[t,state_matrix_0] = ode45(EOMfun, tspan, X0,opts);

%calculate time history of energy and momentum
energy = zeros(1,length(state_matrix_0));
momentum = zeros(1,length(state_matrix_0));
for i=1:length(state_matrix_0)
    momentum(i) = norm(cross(state_matrix_0(i,1:3),state_matrix_0(i,4:6)));
    energy(i) = (dot(state_matrix_0(i,4:6),state_matrix_0(i,4:6)))/2 - ...
        const.mu.Earth/norm(state_matrix_0(i,1:3));
end
mom_range = abs(max(momentum) - min(momentum));
eng_range = abs(max(energy)-min(energy));

%plotting

% figure();
% plot(t,energy);
% title("Specific Mechanical Energy over Time");
% ylim([mean(energy) - 2*eng_range, mean(energy) + 2*eng_range]);
% xlabel("Time (s)");
% ylabel("$\varepsilon$ ($\frac{km^2}{s^2}$)",'interpreter','latex','FontSize',18,'Rotation',0);
% 
% figure();
% plot(t,momentum);
% title("Specific Angular Momentum over Time");
% ylim([mean(momentum) - 2*mom_range, mean(momentum) + 2*mom_range]);
% xlabel("Time (s)");
% ylabel("h ($\frac{km^2}{s^2}$)",'interpreter','latex','FontSize',18,'Rotation',0);
% 
% figure(); hold on; grid on;
% [X,Y,Z] = sphere(100);
% surf((X*const.radius.Earth),(Y*const.radius.Earth),(Z*const.radius.Earth),...
%     'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.15,...
%     'EdgeAlpha',0.15);
% plot3(state_matrix_0(:,1),state_matrix_0(:,2),state_matrix_0(:,3),'b','LineWidth',2);
% scatter3(0,0,0,20,'filled','r');
% arrow([0;0;0],R0_vec);
% arrow(state_matrix_0(50,1:3),state_matrix_0(51,1:3),'Length',30,'Width',30,'Color','b');
% arrow(state_matrix_0(100,1:3),state_matrix_0(101,1:3),'Length',30,'Width',30,'Color','b');
% arrow(state_matrix_0(150,1:3),state_matrix_0(151,1:3),'Length',30,'Width',30,'Color','b');
% title("Numerically Integrated Trajectory");
% xlabel("GCRF X-Position (km)");
% ylabel("GCRF Y-Position (km)");
% zlabel("GCRF Z-Position (km)");
% legend("Earth Surface","Orbital Path","Earth Center","Initial Pos Vector","Direction of Motion");
% axis equal

% adjusting tolerances
tol = [1e-4;1e-6;1e-8;1e-10;1e-12];
times = zeros(1,length(tol));
for i= 1:length(tol)
    tic;
    opts = odeset('RelTol',tol(i),'AbsTol',tol(i));
    [t,state.("trajectory" + string(i))] = ode45(EOMfun, tspan, X0,opts);
    times(i) = toc;
end
%% Problem 2d

for j = 1:length(tol)
    trajectory = state.("trajectory" + string(j));
    state.("deltaR" + string(j)) = norm(trajectory(end,1:3)' - X_ref(1:3));
    state.("deltaV" + string(j)) = norm(trajectory(end,4:6)' - X_ref(4:6));
    momentum = norm(cross(trajectory(end,1:3),trajectory(end,4:6)));
    energy = (dot(trajectory(end,4:6),trajectory(end,4:6)))/2 - ...
        const.mu.Earth/norm(trajectory(end,1:3));
    state.("deltaEps" + string(j)) = energy - energy_0;
    state.("deltaH" + string(j)) = momentum - h0;
end

    
%% functions
 function state_dot = OrbitEOM(t,state,mu)
     R_ddot = -mu/(norm(state(1:3))^3)*state(1:3);
     state_dot = [state(4:6);R_ddot];
 end