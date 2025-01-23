%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% find initial vectors
n_hat_2 = [0.6428;-0.7660;0];       %given vectors
h_hat_2 = [-0.3237;-0.2717;0.9063];
evec_2 = [0.0475;0.3755;0.1295];

e_2 = norm(evec_2);     %find eccentricity
e_hat_2 = evec_2./e_2;  %eccentricity unit vector

%inclination
i_2 = acos(h_hat_2(3));  i_2_deg = rad2deg(i_2);

%RAAN
OMEGA_2 = acos(n_hat_2(1));
if n_hat_2(2) < 0
    OMEGA_2 = -OMEGA_2;
end
OMEGA_2_deg = rad2deg(OMEGA_2);

%argument of periapsis
omega_2 = acos(dot(n_hat_2,e_hat_2));
if e_hat_2(3) < 0
    omega_2 = -omega_2;
end
omega_2_deg = rad2deg(omega_2);

%given radius at ascending node
r_2 = 4070.6;   %km
%at ascending node, arg of perigee is the negative of true anomaly
theta_star_2 = -omega_2;
%use conic eqn
a_2 = (r_2+(r_2*e_2*cos(theta_star_2)))/(1-e_2^2);

r_p_2 = a_2*(1-e_2);
v_p_2 = sqrt(((2*4902.799)/r_p_2)-(4902.799/a_2));
Rvec_p_2 = [r_p_2;0;0];
Vvec_p_2 = [0;v_p_2;0];
DCM_p_2 = AstroRot313(OMEGA_2,i_2,omega_2);
Rvec_p_2_inert = DCM_p_2*Rvec_p_2;
Vvec_p_2_inert = DCM_p_2*Vvec_p_2;

P_2 = 2*pi*sqrt(a_2^3/4902.799);

%% initial values
sat1 = [Rvec_p_2_inert; Vvec_p_2_inert];

%% calculate orbits
orbitfun = @OrbitEOM;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
tspan= 0:60:P_2;
[~, Sat1Full] = ode45(orbitfun,tspan,sat1,options);
%% plotting

[X,Y,Z] = sphere(100);
hold on; grid on;
plot3(Sat1Full(:,1),Sat1Full(:,2),Sat1Full(:,3),'r','LineWidth',2);
surf((X*1738),(Y*1738),(Z*1738),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.15,'EdgeColor','none');
plot3(0,0,0,'r.','MarkerSize',30);
axis('equal');

[x y] = meshgrid(-4000:500:4000); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z,'FaceAlpha',0.35) % Plot the surface

XE = [0;e_hat_2(1)*3000]; YE = [0;e_hat_2(2)*3000]; ZE = [0; e_hat_2(3)*3000];
plot3(XE,YE,ZE, 'k','LineWidth',2);

XN = [0;n_hat_2(1)*4500]; YN = [0;n_hat_2(2)*4500]; ZN= [0; n_hat_2(3)*4500];
plot3(XN,YN,ZN, 'k','LineWidth',2);

plot3([0;3000],[0;0],[0;0], 'b','LineWidth',2);
plot3([0;0],[0;3000],[0;0], 'm','LineWidth',2);
plot3([0;0],[0;0],[0;3000], 'b','LineWidth',2);