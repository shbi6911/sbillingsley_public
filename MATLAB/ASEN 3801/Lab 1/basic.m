tspan = [0 20]; x_ic = [0;0;0;0;20;-20];
opts = odeset('Events',@HitGround,'RelTol',tol,'AbsTol',tol);
[t_out_basic,x_out_basic,~,~,~] = ode45(@basic_EOM,tspan,x_ic,opts);

figure(); view(3);
hold on; grid on;
plot3(x_out_i(:,1),x_out_i(:,2),x_out_i(:,3),'LineWidth',2);
plot3(x_out_basic(:,1),x_out_basic(:,2),x_out_basic(:,3),'LineWidth',2,'Color','r');
set(gca,'ZDir','reverse');  %reverse axes so that -z is on top
title("Initial Trajectory with Zero Wind Speed");
xlabel("X-Position (m)");
ylabel("Y-Position (m)");
zlabel("Z-Position (m)");
hold off;

function x_out = basic_EOM(t,x)

g= 9.81;    a_grav = [0;0;g];
p_prime = x(4:6);
v_prime = a_grav;
x_out = [p_prime;v_prime];
end

function [zpos,isterminal,direction] = HitGround(t,x)
%INPUTS     t       scalar time value
%           x       6-element state vector
%OUTPUTS    zpos    z-coordinate of position vector
%           isterminal  Boolean flag of whether to stop integration or not
%           direction   allowable direction to approach the stop point from
%METHODOLOGY    HitGround is an Event function as required by ode45.  It is
%   checking for when zpos = 0.  Isterminal is set to 1 so integration will
%   stop when this condition is detected.  Direction is set to 0 so that
%   the function is agnostic of direction.  Essentially, HitGround will
%   detect when zpos=0 (i.e. when the projectile has hit the ground) and
%   will stop integration at that point.
zpos = x(3);
isterminal = 1;
direction = 0;
end