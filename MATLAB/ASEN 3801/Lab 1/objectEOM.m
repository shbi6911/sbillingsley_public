function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind_vel)
%
%NOTE:  All vector inputs/outputs should use NED axes
%
%INPUTS     t   scalar value of current integration time
%           x   six-element state vector (position x;y;z;velocity x;y;z)
%           rho scalar value of local air density
%           Cd  scalar value of object drag coefficient
%           A   scalar value of object cross-sectional area
%           m   scalar value of object mass
%           g   scalar value of local gravitational acceleration
%           wind_vel    three element wind velocity vector (x;y;z)
%
%OUTPUTS    xdot    six-element derivative state vector
%
%METHODOLOGY    objectEOM creates a derivative state vector for an object
%   moving through space under Newton's 2nd Law, assuming gravity and drag
%   as the only accelerations.  Wind, air density, gravity and the physical
%   properties of the object are all assumed constant.

a_grav = [0;0;g];   %define a gravitational acceleration vector

v_e = x(4:6) - wind_vel;    %find air relative velocity vector
v_a = norm(v_e);            %find airspeed
D = 0.5*Cd*A*rho*(v_a^2);   %calculate magnitude of drag force
a_drag = (-D.*(v_e./v_a))./m;   %define a drag acceleration vector

%create derivative state vector.  Velocity is derivative of position, and
%calculated acceleration is derivative of velocity
xdot = [x(4:6);(a_grav+a_drag)];

end