%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

%rotate ground station over one day
initialpos = [-4460.49;2682.08;-3674.26]; %initial vector of ground station [km]
tspan= 0:86400; %seconds in a day
rotangle = tspan*(360/86400); %rotation angle in degrees per second
rotanglerad = rotangle.*(pi/180); %rotation angle in radians per second
axis = [0;0;1]; %rotate around the z-axis

grnd_rotations = zeros(3,length(rotanglerad));
for i=1:length(tspan)
    newvec = rotate(rotanglerad(i),axis,initialpos);
    grnd_rotations(:,i)=newvec;
end

%calculate satellite trajectory
% initial position and velocity of satellite [km; km/s]
sat1 = [1986.21;6388.28;-1237.15;-4.931;0.403;-5.829];
% calculate satellite orbit using ode45
orbitfun = @OrbitEOM;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, Sat1Full] = ode45(orbitfun,tspan,sat1,options);

figure (1);             %plot trajectory
[X,Y,Z] = sphere(100);
hold on; grid on;
%ground station
plot3(-4460.49,2682.08,-3674.26,'.','MarkerSize',30,'Color',[0.9290 0.6940 0.1250]);
%satellite
plot3(Sat1Full(:,1),Sat1Full(:,2),Sat1Full(:,3),'r','LineWidth',2);
%Earth surface
surf((X*6378),(Y*6378),(Z*6378),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.35);
%Earth center
plot3(0,0,0,'r.','MarkerSize',30);
title ("Orbital Path of a Satellite Over One Day");
xlabel ("X-Position [km]");
ylabel ("Y-Position [km]");
zlabel ("Z-Position [km]");
legend ("Canberra Ground Station Initial Location","Satellite Orbital Path");

%determine visibility
visibility = zeros(1,length(tspan));
for i=1:length(tspan)
    satpos = (Sat1Full(i,1:3))';        %satellite position
    grndpos= grnd_rotations(:,i);       %ground station position
    relpos=satpos-grndpos;              %vector from ground to sat
    [N_local,E_local,D_local]=NED(grndpos); %local NED axes
    %dot product between relative pos vector and local down vector
    visibility(i)=dot(relpos,D_local);
end

figure (2); hold on; grid on;
plot(tspan,visibility);         %plot dot product
yline(0,'r');                   %line at zero dot product
xlabel ("Time [s]");
ylabel ("Dot Product");
title("Satellite Visibility Over One Day");

function [N_local,E_local,D_local] = NED(invec)
%NED    This function constructs local NED axes for a position
%Inputs:    invec is the COLUMN VECTOR position in Earth Centered
%Outputs:   N_local: local North COLUMN VECTOR
%           E_local: local East COLUMN VECTOR
%           D_local: local Down COLUMN VECTOR

    omega = atan2(invec(2),invec(1));  %longitude [rad]
    a=6378137;                       %Earth dimensions [m]
    b=6356752.3142;
                                     %latitude [rad]
    alpha = atan((a^2/b^2)*(invec(3)/sqrt(invec(1)^2+invec(2)^2)));
    N0=[0;0;1];                     %NED UNIT vectors @ (0,0,0)
    E0=[0;1;0];
    D0=[-1;0;0];
    
    E_local=rotate(omega,N0,E0);        %local East vector
    N_local=rotate(alpha,-E_local,N0);   %local North vector
    D_local=cross(N_local,E_local);      %local Down vector
end

function [outvec] = rotate(theta,axis,invec)
%ROTATE This function rotates a vector by an angle theta
%   Inputs: Angle (this is the SCALAR rotation angle)
%           Axis (this is a unit COLUMN VECTOR in the direction of the axis)
%           Vector (this is the COLUMN VECTOR to be rotated)
%   Outputs:    Outvec (this is the new vector after rotation)    
    
    axis = axis./norm(axis);
    nmatrix = axis*axis';           %n*n^T
    crossmatrix = [0 -axis(3) axis(2); axis(3) 0 -axis(1); %cross product
                        -axis(2) axis(1) 0];
    a = cos(theta);
    b = 1-a;
    c= sin(theta);
    %define axis-angle rotation matrix
    rotatematrix = (b*nmatrix)+(a*eye(3))+(c*crossmatrix);
    outvec = rotatematrix*invec; %mulitply for resultant vector
   
end

function dsdt = OrbitEOM(t,s)
% ODE function Handle for Satellite States
% Function is used by ODE45 to generate the Satellite States.
% Inputs: t    -> time  s -> States (size 6x1) [x;y;z;vx;vy;vz];
% Output: dsdt -> Rate of change of States (size 6x1)

    mu = 3.986e5;               % GM of Earth
    
    r(:,1) = s(1:3);            % Position Vector
    v(:,1) = s(4:6);            % Velocity Vector
    
    normR = norm(r);            % Magnitude of position
    dsdt = [v;-mu*r/normR^3];   % Rate of change of States  
    
end