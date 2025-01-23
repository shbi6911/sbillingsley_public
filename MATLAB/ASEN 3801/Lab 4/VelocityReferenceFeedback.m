%% ASEN 3801, Quadrotor Simulation and Control: Task 3
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniala
% Mohammadi
% Date Created: 4/7/2024
% Last Modified: 4/14/2024

% Control laws for roll and pitch
%       Lc= -K1*p - K2*(phi - phi_ref) - K3*(v - v_ref)
%       Mc= -K1*q - K2*(theta - theta_ref) - K3*(u - u_ref)
%   Simplifying from assumptions in 3.6 and Steady Hover trim state
%       phi_ref = theta_ref = 0
%       ref. velocity achieved instantaneously gives v_ref = u_ref = 0.5m/s
%           initial u or v should be 0.5 m/s
%       all "delta" values are equal to the current state value except zE

function [Fc, Gc] = VelocityReferenceFeedback(t, var, dynStr)
    arguments
        t double
        var double
        dynStr = 'lateral' % default dynamics if no input
    end

latStr= 'lateral';
lonStr= 'longitudinal';
t_max= 2; % [s]
lat = strcmp(dynStr,latStr);
lon = strcmp(dynStr,lonStr);
t_logic= t < t_max;

    % Getting minidrone parameters
        [m, g, ~, ~, ~, ~, ~]= getParams();
    % Extracting state vector values
        phi=var(4,:); theta=var(5,:); % [rad] roll & pitch euler angles
        u=var(7,:); v=var(8,:); % [m/s] body x & y velocities
        p=var(10,:); q=var(11,:); r=var(12,:); % [rad/s] body angular rates
    % Control force vector in body coordinates (constant)
        N=length(t);
        Fc= [0; 0; -m*g]*ones(1,N); % [N]

    % Lateral gains
        Kp= 0.005916; Kphi=0.0116 ; Kv=0.000584;
    % Longitudinal gains
        Kq= 0.00374; Ktheta= 0.0072; Ku= -0.00035936;

if t_logic==1 & lat == 1 % true for 'lateral' input
    v_ref= 0.5; % [m/s]
    u_ref= 0; % [m/s]
elseif t_logic==1 & lon == 1 % true for 'longitudinal' input
    v_ref= 0; % [m/s]
    u_ref= 0.5; % [m/s]
else % commanded velocity zero after t=2s
    v_ref= 0; % [m/s]
    u_ref= 0; % [m/s]
end
    Lc= -Kp*p - Kphi*phi - Kv*(v - v_ref); % [N*m] roll control moment
    Mc= -Kq*q - Ktheta*theta - Ku*(u - u_ref); % [N*m] pitch control moment
    Nc= -0.004*r; % [N*m] yaw control moment

    % Vertically concatenating for final 3x1 control moment vector
        Gc= [Lc; Mc; Nc]; % [N*m]
end





