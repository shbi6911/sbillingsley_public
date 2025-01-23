function [x, t, refMat] = getDataState(filename, t_trim)
    arguments
        filename string
        t_trim (1,1) double = 0
    end
% Loading necessary variables from minidrone data file
    load(filename, 'rt_estim', 'rt_cmd');
    x= rt_estim.signals.values; % matrix of quadrotor states
    t= rt_estim.time; % column vector of times
    cmds= rt_cmd.signals.values; % matrix of commanded positions/angles

% Switching indices 4:6 (Euler angles) to reflect order of roll-pitch-yaw
    x(:,4:6)= x(:,6:-1:4); % [rad]
% Removing data prior to desired input time "t_trim"
    t_logic= t >= t_trim;
    t= t(t_logic); % [s]
    x= x(t_logic, :);
% Shifting time vector by first value so t(1) is 0
    t= t - t(1); % [s]

% Reference positions
    zE_ref= cmds(t_logic,2:4); % [m]
% Reference angles
    euler_ref= cmds(t_logic,6:8); % [rad]
% Concatenating into matrix of ref. positions & angles
    refMat= [zE_ref, euler_ref]; 

% fc= rt_motor.signals.values; % matrix of motor forces (column=motor #)

% % need to figure out what units of fc are (can convert to rad/s but not sure
% % how this is useful)
% % Calculating control force and moments vector ctrls=[Zc; Lc; Mc; Nc]
% A_coeff= [   -1,            -1,          -1,          -1; ...
%         -(d/sqrt(2)),  -(d/sqrt(2)),  (d/sqrt(2)), (d/sqrt(2)); ...
%          (d/sqrt(2)),  -(d/sqrt(2)), -(d/sqrt(2)), (d/sqrt(2)); ...
%              km,            -km,          km,         -km];

end


