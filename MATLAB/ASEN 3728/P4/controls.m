function u = controls(t, x)
% CONTROLS - computes the control inputs for the aircraft given the time and state
%
% INPUTS
%   t = time in seconds
%   x = [x, y, z, phi, theta, psi, u^E, v^E, w^E, p, q, r] in SI units and radians
%   NOTE: For this assignment, when there is wind, the velocity components are with respect to the EARTH, not the wind
%
% OUTPUTS
%   u = [delta_e, delta_a, delta_r, delta_t] radians and number beteween 0 and 1

    delta_e_0 = 0.1068;
    delta_t_0 = 0.2439;

    kp_trial = 10;          %values to adjust gains for fine control
    ka_trial = 0.05;
    kr_trial = -10;
    phi_r_trial = 0.03;
    % Insert your control gains here
    kp = 11+kp_trial;
    ka = -0.279257070194989+ka_trial;
    kr = -7.6+ka_trial;

    % Keep these gains the same
    kth = -3;
    kq = -1;

    % Keep this pitch angle
    theta_r = 0.0278;

    % Insert your roll angle here
    phi_r = 0.221096125573615+phi_r_trial;

    % Keep these elevator and thrust controllers
    delta_e = kth*(theta_r-x(5)) + kq*-x(11) + delta_e_0;
    delta_t = delta_t_0;

    % Insert your lateral controllers here
    delta_a = ka*(kp*(phi_r - x(4))-x(10)); %derived from control laws
    delta_r = -kr*x(12);

    u = [delta_e, delta_a, delta_r, delta_t];
end
