function drdt = threeBP_EOM(t,r,mu)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025

% threeBP_EOM calculates the equations of motion for the third body within
% a Circular Restricted Three Body Problem, using nondimensional quantities
% for time and state vector, and a rotating coordinate frame.  It is
% intended for use with ODE solvers such as ODE45 or ODE89.
%
% INPUTS:   t       nondimensional time unit (from ODE solver)
%           r       nondimensional state vector, as a column vector, of the
%                   form [x;y;z;x_dot;y_dot;z_dot]
%           mu      mass ratio parameter
%
% OUTPUTS:  drdt    nondim derivative state vector, as a column vector, 
%                   of the form [x_dot;y_dot;z_dot;x_ddot;y_ddot;z_ddot]

    drdt = zeros(6,1);      %preallocate output vector
    drdt(1:3) = r(4:6);     % derivative of position is velocity
    %calculate respective radii
    r1 = sqrt((r(1) + mu)^2 + r(2)^2 + r(3)^2);
    r2 = sqrt((r(1) - 1 + mu)^2 + r(2)^2 + r(3)^2);
    %calculate accelerations per CR3BP equations of motion
    %note x,y,z = r(1,2,3)  x_dot,y_dot,z_dot = r(4,5,6)
    drdt(4) = 2*r(5) + r(1) - ((1-mu)*(r(1)+mu))/r1^3 - (mu*(r(1)-1+mu))/r2^3;
    drdt(5) = -2*r(4) + r(2) - ((1-mu)*r(2))/r1^3 - (mu*r(2))/r2^3;
    drdt(6) = -((1-mu)*r(3))/r1^3 - (mu*r(3))/r2^3;
end