function drdt = threeBP_refTraj(t,r,mu,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-4-2025

% threeBP_refTraj calculates the equations of motion for a CR3BP trajectory.
% It also generates a State Transition Matrix at the input time step, so
% that the trajectory being calculated can be used as a reference
% trajectory.  It expects an input column state vector of the form 
% [x;y;z;vx;vy;vz;STM], where the 6x6 STM has been transformed into a
% vector via the command reshape(STM,[],1).  This should transform the
% matrix into a column vector COLUMNWISE, that is [c1;c2;c3;c4;c5;c6]
%
% DEPENDENCY:       this function calls the partialDeriv function
%
% INPUTS:   t       nondimensional time unit (from ODE solver)
%           r       nondimensional state vector, as a column vector, of the
%                   form [x;y;z;vx;vy;vz;STM] (see above)
%           mu      mass ratio parameter of the system of interest
%
% OUTPUTS:  drdt    nondim derivative state vector, as a column vector, 
%                   of the form [x_dot;y_dot;z_dot;x_ddot;y_ddot;z_ddot;STM]
arguments
        t  double
        r(42,1) {mustBeNumeric}
        mu double
        vrb logical =0
    end
    drdt = zeros(42,1);      %preallocate output vector
    drdt(1:3) = r(4:6);     % derivative of position is velocity
    %calculate respective distances to primaries
    r1 = sqrt((r(1) + mu)^2 + r(2)^2 + r(3)^2);
    r2 = sqrt((r(1) - 1 + mu)^2 + r(2)^2 + r(3)^2);
    %calculate accelerations per CR3BP equations of motion
    %note x,y,z = r(1,2,3)  x_dot,y_dot,z_dot = r(4,5,6)
    drdt(4) = 2*r(5) + r(1) - ((1-mu)*(r(1)+mu))/r1^3 - (mu*(r(1)-1+mu))/r2^3;
    drdt(5) = -2*r(4) + r(2) - ((1-mu)*r(2))/r1^3 - (mu*r(2))/r2^3;
    drdt(6) = -((1-mu)*r(3))/r1^3 - (mu*r(3))/r2^3;
    %reshape State Transition Matrix into 6x6 form
    STM = reshape(r(7:42),6,6);
    %construct linearized A matrix at local state
    [Uxx,Uxy,Uxz,Uyy,Uyz,Uzz] = partialDeriv(mu,r(1:6)',vrb);
    C = [Uxx,Uxy,Uxz;Uxy,Uyy,Uyz;Uxz,Uyz,Uzz];
    D = zeros(3,3);     D(1,2) = 2;     D(2,1) = -2;
    A = [zeros(3,3),eye(3);C,D];
    %calculate new STM and append to output vector
    STM_dot = A*STM;    drdt(7:42) = reshape(STM_dot,[],1);
end