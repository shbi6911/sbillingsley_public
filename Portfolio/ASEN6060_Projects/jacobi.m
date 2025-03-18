function jc = jacobi(state,mu)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     1-16-2025
%
% jacobi calculates the jacobi constant of an input state matrix
%INPUTS:    state       a 6xn element nondimensional state matrix
%           mu          the mass ratio of the system within which the state
%                       matrix is defined
%OUTPUTS:   jc          a vector of values of the jacobi constant at all
%                       time points given in the input state
r1 = sqrt((state(:,1) + mu).^2 + state(:,2).^2 + state(:,3).^2);
r2 = sqrt((state(:,1)- 1 + mu).^2 + state(:,2).^2 + state(:,3).^2);
jc = (state(:,1).^2 + state(:,2).^2) + (2*(1-mu))./r1 + (2*mu)./r2 -...
    state(:,4).^2 - state(:,5).^2 - state(:,6).^2;
end