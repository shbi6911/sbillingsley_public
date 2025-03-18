function A = stab(mu,state,vrb)
%stab calculates the linearized stability matrix around an input point
%within the CR3BP
%
% INPUTS:   mu      mass ratio parameter of the system of interest
%           state   state value at which the matrix is calculated
%           vrb     optional verbosity flag
%
% OUTPUTS:  A       6x6 linearized stability matrix
%
arguments
    mu double
    state (1,6)
    vrb logical =0
end
[Uxx,Uxy,Uxz,Uyy,Uyz,Uzz] = partialDeriv(mu,state);
C = [Uxx,Uxy,Uxz;Uxy,Uyy,Uyz;Uxz,Uyz,Uzz];
D = zeros(3,3);     D(1,2) = 2;     D(2,1) = -2;
A = [zeros(3,3),eye(3);C,D];
end