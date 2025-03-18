function [Uxx,Uxy,Uxz,Uyy,Uyz,Uzz] = partialDeriv(mu,states,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     3-04-2025
%
% partialDeriv calculates the second partial derivatives of the 
% pseudo-potential function within the CR3BP, referenced to a given input state.
% partialDeriv can take in nx6 matrix input. It assumes that each row of an
% input matrix is one six-element state vector of the form [x,y,z,vx,vy,vz],
% and that the n rows are n different input states.  It outputs 6 column
% vectors of n elements, which are the partial derivatives at each of the n
% input states.
% INPUTS:       mu          mass ratio value for the system of interest
%               states      nx6 matrix of n different 6-element state vectors
%               vrb         standard verbosity flag
%
% OUTPUTS:      Uxx         PD of x-deriv wrt x
%               Uxy         =Uyx, PD of x-deriv wrt y
%               Uxz         =Uzx, PD of x-deriv wrt z
%               Uyy         PD of y-deriv wrt y
%               Uyz         =Uzy, PD of y-deriv wrt z
%               Uzz         PD of z-deriv wrt z
    arguments
        mu double
        states(:,6) {mustBeNumeric}
        vrb logical =0
    end
    %calculate nx1 vector of distances from input state pts to P1 and P2
    r1 = vecnorm([states(:,1)+mu,states(:,2),states(:,3)],2,2);
    r2 = vecnorm([(states(:,1) - 1 + mu),states(:,2),states(:,3)],2,2);
    %calculate nx1 vector of partial derivatives for each input state
    Uxx = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*((states(:,1)+mu).^2))./r1.^5 +...
        (3.*mu.*(states(:,1) -1 + mu).^2)./r2.^5;
    Uxy = (3.*(1-mu).*(states(:,1)+mu).*states(:,2))./r1.^5 + ...
        (3.*mu.*(states(:,1) -1 + mu).*states(:,2))./r2.^5;
    Uxz = (3.*(1-mu).*(states(:,1)+mu).*states(:,3))./r1.^5 + ...
        (3.*mu.*(states(:,1) -1 + mu).*states(:,3))./r2.^5;
    Uyy = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*(states(:,2).^2))./r1.^5 +...
        (3.*mu.*states(:,2).^2)./r2.^5;
    Uyz = (3.*(1-mu).*(states(:,2).*states(:,3)))./r1.^5 +...
        (3.*mu.*states(:,2).*states(:,3))./r2.^5;
    Uzz = -(1-mu)./r1.^3 - mu./r2.^3 + (3.*(1 - mu).*states(:,3).^2)./r1.^5 +...
        (3.*mu.*states(:,3).^2)./r2.^5;
end