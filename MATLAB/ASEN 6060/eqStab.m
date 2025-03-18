function [inplane,outplane] = eqStab(mu,vrb)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     2-11-2025
%
% eqStab is a function to output the linearized stability matrices for the
% five equilibrium points of a CR3BP system
%DEPENDENCY:        this function calls the function eqPts
%
%INPUTS:    mu      required input value of mass ratio
%           vrb     optional input for verbosity, controls output
%OUTPUTS:   inplane     a THREE dimensional matrix of the size 4,4,5, where
%                   each of the five pages is a 4x4 linearized stability
%                   matrix for the in-plane stability modes of an
%                   equilibrium point, L1....L5
%           outplane    a THREE dimensional matrix of the size 2,2,5 where
%                   each page is a 2x2 stability matrix for the out of
%                   plane modes, for L1....L5
    arguments
        mu double
        vrb logical =0
    end
    eq_pts = eqPts(mu,vrb);     %locate equilibrium points
    inplane = zeros(4,4,5);     outplane = zeros(2,2,5);    %preallocate
    %calculate 5x1 vector of distances from equilibrium pts to P1 and P2
    r1 = vecnorm([eq_pts(:,1)+mu,eq_pts(:,2),eq_pts(:,3)],2,2);
    r2 = vecnorm([(eq_pts(:,1) - 1 + mu),eq_pts(:,2),eq_pts(:,3)],2,2);
    %calculate 5x1 vector of partial derivatives for each equilibrium pt
    Uzz = -(1-mu)./r1.^3 - mu./r2.^3 + (3.*(1 - mu).*eq_pts(:,3).^2)./r1.^5 +...
        (3.*mu.*eq_pts(:,3).^2)./r2.^5;
    Uxx = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*((eq_pts(:,1)+mu).^2))./r1.^5 +...
        (3.*mu.*(eq_pts(:,1) -1 + mu).^2)./r2.^5;
    Uyy = 1 - (1-mu)./r1.^3 - mu./r2.^3 + (3.*(1-mu).*(eq_pts(:,2).^2))./r1.^5 +...
        (3.*mu.*eq_pts(:,2).^2)./r2.^5;
    Uxy = (3.*(1-mu).*(eq_pts(:,1)+mu).*eq_pts(:,2))./r1.^5 + ...
        (3.*mu.*(eq_pts(:,1) -1 + mu).*eq_pts(:,2))./r2.^5;
    %build out out-of-plane matrices
    for ii = 1:size(outplane,3)
        outplane(1,2,ii) = 1;
        outplane(2,1,ii) = Uzz(ii);
    end
    %build in-plane matrices
    for jj = 1:size(inplane,3)
        inplane(1,3,jj) = 1;    inplane(2,4,jj) = 1;
        inplane(3,4,jj) = 2;    inplane(4,3,jj) = -2;
        inplane(3,1,jj) = Uxx(jj);  inplane(4,2,jj) = Uyy(jj);
        inplane(3,2,jj) = Uxy(jj);  inplane(4,1,jj) = Uxy(jj);
    end
end