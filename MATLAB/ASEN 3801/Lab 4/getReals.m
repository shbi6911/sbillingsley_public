% Function to take in eigenvalues and output real, negative gains

% K1 = angular rate gain (p, q)
% K2 = angular gain (phi, theta)
% K3 = inertial velocity gain (v, u)

% "gains" must be Nx3 with columns [K1, K2, K3]

function [realEigs, realGains]= getReals(eigVals, gains)

% Logical indexing (true if negative real)
    eigLogicVec1= eigVals(:,1)==real(eigVals(:,1)) & eigVals(:,1)<=0;
    eigLogicVec2= eigVals(:,2)==real(eigVals(:,2)) & eigVals(:,2)<=0;
    eigLogicVec3= eigVals(:,3)==real(eigVals(:,3)) & eigVals(:,3)<=0;
% Logical indexing ( true if row i has columns 1:3 all == 1 (real) )
    eigLogicVec= eigLogicVec1 & eigLogicVec2 & eigLogicVec3==1 ;

% All real, negative eigenvalues from input eigVals
    realEigs= eigVals(eigLogicVec,:);
% Gains that give real, negative eigenvalues
    realGains= gains(eigLogicVec,:);

end