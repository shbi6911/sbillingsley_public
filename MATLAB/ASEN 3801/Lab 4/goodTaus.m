% Function to take in real, negative eigenvalues and determine which
% give time constants that are all <= tauMax

function [tau, gains] = goodTaus(realEigs, realGains, tauMax)

% Ensuring input "goodEigs" contains only real eigenvalues
if isreal(realEigs) == 0
    tau=0;
    gains=0;
    disp('The input eigenvalues have imaginary components')
    return

else
% Determining which eigenvalues (if any) give time constants < 1.25s
    tauMat= -1./realEigs;
% Logical indexing (true if tau less than 1.25s & non-negative)
    tauLogicVec1= tauMat(:,1)<=tauMax & tauMat(:,1)>=0;
    tauLogicVec2= tauMat(:,2)<=tauMax & tauMat(:,2)>=0;
    tauLogicVec3= tauMat(:,2)<=tauMax & tauMat(:,3)>=0;
% Logical indexing [ true if row i has columns 1:3 == 1 ( 0=<tau<=tauMax ) ]
    tauLogicVec= tauLogicVec1 & tauLogicVec2 & tauLogicVec3==1 ;

% All real, negative eigenvalues from input eigVals
    tau= tauMat(tauLogicVec,:);
% Gains that give real, negative eigenvalues
    gains= realGains(tauLogicVec,:);
end

end