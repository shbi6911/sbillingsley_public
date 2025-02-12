function eq_pts = eqPts(mu,vrb,guess)
% By:       Shane Billingsley
% Class:    ASEN 6060 Advanced Astrodynamics
% Date:     2-11-2025
%
% eqPts is a function to determine the equilibrium points of a system under
% the CR3BP, with the system defined by a value for mass ratio
%
%INPUTS:    mu      required input value of mass ratio
%           vrb     optional input for verbosity, controls output
%           guess   optional initial guesses for L1,L2,L3, without input
%                   these are determined programmatically
%OUTPUTS:   eq_pts  5x6 matrix of state vectors at each of the equilibrium
%                   points, in order from top to bottom (L1;L2;L3;L4;L5)
    arguments
        mu double
        vrb logical =0
        guess (1,3) =[0,0,0]
    end
    %define function for collinear equilibrium point location
    fun = @(x,mu) x - ((1 - mu)*(x + mu))/((abs(x+mu))^3) - ...
        (mu*(x - 1 + mu))/((abs(x - 1 + mu))^3);
    %preallocate output matrix
    eq_pts = zeros(5,6);
    %define initial guesses
    if norm(guess) == 0
        guess = [(1-(mu/3)^(1/3));(1+(mu/3)^(1/3));...
                 (-(1 +(5/12)*mu))];
    end
    %set options for rootfinding function
    if vrb
        opts = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1e-16);
    else
        opts = optimoptions(@fsolve,'Display','none','FunctionTolerance',1e-16);
    end
    
    %find all three collinear equilibrium points using fsolve
    for i = 1:length(guess)
        eq_pts(i,1) = fsolve(@(x)fun(x,mu),guess(i),opts);
    end
    % calculate the locations of L4 and L5 using analytical solutions
    eq_pts(4,1:2) = [0.5 - mu,sqrt(3)/2];
    eq_pts(5,1:2) = [0.5 - mu,-sqrt(3)/2];
end