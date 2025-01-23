%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function [PHI_rr,PHI_rv,PHI_vr,PHI_vv] = CWmatrix(n,t)
%CWmatrix generates matrices at a specified time point using the solution
%to the Clohessy-Wiltshire equations for relative motion in the two-body
%relative motion problem
%
%INPUTS:    n       mean motion of the 2BP orbit
%           t       time past t=0
%
%OUTPUTS    PHI_rr, PHI_rv  position matrices for CW solution
%           PHI_vr, PHI_vv  velocity matrices for CW solution
%
    PHI_rr = [4-(3*cos(n*t)), 0, 0;
              6*(sin(n*t)-(n*t)), 1, 0;
              0, 0, cos(n*t)];
    PHI_rv = [(1/n)*(sin(n*t)), (2/n)*(1-cos(n*t)), 0;
              (2/n)*(cos(n*t)-1), (1/n)*(4*sin(n*t)-(3*n*t)), 0;
              0, 0, (1/n)*sin(n*t)];
    PHI_vr = [3*n*sin(n*t), 0, 0;
              6*n*(cos(n*t)-1), 0, 0;
              0, 0, -n*sin(n*t)];
    PHI_vv = [cos(n*t), 2*sin(n*t), 0;
              -2*sin(n*t), (4*cos(n*t))-3, 0;
              0, 0, cos(n*t)];
end