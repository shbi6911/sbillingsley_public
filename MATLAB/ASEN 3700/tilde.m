%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function [a_tilde] = tilde(a)
% tilde converts a vector into its corresponding cross product matrix to
% perform cross products directly
%
%INPUTS     a           a 3-element vector in Euclidean space
%
%OUTPUTS    a_tilde     the corresponding cross product matrix

a_tilde = [0 -a(3) a(2);a(3) 0 -a(1);-a(2) a(1) 0];
end