function [H_an] = slope_an(k,A,V,I)
%slope_an   This function calculates an analytical slope value for heat
%transfer
%INPUTS     k   :thermal conductivity
%           A   :cross-sectional area
%           V   :voltage
%           I   :current
%OUTPUTS    H_an    :scalar value of analytical slope
P=I*V;
H_an=P/(k*A);
end