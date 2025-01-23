function [H_exp,T0_exp,H_exp_delta] = slope_exp(x,y)
%slope_exp      This function calculates slope of heat transfer from
%experimental data
%INPUTS     x   :a vector of x-values
%           y   :a vector of temperature values
%OUTPUTS    H_exp       :slope of the polyfit line between x and y
%           T0_exp      :y-intercept of the polyfit line
%           H_exp_delta :error struct output by polyfit
[fit,s1]= polyfit(x,y,1);
H_exp=fit(1);
T0_exp=fit(2);
H_exp_delta=s1;
end