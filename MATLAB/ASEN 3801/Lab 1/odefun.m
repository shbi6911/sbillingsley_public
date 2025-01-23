function [y_prime] = odefun(t,y_vec,var)
%INPUTS     t   scalar time value given by ode45
%           y   4-element state vector
%
%OUTPUTS    y_prime     4-element derivative state vector
%
%METHODOLOGY    odefun takes in a state vector and outputs the derivate 
% of that state vector according to given equations

w = y_vec(1);   x=y_vec(2);     y=y_vec(3);     z=y_vec(4);

y_prime = [0;0;0;0];
y_prime(1)= (-9*w) + y;
y_prime(2)= (4*w*x*y) - x^2;
y_prime(3)= (2*w) - x -(2*z);
y_prime(4)= (x*y) - y^2 - (3*z^3);

end