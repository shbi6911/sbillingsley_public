function [shear,mom] = moment(w0,w1,L,mag,loc)
C1 = sum(mag);
C2 = sum(mag.*loc);
V = @(x) -((w1-w0)/(2*L))*x.^2 - (w0*x) + C1;
M = @(x) ((w1-w0)/(6*L))*x.^3 + ((w0/2)*x.^2) - (C1*x) + C2;
shear = V(loc);
%disp(shear);
mom = M(loc);
% disp(mom);
end