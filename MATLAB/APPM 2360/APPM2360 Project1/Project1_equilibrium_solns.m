%By:        Shane Billingsley
%Class:     APPM 2360 Intro to Diff Eqns with Linear Algebra
%Date:      Summer 2022

% fun = @(y)((0.75*y) - (0.005*y^2) - ((1.5*y^3)/(y^3 + 1.25)));
% y0 = [1 2];
% x = fzero(fun,y0);
% print ('zero for 0.005 is ',x);

% gun = @(y)((0.75*y) - (0.05*y^2) - ((1.5*y^3)/(y^3 + 1.25)));
% y0 = [0.5 1.5];
% y1 = [1.5 2.5];
% y2 = [12 13];
% x = fzero(gun,y0)
% z = fzero(gun,y1)
% q = fzero(gun,y2)


hun = @(y)((0.75*y) - (0.1*y^2) - ((1.5*y^3)/(y^3 + 1.25)));
y0 = [0.5 1.5];
x = fzero(hun,y0)
