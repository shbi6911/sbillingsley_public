%By:        Shane Billingsley
%Class:     APPM 2350 Calculus 3 for Engineers
%Date:      Fall 2022

% f(x) = 980.6160*(1-(0.0026373*cos(2*x))+5.9e-6*(cos(2*x))^2);
% h(x,y) = y*(3.085462e-4+(2.27e-7*cos(2*x)));
% j(x,y) = (y^2)*(7.254e-11+(1.0e-18*(cos(2*x))));
% k(x,y) = (y^3)*(1.517e-17+(6e-20*cos(2*x)));
% g(x,y) = f(x)-h(x,y)+j(x,y)-k(x,y);
g = @(x,y) (980.6160*(1-(0.0026373*cos(2*x))+5.9e-6*(cos(2*x))^2))-...
    (y*(3.085462e-4+(2.27e-7*cos(2*x))))+...
    ((y^2)*(7.254e-11+(1.0e-18*(cos(2*x)))))-...
    ((y^3)*(1.517e-17+(6e-20*cos(2*x))));
format long;
disp("g at Baseline Road");
disp(g(40,1609.34));
disp("g 10km above North Pole");
disp(g(90,10000));
disp("halfway between North Pole and equator at mean sea level");
disp(g(45,0));
disp("10000 feet above -66.5 degrees latitude");
disp(g(-66.5,3048));

