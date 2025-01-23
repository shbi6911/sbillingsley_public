function linearShooting(L)

%Solves the BVP y'' = p(x)y' + q(x)y + r(x), for a<x<b, with the boundary
%conditions y(a)=alpha and y(b)=beta.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS.  Change these to adjust for the problem you are solving.

a = 0;  b = L;             %the endpoints of the interval, a<x<b.
h = 12;                    %space between points on x axis.
alpha = 0;  beta = 0;       %boundary values.  y(a)=alpha, y(b)=beta.
p = @(x) 0;     %continuous function
q = @(x) 3e-7;      %positive continuous function
r = @(x) (25/3)*10^-9*x^2-((50*L/6e9)*x);       %continuous function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Main part of the code.  Solves numerically the two IVP systems with
%ode45, and then combines the results to form the solution y to the BVP.

t = a:h:b;

[ ~, y1 ] = ode45( @odefun1, t, [alpha,0] );
[ ~, y2 ] = ode45( @odefun2, t,     [0,1] );

y1 = y1(:,1);  y2 = y2(:,1);

y = y1 + (beta-y1(end)) / y2(end) * y2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots the numerical solution y
maxy = max(y);
maxt = t(find(y==max(y)));
meany = mean(y);


figure(1), clf, hold('on')
plot( t, y, 'k', 'lineWidth', 2, color='blue' )
plot( t, y, 'k.', 'markerSize', 20 )
yline(maxy,color='red', lineWidth=2)
xline(maxt,color='red',LineWidth=2)
plot (maxt, maxy,'p', 'MarkerFaceColor','red','markersize', 15)
set( gca, 'fontSize', 15 )
xlabel('x (in)'), ylabel('y(x) (in)')
title('Deflection of Beam with Two End Supports')
grid('on')
drawnow, hold('off')
disp(meany)
disp(maxy)
disp(maxt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The two ODE functions that are passed into ode45

    function u = odefun1(t,y)
        u = zeros(2,1);
        u(1) = y(2);
        u(2) = p(t)*y(2) + q(t)*y(1) + r(t);
    end

    function u = odefun2(t,y)
        u = zeros(2,1);
        u(1) = y(2);
        u(2) = p(t)*y(2) + q(t)*y(1);
    end

end