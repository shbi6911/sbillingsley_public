%By:        Shane Billingsley
%Class:     APPM 2360 Intro to Diff Eqns with Linear Algebra
%Date:      Summer 2022

%define sinusoidal functions
f = @(t,y) (0.75*y) - (0.005*y^2) - ((1.5*y^3)/(y^3 + 1.25)*(abs(sin((pi/365)*t))));
g = @(t,y) (0.75*y) - (0.05*y^2) - ((1.5*y^3)/(y^3 + 1.25)*(abs(sin((pi/365)*t))));
h = @(t,y) (0.75*y) - (0.1*y^2) - ((1.5*y^3)/(y^3 + 1.25)*(abs(sin((pi/365)*t))));

figure(1);  %figure for b=0.005
hold on;
title('Seasonal Variation in Klingon Hunting of Brown Tribbles (b=0.005');
xlabel ('Time in Days');
ylabel ('Hundreds of Tribbles')
ax = gca;
ax.YLim = [-1 200];
ax.XLim = [-1 2000];

tspan = [0,2000];  %tspan for all solutions (2000 days)
[t1,y1] = ode45(f,tspan,1.5);  %solution for initial pt 1.5
plot(t1,y1,'blue',LineWidth=1);

[t2,y2] = ode45(f,tspan,0.5);   %solution for initial pt 0.5
plot(t2,y2,'magenta',LineWidth=1);

[t3,y3] = ode45(f,tspan,0.05); %solution for initial pt 0.05
plot(t3,y3,'green',LineWidth=1);

[ta,ya] = ode45(f,tspan,200);  %solution for initial pt 200
plot(ta,ya,'black',LineWidth=1);
legend('Solution for initial point (0,1.5)', ...
    'Solution for initial point (0,0.5)', ...
    'Solution for initial point (0,0.05)', ...
    'Solution for initial point (0,200');

figure(2);  %figure for b=0.05
hold on;
title('Seasonal Variation in Klingon Hunting of White Tribbles (b=0.05');
xlabel ('Time in Days');
ylabel ('Hundreds of Tribbles')
ax = gca;
ax.YLim = [-1 22];
ax.XLim = [-1 2000];

[t4,y4] = ode45(g,tspan,14);    %solution for initial pt 14
plot(t4,y4,'blue',LineWidth=1);

[t5,y5] = ode45(g,tspan,2.5);   %solution for initial pt 2.5
plot(t5,y5,'green',LineWidth=1);

[t6,y6] = ode45(g,tspan,1.8);   %solution for initial pt 1.8
plot(t6,y6,'magenta',LineWidth=1);

[t7,y7] = ode45(g,tspan,0.05);  %solution for initial pt 0.05
plot(t7,y7,'red',LineWidth=1);

[tb,yb] = ode45(g,tspan,20);    %solution for initial pt 20
plot(tb,yb,'black',LineWidth=1);
legend('Solution for initial point (0,14)', ...
    'Solution for initial point (0,2.5)', ...
    'Solution for initial point (0,1.8)', ...
    'Solution for initial point (0,0.05)', ...
    'Solution for initial point (0,20');

figure(3);   %figure for b=0.1
hold on;
title('Seasonal Variation in Klingon Hunting of Grey Tribbles (b=0.1');
xlabel ('Time in Days');
ylabel ('Hundreds of Tribbles')
ax = gca;
ax.YLim = [-1 12];
ax.XLim = [-1 2000];

[t8,y8] = ode45(h,tspan,0.05);    %solution for initial pt 0.05
plot(t8,y8,'blue',LineWidth=1);

[t9,y9] = ode45(h,tspan,3);       %solution for initial pt 3
plot(t9,y9,'magenta',LineWidth=1);

[t10,y10] = ode45(h,tspan,6.5);   %solution for initial pt 6.5
plot(t10,y10,'green',LineWidth=1);

[tc,yc] = ode45(h,tspan,20);     %solution for initial pt 20
plot(tc,yc,'black',LineWidth=1);
legend('Solution for initial point (0,0.05)', ...
    'Solution for initial point (0,3)', ...
    'Solution for initial point (0,6.5)', ...
    'Solution for initial point (0,20');