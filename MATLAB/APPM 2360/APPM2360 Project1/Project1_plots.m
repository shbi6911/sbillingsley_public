%By:        Shane Billingsley
%Class:     APPM 2360 Intro to Diff Eqns with Linear Algebra
%Date:      Summer 2022

f = @(y) (0.75*y) - (0.005*y.^2) - ((1.5*y.^3)/(y.^3 + 1.25));
g = @(y) (0.75*y) - (0.05*y.^2) - ((1.5*y.^3)/(y.^3 + 1.25));
h = @(y) (0.75*y) - (0.1*y.^2) - ((1.5*y.^3)/(y.^3 + 1.25));
x = linspace(-1,150,16);
y = zeros([1 16]);

figure(1);
hold on;
fplot(f,[0,150]);
ax = gca;
ax.YLim = [-1 80];
ax.XLim = [-1 150];
plot(x,y);
z = fzero(f,[147,148]);
r = 0;
plot(z, r,'.','MarkerSize',15);
text (z-0.3,r+0.4,'147.973,0');

figure(2);
hold on;
fplot(f,[0,15]);
ax = gca;
ax.YLim = [0 0.5];
ax.XLim = [-0.1 2];
plot(x,y);

figure(3);
hold on;
fplot(g,[0,15]);
ax = gca;
ax.YLim = [-1 3];
ax.XLim = [-1 15];
plot(x,y);
z = fzero(g,[1,1.2]);
r = 0;
plot(z, r,'.','MarkerSize',15);
text (z-0.3,r+0.4,'1.077,0');

z = fzero(g,[1.9,2.1]);
plot(z, r,'.','MarkerSize',15);
text (z-0.3,r-0.4,'1.991,0');

z = fzero(g,[12.5,12.7]);
plot(z, r,'.','MarkerSize',15);
text (z-0.3,r+0.4,'12.625,0');

figure(4);
hold on;
fplot(h,[0,15]);
ax = gca;
ax.YLim = [-4 1];
ax.XLim = [-1 10];
plot(x,y);

z = fzero(h,[0.9,1.0]);
plot(z, r,'.','MarkerSize',15);
text (z-0.3,r+0.4,'0.971,0');