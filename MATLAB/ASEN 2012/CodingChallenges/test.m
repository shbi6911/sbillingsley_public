alpha = 0.008;
beta = 0.0002;
gamma = -0.26;

f = @(x) ((gamma*x)./(1-(beta*x)))-(alpha*(x.^2));

x = linspace(-10,10,100);
y = f(x);

plot(x,y);
hold on;
yline(0);