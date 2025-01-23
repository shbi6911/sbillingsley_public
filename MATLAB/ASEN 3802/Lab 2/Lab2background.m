y=[17.6 21.61 25.61 25.13 29.22 34.92 38.10 45.21 47.01];
x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875 5.375];

hold on;
plot(x,y)

[fit, s1]=polyfit(x,y,1);

[fitline,delta] = polyval(fit,x,s1);

plot (x,fitline);
plot (x,fitline+2*delta);
plot (x,fitline-2*delta);

T_nought = polyval(fit,0);
H=fit(1);
