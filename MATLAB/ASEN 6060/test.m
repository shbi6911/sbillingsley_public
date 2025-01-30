x = linspace(0,2*pi,100);
y = sin(x);

figure();
axis equal
plot(x,y);
arrow([x(length(x)/2),y(length(x)/2)],[x(1+length(x)/2),y(1+length(x)/2)]);