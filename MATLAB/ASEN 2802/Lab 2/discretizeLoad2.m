function [magnitude,location] = discretizeLoad2(f,N,L)
h = L/N;
x = h * (0:N);
y = f(x);
magnitude = zeros(N,1);
location = zeros(N,1);
for i = 1:N
    magnitude(i) = mean([y(i) y(i+1)])*h;
    if y(i) >= y(i+1)
        x1 = (x(i)+(h/2))*(h*y(i+1));
        x2 = (x(i)+(h/3))*((1/2)*h*(y(i)-y(i+1)));
        location(i) = (x1+x2)/magnitude(i);
    else
        x1 = (x(i)+(h/2))*(h*y(i));
        x2 = (x(i)+((2*h)/3))*((1/2)*h*(y(i+1)-y(i)));
        location(i) = (x1+x2)/magnitude(i);
    end
end
end