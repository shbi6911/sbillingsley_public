%% define constants
w0 = 0.3227;
w1 = 0.2225;
% w0 = 0.51;
% w1 = 0.17;
L = 10;
N = 4;
E = 10^7;
I = (0.5*(0.125^3))/12;

%% computation
% build load function
f = loading(w0,w1,L);
% find dicretized load
[mag,loc]=discretizeLoad2(f,N,L);
%calculate whiffletree
[spans] = whiffletree4(mag,loc);
%calculate moment equations
[V,M] = moment(w0,w1,L,mag,loc);
%calculate deflection
v = deflection(w0,w1,L,E,I,L);

%%output values
disp ("force magnitudes");
disp (mag);
disp ("locations");
disp (loc);
disp ("whiffletree dimensions");
disp (spans);
disp ("shear @loc")
disp (V);
disp ("moment @loc");
disp (M);
disp ("deflection at end")
disp (v);

function f = loading(w0,w1,L)
f = @(x) ((w1-w0)/L)*x + w0;
end

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

function [spars] = whiffletree4(magnitude,location)
matrix = [0,0,0,0,1,1;...
          0,0,1,1,0,0;...
          0,0,0,0,magnitude(3),-magnitude(4);...
          0,0,magnitude(1),-magnitude(2),0,0;...
          (magnitude(1)+magnitude(2)),-(magnitude(3)+magnitude(4)),0,0,0,0;...
          1,1,0,-1,-1,0];
b = [location(4)-location(3);...
         location(2)-location(1);...
         0;...
         0;...
         0;...
         location(3)-location(2)];
spars = inv(matrix)*b;
end

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

function v = deflection(w0,w1,L,E,I,x)
C1 = (w1-w0)/(120*L);
C2 = w0/24;
C3 = ((w0+w1)*L)/12;
C4 = (L^2)*(w0+(2*w1))/12;
C5 = 1/(E*I);
v = C5*((C1*x^5)+(C2*x^4)-(C3*x^3)+(C4*x^2));
end