%% NACA_Aifroils Function
%Compute x_b =  [Xmin,Xmax] and y_b = [Ymin,Ymax] for NACA airfoils based
%on variables m = maximum camber in 
%

clc
clear 
close all

NACA =  input("Input NACA Designation: (XXXX)\n","s" );


m = str2double(NACA(1))/100;
p = str2double(NACA(2))/10;
t = str2double(NACA(3:4))/100;
c = input("What is the chord length: (X)\n");
N = input("How many panels do you want to use: (X)\n");

[x_b,y_b] = NACA_Airfoils(m,p,t,c,N);

writematrix([x_b,y_b],"NACA "+ NACA + " " + num2str(N))

plot(x_b,y_b,"X")
axis equal;
hold on 


function [x_b,y_b] = NACA_Airfoils(m,p,t,c,N)

yt_equation = @(t,x_c) (t.*c)./0.2 .* ((0.2969.*sqrt(x_c)) - (.1260.*x_c) - (0.3516.*(x_c.^(2))) + (0.2843.*(x_c.^(3))) - (0.1036.*(x_c.^(4))) );

%yc_equation = @(g) piecewise(0<= g &  g<= p*c, ((m.*g)./(p^2) .*(2*p - g./c)), p*c < g& g <= c, ((m.*(c - g))/((1 - p)^2) .* (1 + g./c - 2*p)))
%diff_yc = @(g) piecewise(0<= g&  g<=p*c,((-2*m)/(c*p^2))*(g-c*p),p*c < g& g <= c, ((-2*m)/(c*(p-1)^2))*(g-c*p))

theta = linspace(pi,0,N);

x = c/2.*(1-cos(theta));

x_c = x./c;

yt = yt_equation(t,x_c);
if m > 0 & p > 0
    yc = yc_function(m,c,p,x);
    yc(isnan(yc))=0;
    angle = atan(Diff_yc(m,c,p,x));
    angle(isnan(angle))=0;
    sin_angle = sin(angle);
    cos_angle = cos(angle);
elseif m == 0 & p == 0
    yc = zeros([1 N]);
    angle = atan(yc);
    sin_angle = sin(angle);
    cos_angle = cos(angle);
end


x_u = x - yt.*sin_angle
x_l = x + yt.*sin_angle;

y_u = yc + yt.*cos_angle;
y_l = yc - yt.*cos_angle;

x_b = [x_l(1:end -1)';flip(x_u)']
y_b = [y_l(1:end - 1)';flip(y_u)'];


    function [yc] = yc_function(m,c,p,x)
        %% Equation for Yc as outlined in the lab document
        r = 1;
        for g = x
            if g >=0 && g<= p*c
                yc(r) = ((m.*g)./(p^2) .*(2*p - g./c));
                fprintf("G lil: Loop" + num2str(r) +"\n")
                r = r + 1;

            elseif g > p*c && g <=c

                yc(r) = ((m.*(c - g))/((1 - p)^2) .* (1 + g./c - 2*p));
                fprintf("G big: Loop" + num2str(r) + "\n")
                r = r+1;
            end
        end
    end
    
    function [dyc] = Diff_yc(m,c,p,x)
        %% Equation for dyc/dx 
        r = 1;
        for g = x
             if g >=0 && g<= p*c
                dyc(r) = ((-2*m)/(c*p^2))*(g-c*p);
                r = r + 1;
            elseif g > p*c && g <=c
                dyc(r) = ((-2*m)/(c*(p-1)^2))*(g-c*p);
                r = r+1;
            end
        end
    end
end