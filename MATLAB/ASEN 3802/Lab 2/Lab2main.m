y=[17.6 21.61 25.13 29.22 34.92 38.10 45.21 47.01]; %initial given temps
%   positions of thermocouples along the bar converted to metric
x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875] ./ 39.37;
L = .14925; %total length of the bar in meters
 

figure(1);
hold on;
%plot(x,y)
 
[fit, s1]=polyfit(x,y,1);           %fit a line to given data
 
[fitline,delta] = polyval(fit,x,s1);
 
%plot (x,fitline);
%plot (x,fitline+2*delta);
%plot (x,fitline-2*delta);
%xlabel("Distance Along the Bar (m)");
%ylabel("Temperature (deg C)");

hold off;
 
T_nought = polyval(fit,0);  %T_nought is the value of the fit line at x=zero
H=fit(1);                   %H is the slope of the fit line

%set material constants (in metric)
rho = 2810;
c_p=960;
k=130;
alpha1 = k/(c_p*rho);

t = linspace(0,3500,350); %range of time being evaluated
 
% prepare coefficients for the summation within U(x,t)
N = 10; %value of N
Nvec=linspace(1,N,N);   %array of values of N
b_n=bncoeff(Nvec,H,L);  %array of values of coefficient b_n
lambda_n=(((2.*Nvec)-1).*pi)./(2*L);    %array of values of coefficient lambda_n

%initialize arrays for summation values and values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
Z_sum=zeros(length(x),length(t));
U=zeros(length(x),length(t));

for j = 1:length(x)     %loop through x-values (positions along the bar)

    X=x(j); %pull out the current value of x

    for i = 1:N     %loop through N for the summation
        
        %generate a term of the summation for the current value of n
        %(across all time values)
        Z = b_n(i) .* sin(lambda_n(i) .* X).*exp(-(lambda_n(i).^2) .* alpha1 .* t); 
        
        %add this new term to the summation
        Z_sum(j,:) = Z_sum(j,:)+Z;

    end
    %calculate U(x,t) for the current x-value across time
    U(j,:) = T_nought + (H*X) + Z_sum(j,:);
end

%plot a sample value
figure(2);
hold on;
plot(t,U(8,:));
title ("Thermocouple 8 across Time");
xlabel("Time in Seconds");
ylabel("Temperature (deg C)");


function [b_n] = bncoeff(n,H,L)
%bncoeff calculates the value of the b_n coefficient
%   INPUT:  values of n as a scalar or as an array
%           constant value of H (slope of the temperature line)
%           constant value of L (distance from T_0 to heater)
%   OUTPUT: values of b_n as a scalar or as an array

coeff1=(-2*H)/L;
coeff2=((2.*n)-1).^2;
coeff3=((-1).^(n+1));
b_n=coeff1.*coeff3.*((4*L^2)./(coeff2.*(pi^2)));

end

