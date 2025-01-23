steel = "ST";   aluminum="AL";  brass="BR";

[H,T_nought,delta]=slope_exp(x,y);
const = getConst();

a_delta=const.(aluminum).alpha1*0.1;
alpha2_AL = linspace((const.AL.alpha1-a_delta),(const.AL.alpha1+a_delta),10);

t = linspace(0,3500,350); %range of time being evaluated
 
% prepare coefficients for the summation within U(x,t)
N = 10; %value of N

U=temperature(t,N,aluminum);

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

function [H_exp,T0_exp,H_exp_delta] = slope_exp(x,y)
%slope_exp      This function calculates slope of heat transfer from
%experimental data
%INPUTS     x   :a vector of x-values
%           y   :a vector of temperature values
%OUTPUTS    H_exp       :slope of the polyfit line between x and y
%           T0_exp      :y-intercept of the polyfit line
%           H_exp_delta :error struct output by polyfit
[fit,s1]= polyfit(x,y,1);
H_exp=fit(1);
T0_exp=fit(2);
H_exp_delta=s1;
end

function [H_an] = slope_an(k,A,V,I)
%slope_an   This function calculates an analytical slope value for heat
%transfer
%INPUTS     k   :thermal conductivity
%           A   :cross-sectional area
%           V   :voltage
%           I   :current
%OUTPUTS    H_an    :scalar value of analytical slope
P=I*V;
H_an=P/(k*A);
end

function [U] = temperature(t,N,mat)
%temperature This function calculates temperature in a bar across position
%and time
%INPUTS     t   :a 1D vector of time values
%           N   :a scalar value for iterations of a Fourier series
%           mat :a string designating which material is being used
%               : AL=aluminum, BR=brass, ST=steel
%OUTPUTS    U   : a 2D array of temperature values ordered by x-position on
%the first direction and time value on the second

const=getConst();   %define constants
[H,T0,~]=slope_exp(const.x,const.y);    %find slope of experimental data

Nvec=1:N;   %array of values of N
b_n=bncoeff(Nvec,H,const.L);  %array of values of coefficient b_n
lambda_n=(((2.*Nvec)-1).*pi)./(2*const.L);    %array of values of coefficient lambda_n

%initialize arrays for summation values and values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
Z_sum=zeros(length(const.x),length(t));
U=zeros(length(const.x),length(t));

for j = 1:length(const.x)     %loop through x-values (positions along the bar)

    X=const.x(j); %pull out the current value of x

    for i = 1:N     %loop through N for the summation
        
        %generate a term of the summation for the current value of n
        %(across all time values)
        Z = b_n(i) .* sin(lambda_n(i) .* X).*exp(-(lambda_n(i).^2) .* const.(mat).alpha1 .* t); 
        
        %add this new term to the summation
        Z_sum(j,:) = Z_sum(j,:)+Z;

    end
    %calculate U(x,t) for the current x-value across time
    U(j,:) = T0 + (H*X) + Z_sum(j,:);
end

end

function [const] = getConst()
%getConst defines material and dimensional constants
%   INPUTS  None
%   OUTPUTS const:  a data structure containing contants

const.L = .14925; %total length of the bar in meters

%set material constants (in metric)
const.AL.rho = 2810;    %densities
const.BR.rho = 8500;
const.ST.rho = 8000;

const.AL.cp=960;        %specific heats
const.BR.cp=380;
const.ST.cp=500;

const.AL.k=130;         %thermal conductivities
const.BR.k=115;
const.ST.k=16.2;

%calculate thermal diffusivities
const.AL.alpha1 = const.AL.k/(const.AL.cp*const.AL.rho);
const.BR.alpha1 = const.BR.k/(const.BR.cp*const.BR.rho);
const.ST.alpha1 = const.ST.k/(const.ST.cp*const.ST.rho);

const.y=[17.6 21.61 25.13 29.22 34.92 38.10 45.21 47.01]; %initial given temps
%   positions of thermocouples along the bar converted to metric
const.x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875] ./ 39.37;
end