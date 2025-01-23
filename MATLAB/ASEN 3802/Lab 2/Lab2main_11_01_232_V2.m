%define reference strings for materials
steel="ST"; aluminum="AL"; brass="BR" ;

const = getConst(); %define constants

[HSteel,T_noughtSt,deltaST]=slope_exp(const.x,const.data.steel.steady);  %find experimental slope
[HAL25,T_noughtBR25,deltaAL25]=slope_exp(const.x,const.data.Alum25.steady);  %find experimental slope
[HAL28,T_noughtAL28,deltaAL28]=slope_exp(const.x,const.data.Alum28.steady);  %find experimental slope
[HBR26,T_noughtBR26,deltaBR26]=slope_exp(const.x,const.data.Brass26.steady);  %find experimental slope
[HBR29,T_noughtBR29,deltaBR29]=slope_exp(const.x,const.data.Brass29.steady);  %find experimental slope


t = linspace(0,3500,350); %range of time being evaluated
 
% prepare coefficients for the summation within U(x,t)
N = 10; %value of N

U=temperature(t,N,aluminum);  %calculate temperature

figure(1);
plot(const.x,const.data.steel.steady,".")
hold on
x = linspace(0,const.x(end));
p = [H T_nought];
y = polyval(p,x);
plot(x,y);
hold off
%plot a sample value
figure(2);
hold on;
plot(t,U.alpha6(8,:));
title ("Thermocouple 8 across Time");
xlabel("Time in Seconds");
ylabel("Temperature (deg C)");
hold off

figure(3);
hold on;
for kk = 1:length(const.AL.alpha2)
    plot(t,U.("alpha"+string(kk))(8,:));
end
title ("Thermocouple 8 across Time, various alpha values");
xlabel("Time in Seconds");
ylabel("Temperature (deg C)");
hold off











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

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.("alpha"+string(ii))=zeros(length(const.x),length(t));
end

for jj=1:length(const.(mat).alpha2)
    
    %create summation array
        Z_sum=zeros(length(const.x),length(t));

    for j = 1:length(const.x)     %loop through x-values (positions along the bar)
    
        X=const.x(j); %pull out the current value of x
        
        for i = 1:N     %loop through N for the summation
            
            %generate a term of the summation for the current value of n
            %(across all time values)
            Z = b_n(i) .* sin(lambda_n(i) .* X).*exp(-(lambda_n(i).^2) .* const.(mat).alpha2(jj) .* t); 
            
            %add this new term to the summation
            Z_sum(j,:) = Z_sum(j,:)+Z;
    
        end
        %calculate U(x,t) for the current x-value across time
        U.("alpha"+string(jj))(j,:) = T0 + (H*X) + Z_sum(j,:);
    end 
end
end

function [steady] = steadyState(th1,th2,th3,th4,th5,th6,th7,th8)
slope1 = (th1(end)-th1);
slope2 = (th2(end)-th2);
slope3 = (th3(end)-th3);
slope4 = (th4(end)-th4);
slope5 = (th5(end)-th5);
slope6 = (th6(end)-th6);
slope7 = (th7(end)-th7);
slope8 = (th8(end)-th8);
sth1 = th1(find(slope1==0,1)+1);
sth2 = th2(find(slope2==0,1)+1);
sth3 = th3(find(slope3==0,1)+1);
sth4 = th4(find(slope4==0,1)+1);
sth5 = th5(find(slope5==0,1)+1);
sth6 = th6(find(slope6==0,1)+1);
sth7 = th7(find(slope7==0,1)+1);
sth8 = th8(find(slope8==0,1)+1);
steady =[sth1 sth2 sth3 sth4 sth5 sth6 sth7 sth8];
end

function [const] = getConst()
%getConst defines material and dimensional constants
%   INPUTS  None
%   OUTPUTS const:  a data structure containing contants

const.L = 5.875/39.37; %total length of the bar in meters

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

%create a range of thermal diffusivities
a_delta=const.AL.alpha1*0.1;
const.AL.alpha2 = linspace((const.AL.alpha1-a_delta),(const.AL.alpha1+a_delta),11);
a_delta=const.BR.alpha1*0.1;
const.BR.alpha2 = linspace((const.BR.alpha1-a_delta),(const.BR.alpha1+a_delta),11);
a_delta=const.ST.alpha1*0.1;
const.ST.alpha2 = linspace((const.ST.alpha1-a_delta),(const.ST.alpha1+a_delta),11);

steelAll = load('Steel_21V_192mA');
Brass26All = load("Brass_26V_245mA");
Brass29All = load("Brass_29V_273mA");
Aluminum25All = load('Aluminum_25V_240mA');
Aluminum28All = load('Aluminum_28V_270mA');
const.data.steel.th1 = steelAll(:,2);
const.data.steel.th2 = steelAll(:,3);
const.data.steel.th3 = steelAll(:,4);
const.data.steel.th4 = steelAll(:,5);
const.data.steel.th5 = steelAll(:,6);
const.data.steel.th6 = steelAll(:,7);
const.data.steel.th7 = steelAll(:,8);
const.data.steel.th8 = steelAll(:,9);

const.data.Brass26.th1 = Brass26All(:,2);
const.data.Brass26.th2 = Brass26All(:,3);
const.data.Brass26.th3 = Brass26All(:,4);
const.data.Brass26.th4 = Brass26All(:,5);
const.data.Brass26.th5 = Brass26All(:,6);
const.data.Brass26.th6 = Brass26All(:,7);
const.data.Brass26.th7 = Brass26All(:,8);
const.data.Brass26.th8 = Brass26All(:,9);


const.data.Brass29.th1 = Brass29All(:,2);
const.data.Brass29.th2 = Brass29All(:,3);
const.data.Brass29.th3 = Brass29All(:,4);
const.data.Brass29.th4 = Brass29All(:,5);
const.data.Brass29.th5 = Brass29All(:,6);
const.data.Brass29.th6 = Brass29All(:,7);
const.data.Brass29.th7 = Brass29All(:,8);
const.data.Brass29.th8 = Brass29All(:,9);

const.data.Alum25.th1 = Aluminum25All(:,2);
const.data.Alum25.th2 = Aluminum25All(:,3);
const.data.Alum25.th3 = Aluminum25All(:,4);
const.data.Alum25.th4 = Aluminum25All(:,5);
const.data.Alum25.th5 = Aluminum25All(:,6);
const.data.Alum25.th6 = Aluminum25All(:,7);
const.data.Alum25.th7 = Aluminum25All(:,8);
const.data.Alum25.th8 = Aluminum25All(:,9);


const.data.Alum28.th1 = Aluminum28All((200:end),2);
const.data.Alum28.th2 = Aluminum28All((200:end),3);
const.data.Alum28.th3 = Aluminum28All((200:end),4);
const.data.Alum28.th4 = Aluminum28All((200:end),5);
const.data.Alum28.th5 = Aluminum28All((200:end),6);
const.data.Alum28.th6 = Aluminum28All((200:end),7);
const.data.Alum28.th7 = Aluminum28All((200:end),8);
const.data.Alum28.th8 = Aluminum28All((200:end),9);


const.data.steel.steady =steadyState(const.data.steel.th1,const.data.steel.th2,const.data.steel.th3,const.data.steel.th4,const.data.steel.th5,const.data.steel.th6,const.data.steel.th7,const.data.steel.th8);

const.data.Brass26.steady =steadyState(const.data.Brass26.th1,const.data.Brass26.th2,const.data.Brass26.th3,const.data.Brass26.th4,const.data.Brass26.th5,const.data.Brass26.th6,const.data.Brass26.th7,const.data.Brass26.th8);

const.data.Brass29.steady =steadyState(const.data.Brass29.th1,const.data.Brass29.th2,const.data.Brass29.th3,const.data.Brass29.th4,const.data.Brass29.th5,const.data.Brass29.th6,const.data.Brass29.th7,const.data.Brass29.th8);


const.data.Alum25.steady =steadyState(const.data.Alum25.th1,const.data.Alum25.th2,const.data.Alum25.th3,const.data.Alum25.th4,const.data.Alum25.th5,const.data.Alum25.th6,const.data.Alum25.th7,const.data.Alum25.th8);


const.data.Alum28.steady =steadyState(const.data.Alum28.th1,const.data.Alum28.th2,const.data.Alum28.th3,const.data.Alum28.th4,const.data.Alum28.th5,const.data.Alum28.th6,const.data.Alum28.th7,const.data.Alum28.th8);

const.y=[17.6 21.61 25.13 29.22 34.92 38.10 45.21 47.01]; %initial given temps
%   positions of thermocouples along the bar converted to metric
const.x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875] ./ 39.37;
end