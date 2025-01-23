
const = getConst(); %define constants and read in data

%[H,T_nought,delta]=slope_exp(const.x,const.y);  %find experimental slope

%calculate 
%t = linspace(0,3500,350); %range of time being evaluated
 
% prepare coefficients for the summation within U(x,t)
N = 10; %value of N

U.steel=temperature(const.data.steel.time,N,"steel");  %calculate temperature for all datasets
U.Brass29=temperature(const.data.Brass29.time,N,"Brass29");
U.Brass26=temperature(const.data.Brass26.time,N,"Brass26");
U.Alum25=temperature(const.data.Alum25.time,N,"Alum25");
U.Alum28=temperature(const.data.Alum28.time,N,"Alum28");

% %plot a sample value
% figure(2);
% hold on;
% plot(t,U.alpha6(8,:));
% title ("Thermocouple 8 across Time");
% xlabel("Time in Seconds");
% ylabel("Temperature (deg C)");
% 
% figure(3);
% hold on;
% for kk = 1:length(const.AL.alpha2)
%     plot(t,U.("alpha"+string(kk))(8,:));
% end
% title ("Thermocouple 8 across Time, various alpha values");
% xlabel("Time in Seconds");
% ylabel("Temperature (deg C)");

figure(4);
hold on;
for ll=1:8
    plot(const.data.steel.time,const.data.steel.("th"+string(ll)));
end
title("Temperature Progression of Steel");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');

figure(5);
hold on;
for ll=1:8
    plot(const.data.Brass26.time,const.data.Brass26.("th"+string(ll)));
end
title("Temperature Progression of Brass at 26 V");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');

figure(6);
hold on;
for ll=1:8
    plot(const.data.Brass29.time,const.data.Brass29.("th"+string(ll)));
end
title("Temperature Progression of Brass at 29 V");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');

figure(7);
hold on;
for ll=1:8
    plot(const.data.Alum25.time,const.data.Alum25.("th"+string(ll)));
end
title("Temperature Progression of Aluminum at 25 V");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');

figure(8);
hold on;
for ll=1:8
    plot(const.data.Alum28.time,const.data.Alum28.("th"+string(ll)));
end
title("Temperature Progression of Aluminum at 28 V");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
hold off;
legend (const.legend,'location','eastoutside');

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

%ensure t is the correct vector type for operations
if (isrow(t))==0
    t=t';
end

const=getConst();   %define constants
[H,T0,~]=slope_exp(const.x,const.data.(mat).ss);    %find slope of experimental data

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

function [tss,steadyvalue] = steadystate(data_x,data_y,window,threshold)
%Steadystate:   This function finds the time to steady state and the steady
%state value of input data
%   INPUTS:     data_x: a vector of independent variable values
%               data_y: a vector of dependent variable values
%               window: a scalar value of the size of the window over which
%               to smooth the data and find a forward difference
%               threshold: when the slope drops below this value this is
%               considered steady-state
%   The function smooths the input y-data, and then finds the slope of the
%   smoothed data.  The minimum value of slope is taken to be the steady
%   state value.

smoothed = smoothdata(data_y,'gaussian',window);    %smooth the data
slope=slopefinder(data_x,smoothed,window);      %find slope
index=find(slope < threshold);              %locate index where slope drops below threshold
tss=data_x(index(1));                   %find x-value at that index
steadyvalue=data_y(index(1));           %find y-value at that index
end

function [slope] = slopefinder(data_x,data_y,window)
%Slopefinder:   This function finds the slope of a dataset using a forward
%difference scheme.
%   INPUTS:     data_x: a vector of independent variable values
%               data_y: a vector of dependent variable values
%               window: a scalar value of the size of the window over which
%               to find a forward difference
%   Similar to the diff function, the function creates a vector of forward
%   difference values one shorter than the original vector.

slope=zeros(1,length(data_y)-1); %initialize output vector
for i=1:length(data_y)-1    %loop over input vector
    
    y0=data_y(i);           %find y-values
    if (i+window) > length(data_y)
        y1=data_y(end);
    else
        y1=data_y(i+window);
    end

    x0=data_x(i);           %find x-values
    if (i+window) > length(data_x)
        x1=data_x(end);
    else
        x1=data_x(i+window);
    end
slope(i)=(y1-y0)/(x1-x0);   %change in y over change in x
end
end

function [const] = getConst()
%getConst defines material and dimensional constants
%   INPUTS  None
%   OUTPUTS const:  a data structure containing contants

const.L = 5.875/39.37; %total length of the bar in meters
const.A = pi*(0.5/39.37)^2;

%set material constants (in metric)
const.Alum28.rho = 2810;    %densities
const.Alum25.rho = 2810;
const.Brass26.rho = 8500;
const.Brass29.rho = 8500;
const.steel.rho = 8000;

const.Alum28.cp=960;        %specific heats
const.Alum25.cp=960;
const.Brass26.cp=380;
const.Brass29.cp=380;
const.steel.cp=500;

const.Alum28.k=130;         %thermal conductivities
const.Alum25.k=130;
const.Brass26.k=115;
const.Brass29.k=115;
const.steel.k=16.2;

%calculate thermal diffusivities
const.Alum28.alpha1 = const.Alum28.k/(const.Alum28.cp*const.Alum28.rho);
const.Alum25.alpha1 = const.Alum28.alpha1;
const.Brass26.alpha1 = const.Brass26.k/(const.Brass26.cp*const.Brass26.rho);
const.Brass29.alpha1 = const.Brass26.alpha1;
const.steel.alpha1 = const.steel.k/(const.steel.cp*const.steel.rho);

%create a range of thermal diffusivities
a_delta=const.Alum28.alpha1*0.1;
const.Alum28.alpha2 = linspace((const.Alum28.alpha1-a_delta),(const.Alum28.alpha1+a_delta),11);
const.Alum25.alpha2 = const.Alum28.alpha2;
a_delta=const.Brass26.alpha1*0.1;
const.Brass26.alpha2 = linspace((const.Brass26.alpha1-a_delta),(const.Brass26.alpha1+a_delta),11);
const.Brass29.alpha2 = const.Brass26.alpha2;
a_delta=const.steel.alpha1*0.1;
const.steel.alpha2 = linspace((const.steel.alpha1-a_delta),(const.steel.alpha1+a_delta),11);

steelAll = load('Steel_21V_192mA');
Brass26All = load("Brass_26V_245mA");
Brass29All = load("Brass_29V_273mA");
Aluminum25All = load('Aluminum_25V_240mA');
Aluminum28All = load('Aluminum_28V_270mA');

grad=gradient(steelAll(:,9));
index = find(grad > 0.075,1)-1;

time=steelAll(:,1);
th1=steelAll(:,2);
th2=steelAll(:,3);
th3=steelAll(:,4);
th4=steelAll(:,5);
th5=steelAll(:,6);
th6=steelAll(:,7);
th7=steelAll(:,8);
th8=steelAll(:,9);

const.data.steel.time=time(index:length(time))-time(index);
const.data.steel.th1=th1(index:length(th1));
const.data.steel.th2=th2(index:length(th2));
const.data.steel.th3=th3(index:length(th3));
const.data.steel.th4=th4(index:length(th4));
const.data.steel.th5=th5(index:length(th5));
const.data.steel.th6=th6(index:length(th6));
const.data.steel.th7=th7(index:length(th7));
const.data.steel.th8=th8(index:length(th8));

grad=gradient(Brass26All(:,9));
index = find(grad > 0.075,1)-1;

time=Brass26All(:,1);
th1=Brass26All(:,2);
th2=Brass26All(:,3);
th3=Brass26All(:,4);
th4=Brass26All(:,5);
th5=Brass26All(:,6);
th6=Brass26All(:,7);
th7=Brass26All(:,8);
th8=Brass26All(:,9);

const.data.Brass26.time=time(index:length(time))-time(index);
const.data.Brass26.th1=th1(index:length(th1));
const.data.Brass26.th2=th2(index:length(th2));
const.data.Brass26.th3=th3(index:length(th3));
const.data.Brass26.th4=th4(index:length(th4));
const.data.Brass26.th5=th5(index:length(th5));
const.data.Brass26.th6=th6(index:length(th6));
const.data.Brass26.th7=th7(index:length(th7));
const.data.Brass26.th8=th8(index:length(th8));

grad=gradient(Brass29All(:,9));
index = find(grad > 0.075,1)-1;

time=Brass29All(:,1);
th1=Brass29All(:,2);
th2=Brass29All(:,3);
th3=Brass29All(:,4);
th4=Brass29All(:,5);
th5=Brass29All(:,6);
th6=Brass29All(:,7);
th7=Brass29All(:,8);
th8=Brass29All(:,9);

const.data.Brass29.time=time(index:length(time))-time(index);
const.data.Brass29.th1=th1(index:length(th1));
const.data.Brass29.th2=th2(index:length(th2));
const.data.Brass29.th3=th3(index:length(th3));
const.data.Brass29.th4=th4(index:length(th4));
const.data.Brass29.th5=th5(index:length(th5));
const.data.Brass29.th6=th6(index:length(th6));
const.data.Brass29.th7=th7(index:length(th7));
const.data.Brass29.th8=th8(index:length(th8));

grad=gradient(Aluminum25All(:,9));
index = find(grad > 0.075,1)-1;

time=Aluminum25All(:,1);
th1=Aluminum25All(:,2);
th2=Aluminum25All(:,3);
th3=Aluminum25All(:,4);
th4=Aluminum25All(:,5);
th5=Aluminum25All(:,6);
th6=Aluminum25All(:,7);
th7=Aluminum25All(:,8);
th8=Aluminum25All(:,9);

const.data.Alum25.time=time(index:length(time))-time(index);
const.data.Alum25.th1=th1(index:length(th1));
const.data.Alum25.th2=th2(index:length(th2));
const.data.Alum25.th3=th3(index:length(th3));
const.data.Alum25.th4=th4(index:length(th4));
const.data.Alum25.th5=th5(index:length(th5));
const.data.Alum25.th6=th6(index:length(th6));
const.data.Alum25.th7=th7(index:length(th7));
const.data.Alum25.th8=th8(index:length(th8));

grad=gradient(Aluminum28All(:,9));
index = find(grad > 0.075,1)-1;

time=Aluminum28All(:,1);
th1=Aluminum28All(:,2);
th2=Aluminum28All(:,3);
th3=Aluminum28All(:,4);
th4=Aluminum28All(:,5);
th5=Aluminum28All(:,6);
th6=Aluminum28All(:,7);
th7=Aluminum28All(:,8);
th8=Aluminum28All(:,9);

const.data.Alum28.time=time(index:length(time))-time(index);
const.data.Alum28.th1=th1(index:length(th1));
const.data.Alum28.th2=th2(index:length(th2));
const.data.Alum28.th3=th3(index:length(th3));
const.data.Alum28.th4=th4(index:length(th4));
const.data.Alum28.th5=th5(index:length(th5));
const.data.Alum28.th6=th6(index:length(th6));
const.data.Alum28.th7=th7(index:length(th7));
const.data.Alum28.th8=th8(index:length(th8));

%initialize arrays for steady-state times and values
const.data.steel.tss=zeros(1,8);
const.data.steel.ss=zeros(1,8);
const.data.Brass29.tss=zeros(1,8);
const.data.Brass29.ss=zeros(1,8);
const.data.Brass26.tss=zeros(1,8);
const.data.Brass26.ss=zeros(1,8);
const.data.Alum25.tss=zeros(1,8);
const.data.Alum25.ss=zeros(1,8);
const.data.Alum28.tss=zeros(1,8);
const.data.Alum28.ss=zeros(1,8);

%find steady-state values for all datasets
for qq=1:8
    str=string(qq);
    [const.data.steel.tss(qq),const.data.steel.ss(qq)]=steadystate(const.data.steel.time,...
        const.data.steel.("th"+str),50,0.00001);
    [const.data.Brass29.tss(qq),const.data.Brass29.ss(qq)]=steadystate(const.data.Brass29.time,...
        const.data.Brass29.("th"+str),25,0.00001);
    [const.data.Brass26.tss(qq),const.data.Brass26.ss(qq)]=steadystate(const.data.Brass26.time,...
        const.data.Brass26.("th"+str),25,0.00001);
    [const.data.Alum25.tss(qq),const.data.Alum25.ss(qq)]=steadystate(const.data.Alum25.time,...
        const.data.Alum25.("th"+str),25,0.00001);
    [const.data.Alum28.tss(qq),const.data.Alum28.ss(qq)]=steadystate(const.data.Alum28.time,...
        const.data.Alum28.("th"+str),25,0.00001);
end

const.legend=["Thermocouple 1","Thermocouple 2","Thermocouple 3",...
    "Thermocouple 4","Thermocouple 5","Thermocouple 6","Thermocouple 7",...
    "Thermocouple 8"];

const.y=[17.6 21.61 25.13 29.22 34.92 38.10 45.21 47.01]; %initial given temps
%   positions of thermocouples along the bar converted to metric
const.x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875] ./ 39.37;
end