function [const] = getConst()
%getConst defines material and dimensional constants
%   INPUTS  None
%   OUTPUTS const:  a data structure containing contants

const.L = 5.875/39.37; %total length of the bar in meters
const.A = pi*(0.5/39.37)^2;
%   positions of thermocouples along the bar converted to metric
const.x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875] ./ 39.37;

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

%find slope values for all datasets
[const.steel.H.exp,const.steel.T0]=slope_exp(const.x,const.data.steel.ss);
const.steel.H.an=slope_an(const.steel.k,const.A,21,(192/1000));

[const.Brass26.H.exp,const.Brass26.T0]=slope_exp(const.x,const.data.Brass26.ss);
const.Brass26.H.an=slope_an(const.Brass26.k,const.A,26,(245/1000));

[const.Brass29.H.exp,const.Brass29.T0]=slope_exp(const.x,const.data.Brass29.ss);
const.Brass29.H.an=slope_an(const.Brass29.k,const.A,29,(273/1000));

[const.Alum25.H.exp,const.Alum25.T0]=slope_exp(const.x,const.data.Alum25.ss);
const.Alum25.H.an=slope_an(const.Alum25.k,const.A,25,(240/1000));

[const.Alum28.H.exp,const.Alum28.T0]=slope_exp(const.x,const.data.Alum28.ss);
const.Alum28.H.an=slope_an(const.Alum28.k,const.A,28,(270/1000));

const.legend=["Thermocouple 1","Thermocouple 2","Thermocouple 3",...
    "Thermocouple 4","Thermocouple 5","Thermocouple 6","Thermocouple 7",...
    "Thermocouple 8"];

%const.y=[17.6 21.61 25.13 29.22 34.92 38.10 45.21 47.01]; %initial given temps
end