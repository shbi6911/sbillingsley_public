%% Main Code

const = getConst(); %define constants and read in data

% designate N-value for Fourier series summation
N = 10; %value of N

%calculate temperature for all datasets for models 1A and 1B
U.steel=temperature(const.data.steel.time,N,"steel");
U.Brass29=temperature(const.data.Brass29.time,N,"Brass29");
U.Brass26=temperature(const.data.Brass26.time,N,"Brass26");
U.Alum25=temperature(const.data.Alum25.time,N,"Alum25");
U.Alum28=temperature(const.data.Alum28.time,N,"Alum28");

%calculate initial state distributions from experimental data
fitlines=initialstate();

%calculate temperature for all datasets using Model II
U2.steel=temperature_model2(const.data.steel.time,N,"steel",fitlines);  
U2.Brass29=temperature_model2(const.data.Brass29.time,N,"Brass29",fitlines);
U2.Brass26=temperature_model2(const.data.Brass26.time,N,"Brass26",fitlines);
U2.Alum25=temperature_model2(const.data.Alum25.time,N,"Alum25",fitlines);
U2.Alum28=temperature_model2(const.data.Alum28.time,N,"Alum28",fitlines);


%% Plot Analytical Slope vs. Experimental Data (Model 1A)
figure(4);
%figure('units','inch','position',[1,1,4.5,4]);
hold on;

for ll=1:8
    h1 = plot(const.data.steel.time,const.data.steel.("th"+string(ll)),'r');
    h2 = plot(const.data.steel.time,U.steel.an.alpha6(ll,:),'b');
end
title("Temperature Progression of Steel (Model 1A Analytical)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1A"],'location','southeast');
hold off;
%print('Steel_Model1A','-dpng','-r300');

figure(5);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass26.time,const.data.Brass26.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass26.time,U.Brass26.an.alpha6(ll,:),'b');
end
title("Temperature Progression of Brass at 26 V (Model 1A Analytical");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1A"],'location','southeast');
hold off;
%print('Brass26_Model1A','-dpng','-r300');


figure(6);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass29.time,const.data.Brass29.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass29.time,U.Brass29.an.alpha6(ll,:),'b');
end
title("Temperature Progression of Brass at 29 V (Model 1A (Analytical)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1A"],'location','southeast');
hold off;
%print('Brass29_Model1A','-dpng','-r300');


figure(7);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum25.time,const.data.Alum25.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum25.time,U.Alum25.an.alpha6(ll,:),'b');
end
title("Temperature Progression of Aluminum at 25 V (Model 1A Analytical)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1A"],'location','southeast');
hold off;
%print('Alum25_Model1A','-dpng','-r300');


figure(8);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum28.time,const.data.Alum28.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum28.time,U.Alum28.an.alpha6(ll,:),'b');
end
title("Temperature Progression of Aluminum at 28 V (Model 1A Analytical)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1A"],'location','southeast');
hold off;
%print('Alum28_Model1A','-dpng','-r300');

%% Plot Experimental Slope vs. Experimental Data (Model 1B)
figure(9);
hold on;
for ll=1:8
    h1 = plot(const.data.steel.time,const.data.steel.("th"+string(ll)),'r');
    h2 = plot(const.data.steel.time,U.steel.exp.alpha6(ll,:),'k');
end
title("Temperature Progression of Steel (Model 1B Experimental)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1B"],'location','southeast');
hold off;
%print('Steel_Model1B','-dpng','-r300');


figure(10);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass26.time,const.data.Brass26.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass26.time,U.Brass26.exp.alpha6(ll,:),'k');
end
title("Temperature Progression of Brass at 26 V (Model 1B Experimental");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1B"],'location','southeast');
hold off;
%print('Brass26_Model1B','-dpng','-r300');

figure(11);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass29.time,const.data.Brass29.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass29.time,U.Brass29.exp.alpha6(ll,:),'k');
end
title("Temperature Progression of Brass at 29 V (Model 1B (Experimental)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1B"],'location','southeast');
hold off;
%print('Brass29_Model1B','-dpng','-r300');

figure(12);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum25.time,const.data.Alum25.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum25.time,U.Alum25.exp.alpha6(ll,:),'k');
end
title("Temperature Progression of Aluminum at 25 V (Model 1B Experimental)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1B"],'location','southeast');
hold off;
%print('Alum25_Model1B','-dpng','-r300');

figure(13);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum28.time,const.data.Alum28.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum28.time,U.Alum28.exp.alpha6(ll,:),'k');
end
title("Temperature Progression of Aluminum at 28 V (Model 1B Experimental)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model 1B"],'location','southeast');
hold off;
%print('Alum28_Model1B','-dpng','-r300');

%% Plot Results of Model II

figure(19);
hold on;
for ll=1:8
    h1 = plot(const.data.steel.time,const.data.steel.("th"+string(ll)),'r');
    h2 = plot(const.data.steel.time,U2.steel.exp.model2.alpha6(ll,:),'k');
end
title("Temperature Progression of Steel (Model II)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model II"],'location','southeast');
hold off;
%print('Steel_Model2','-dpng','-r300');

figure(20);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass26.time,const.data.Brass26.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass26.time,U2.Brass26.exp.model2.alpha6(ll,:),'k');
end
title("Temperature Progression of Brass at 26 V (Model II");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model II"],'location','southeast');
hold off;
%print('Brass26_Model2','-dpng','-r300');

figure(21);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass29.time,const.data.Brass29.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass29.time,U2.Brass29.exp.model2.alpha6(ll,:),'k');
end
title("Temperature Progression of Brass at 29 V (Model II");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model II"],'location','southeast');
hold off;
%print('Brass29_Model2','-dpng','-r300');

figure(22);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum25.time,const.data.Alum25.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum25.time,U2.Alum25.exp.model2.alpha6(ll,:),'k');
end
title("Temperature Progression of Aluminum at 25 V (Model II");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model II"],'location','southeast');
hold off;
%print('Alum25_Model2','-dpng','-r300');

figure(23);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum28.time,const.data.Alum28.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum28.time,U2.Alum28.exp.model2.alpha6(ll,:),'k');
end
title("Temperature Progression of Aluminum at 28 V (Model II)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model II"],'location','southeast');
hold off;
%print('Alum28_Model2','-dpng','-r300');

%% Plot Results of Model III
figure(24);
hold on;
for ll=1:8
    h1 = plot(const.data.steel.time,const.data.steel.("th"+string(ll)),'r');
    h2 = plot(const.data.steel.time,U2.steel.exp.model2.alpha5(ll,:),'k');
end
title("Temperature Progression of Steel (Model III)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model III"],'location','southeast');
hold off;
%print('Steel_Model3','-dpng','-r300');

figure(25);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass26.time,const.data.Brass26.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass26.time,U2.Brass26.exp.model2.alpha2(ll,:),'k');
end
title("Temperature Progression of Brass at 26 V (Model III");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model III"],'location','southeast');
hold off;
%print('Brass26_Model3','-dpng','-r300');

figure(26);
hold on;
for ll=1:8
    h1 = plot(const.data.Brass29.time,const.data.Brass29.("th"+string(ll)),'r');
    h2 = plot(const.data.Brass29.time,U2.Brass29.exp.model2.alpha2(ll,:),'k');
end
title("Temperature Progression of Brass at 29 V (Model III");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model III"],'location','southeast');
hold off;
%print('Brass29_Model3','-dpng','-r300');

figure(27);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum25.time,const.data.Alum25.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum25.time,U2.Alum25.exp.model2.alpha3(ll,:),'k');
end
title("Temperature Progression of Aluminum at 25 V (Model III");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model III"],'location','southeast');
hold off;
%print('Alum25_Model3','-dpng','-r300');

figure(28);
hold on;
for ll=1:8
    h1 = plot(const.data.Alum28.time,const.data.Alum28.("th"+string(ll)),'r');
    h2 = plot(const.data.Alum28.time,U2.Alum28.exp.model2.alpha3(ll,:),'k');
end
title("Temperature Progression of Aluminum at 28 V (Model III)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ([h1(1) h2(1)],["Experimental Data","Model III"],'location','southeast');
hold off;
%print('Alum28_Model3','-dpng','-r300');

%% Plot Steady State Values

figure(29);
hold on;
scatter(const.x,const.data.steel.ss);
fit = [const.steel.H.exp const.steel.T0];
plot(const.x,polyval(fit,const.x), 'Linewidth' , 2);
title("Steady State Values of Steel");
xlabel ("Position Along Rod (m)");
ylabel ("Temperature (degrees C)");
legend ("Steady State Data","Linear Fit Line",'location','southeast');
hold off;
%print('Steel_Steady','-dpng','-r300');

figure(30);
hold on;
scatter(const.x,const.data.Brass29.ss);
fit = [const.Brass29.H.exp const.Brass29.T0];
plot(const.x,polyval(fit,const.x), 'Linewidth' , 2);
title("Steady State Values of Brass @ 29 V");
xlabel ("Position Along Rod (m)");
ylabel ("Temperature (degrees C)");
legend ("Steady State Data","Linear Fit Line",'location','southeast');
hold off;
%print('Brass29_Steady','-dpng','-r300');

figure(31);
hold on;
scatter(const.x,const.data.Brass26.ss);
fit = [const.Brass26.H.exp const.Brass26.T0];
plot(const.x,polyval(fit,const.x), 'Linewidth' , 2);
title("Steady State Values of Brass @ 26 V");
xlabel ("Position Along Rod (m)");
ylabel ("Temperature (degrees C)");
legend ("Steady State Data","Linear Fit Line",'location','southeast');
hold off;
%print('Brass29_Steady','-dpng','-r300');

figure(32);
hold on;
scatter(const.x,const.data.Alum28.ss);
fit = [const.Alum28.H.exp const.Alum28.T0];
plot(const.x,polyval(fit,const.x), 'Linewidth' , 2);
title("Steady State Values of Aluminum @ 28 V");
xlabel ("Position Along Rod (m)");
ylabel ("Temperature (degrees C)");
legend ("Steady State Data","Linear Fit Line",'location','southeast');
hold off;
%print('Alum28_Steady','-dpng','-r300');

figure(33);
hold on;
scatter(const.x,const.data.Alum25.ss);
fit = [const.Alum25.H.exp const.Alum25.T0];
plot(const.x,polyval(fit,const.x), 'Linewidth' , 2);
title("Steady State Values of Aluminum @ 25 V");
xlabel ("Position Along Rod (m)");
ylabel ("Temperature (degrees C)");
legend ("Steady State Data","Linear Fit Line",'location','southeast');
hold off;
%print('Alum28_Steady','-dpng','-r300');

%% Plot Error Bar plots for Thermocouple 8

figure(34);
hold on;
err=NaN(1,length(const.data.steel.th8));
err(1:25:length(const.data.steel.th8))=2;
errorbar(const.data.steel.time,const.data.steel.th8,err,'k');
plot(const.data.steel.time,U2.steel.exp.model2.alpha6(8,:),'b');
plot(const.data.steel.time,U2.steel.exp.model2.alpha5(8,:),'r');

title("Temperature Progression of Steel (Thermocouple 8)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ("Experimental Data w/Errorbars","Model II","Model III",'location','southeast');
hold off;
%print('Steel_Error','-dpng','-r300');

figure(35);
hold on;
err=NaN(1,length(const.data.Brass29.th8));
err(1:25:length(const.data.Brass29.th8))=2;
errorbar(const.data.Brass29.time,const.data.Brass29.th8,err,'k');
plot(const.data.Brass29.time,U2.Brass29.exp.model2.alpha6(8,:),'b');
plot(const.data.Brass29.time,U2.Brass29.exp.model2.alpha2(8,:),'r');

title("Temperature Progression of Brass @ 29 V (Thermocouple 8)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ("Experimental Data w/Errorbars","Model II","Model III",'location','southeast');
hold off;
%print('Brass29_Error','-dpng','-r300');

figure(36);
hold on;
err=NaN(1,length(const.data.Brass26.th8));
err(1:25:length(const.data.Brass26.th8))=2;
errorbar(const.data.Brass26.time,const.data.Brass26.th8,err,'k');
plot(const.data.Brass26.time,U2.Brass26.exp.model2.alpha6(8,:),'b');
plot(const.data.Brass26.time,U2.Brass26.exp.model2.alpha2(8,:),'r');

title("Temperature Progression of Brass @ 26 V (Thermocouple 8)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ("Experimental Data w/Errorbars","Model II","Model III",'location','southeast');
hold off;
%print('Brass26_Error','-dpng','-r300');

figure(37);
hold on;
err=NaN(1,length(const.data.Alum28.th8));
err(1:25:length(const.data.Alum28.th8))=2;
errorbar(const.data.Alum28.time,const.data.Alum28.th8,err,'k');
plot(const.data.Alum28.time,U2.Alum28.exp.model2.alpha6(8,:),'b');
plot(const.data.Alum28.time,U2.Alum28.exp.model2.alpha3(8,:),'r');

title("Temperature Progression of Aluminum @ 28 V (Thermocouple 8)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ("Experimental Data w/Errorbars","Model II","Model III",'location','southeast');
hold off;
%print('Alum28_Error','-dpng','-r300');

figure(38);
hold on;
err=NaN(1,length(const.data.Alum25.th8));
err(1:25:length(const.data.Alum25.th8))=2;
errorbar(const.data.Alum25.time,const.data.Alum25.th8,err,'k');
plot(const.data.Alum25.time,U2.Alum25.exp.model2.alpha6(8,:),'b');
plot(const.data.Alum25.time,U2.Alum25.exp.model2.alpha3(8,:),'r');

title("Temperature Progression of Aluminum @ 25 V (Thermocouple 8)");
xlabel ("Time (s)");
ylabel ("Temperature (degrees C)");
legend ("Experimental Data w/Errorbars","Model II","Model III",'location','southeast');
hold off;
%print('Alum25_Error','-dpng','-r300');


%% Generate Fourier numbers
fourier.steel.tss = mean(const.data.steel.tss);
fourier.Alum28.tss = mean(const.data.Alum28.tss);
fourier.Alum25.tss = mean(const.data.Alum25.tss);
fourier.Brass29.tss = mean(const.data.Brass29.tss);
fourier.Brass26.tss = mean(const.data.Brass26.tss);
fourier.steel.alpha = const.steel.alpha2(5);
fourier.Alum28.alpha = const.Alum28.alpha2(3);
fourier.Alum25.alpha = const.Alum25.alpha2(3);
fourier.Brass29.alpha = const.Brass29.alpha2(2);
fourier.Brass26.alpha = const.Brass26.alpha2(2);
fourier.steel.F0 = (fourier.steel.tss*fourier.steel.alpha)/const.L^2;
fourier.Alum28.F0 = (fourier.Alum28.tss*fourier.Alum28.alpha)/const.L^2;
fourier.Alum25.F0 = (fourier.Alum25.tss*fourier.Alum25.alpha)/const.L^2;
fourier.Brass29.F0 = (fourier.Brass29.tss*fourier.Brass29.alpha)/const.L^2;
fourier.Brass26.F0 = (fourier.Brass26.tss*fourier.Brass26.alpha)/const.L^2;

%% functions

function [fitlines]=initialstate()
%initialstate calculates and plots the initial temperature distribution of
%the bar, and outputs the fit line coefficients
%   INPUT:  None
%   OUTPUT: fitlines is a matrix of five linear best fit vectors

%load in data set
data1 = load('Aluminum_25V_240mA');
%remove the first column - time
data1new = data1(:, 2:end);
%make an array of the first row - inital temp
TempAL1 = data1new(1,:);
%x axis - position of thermocouple
x=[1.375 1.875 2.375 2.875 3.375 3.875 4.375 4.875]./ 39.37 ; 

% Get coefficients of a line fit through the data
coefficients1 = polyfit(x, TempAL1, 1);
% Create a new x axis with exactly 8 points 
xFit1 = linspace(min(x), max(x), 8);
% Get the estimated yFit value for each of those 8 new x locations
yFit1 = polyval(coefficients1 , xFit1);

figure(14)
%dot plot of experiment data
scatter(x,TempAL1);
hold on
%plot line of best fit
plot(xFit1,yFit1,'LineWidth',2); 
title("Initial State Distribution of Aluminum at 25V")
xlabel("Location of Thermocouple (Meters)")
ylabel("Initial Temperature (Degrees C)")
legend("Initial State Data","Linear Fit Line",'location','southeast');
%re fit plot
ylim([10.5,13.5]); 
%print('Alum25_Initial','-dpng','-r300');

%load in data set
data2 = load('Aluminum_28V_270mA');
%remove the first column - temprature
data2new = data2(:, 2:end);
%make an array of the first row - initial temprature 
TempAL2 = data2new(205,:);

% Get coefficients of a line fit through the data
coefficients2 = polyfit(x, TempAL2, 1);
% Create a new x axis with exactly 8 points 
xFit2 = linspace(min(x), max(x), 8);
% Get the estimated yFit value for each of those 8 new x locations
yFit2 = polyval(coefficients2 , xFit2);

figure(15)
%dot plot of experiment data
scatter(x,TempAL2);
hold on
%plot line of best fit
plot(xFit2,yFit2,'LineWidth',2); 
title("Initial State Distribution of Aluminum at 28V")
xlabel("Location of Thermocouple (Meters)")
ylabel("Initial Temperature (Degrees C)")
legend("Initial State Data","Linear Fit Line",'location','southeast');
%re fit plot
ylim([8,11]); 
%print('Alum28_Initial','-dpng','-r300');

%load in data set 
data3 = load('Brass_26V_245mA');
%remove the first column - time
data3new = data3(:, 2:end);
%make an array of the first row - initial temp
TempBR1 = data3new(1,:);

% Get coefficients of a line fit through the data
coefficients3 = polyfit(x, TempBR1, 1);
% Create a new x axis with exactly 8 points 
xFit3 = linspace(min(x), max(x), 8);
% Get the estimated yFit value for each of those 8 new x locations
yFit3 = polyval(coefficients3 , xFit3);

figure(16)
%dot plot of experiment data
scatter(x,TempBR1);
hold on
%plot line of best fit
plot(xFit3,yFit3,'LineWidth',2); 
title("Initial State Distribution of Brass at 26V")
xlabel("Location of Thermocouple (Meters)")
ylabel("Initial Temperature (Degrees C)")
legend("Initial State Data","Linear Fit Line",'location','southeast');
%re fit plot
ylim([11.8,12.8]); 
%print('Brass26_Initial','-dpng','-r300');

%load in data set
data4 = load('Brass_29V_273mA');
%remove the first column - time
data4new = data4(:, 2:end);
%make an array of the first row - initial temp
TempBR2 = data4new(1,:);

% Get coefficients of a line fit through the data
coefficients4 = polyfit(x, TempBR2, 1);
% Create a new x axis with exactly 8 points 
xFit4 = linspace(min(x), max(x), 8);
% Get the estimated yFit value for each of those 8 new x locations
yFit4 = polyval(coefficients4 , xFit4);

figure(17)
%dot plot of experiment data
scatter(x,TempBR2);
hold on
%plot line of best fit
plot(xFit4,yFit4, 'LineWidth',2); 
title("Initial State Distribution of Brass at 29V")
xlabel("Location of Thermocouple (Meters)")
ylabel("Initial Temperature (Degrees C)")
legend("Initial State Data","Linear Fit Line",'location','southeast');
%re fit plot
ylim([11.6,12.6]); 
%print('Brass29_Initial','-dpng','-r300');

%load in data set
data5 = load('Steel_21V_192mA');
%remove the first column - time
data5new = data5(:, 2:end);
%make an array of the first row - inital temp
TempSteel = data5new(1,:);


% Get coefficients of a line fit through the data
coefficients5 = polyfit(x, TempSteel, 1);
% Create a new x axis with exactly 8 points 
xFit5 = linspace(min(x), max(x), 8);
% Get the estimated yFit value for each of those 8 new x locations
yFit5 = polyval(coefficients5 , xFit4);

figure(18)
%dot plot of experiment data
scatter(x,TempSteel);
hold on
%plot line of best fit
plot(xFit5,yFit5, 'LineWidth',2); 
title("Initial State Distribution of Steel")
xlabel("Location of Thermocouple (Meters)")
ylabel("Initial Temperature (Degrees C)")
legend("Initial State Data","Linear Fit Line",'location','southeast');
%re fit plot
ylim([13,17]);
%print('Steel_Initial','-dpng','-r300');

fitlines.Alum25=coefficients1;
fitlines.Alum28=coefficients2;
fitlines.Brass26=coefficients3;
fitlines.Brass29=coefficients4;
fitlines.steel=coefficients5;

end

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

function [b_n] = bncoeff2(n,H,L,M)
%bncoeff calculates the value of the b_n coefficient
%   INPUT:  values of n as a scalar or as an array
%           constant value of H (slope of the temperature line)
%           constant value of L (distance from T_0 to heater)
%           scalar initial state linear fit slope M
%   OUTPUT: values of b_n as a scalar or as a vector

coeff1=(2*(M-H))/L;
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
%[H,T0,~]=slope_exp(const.x,const.data.(mat).ss);    %find slope of experimental data

Nvec=1:N;   %array of values of N
b_n=bncoeff(Nvec,const.(mat).H.exp,const.L);  %array of values of coefficient b_n
lambda_n=(((2.*Nvec)-1).*pi)./(2*const.L);    %array of values of coefficient lambda_n

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.exp.("alpha"+string(ii))=zeros(length(const.x),length(t));
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
        U.exp.("alpha"+string(jj))(j,:) = const.(mat).T0 + (const.(mat).H.exp*X) + Z_sum(j,:);
    end 
end

b_n=bncoeff(Nvec,const.(mat).H.an,const.L);  %array of values of coefficient b_n

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.an.("alpha"+string(ii))=zeros(length(const.x),length(t));
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
        U.an.("alpha"+string(jj))(j,:) = const.(mat).T0 + (const.(mat).H.an*X) + Z_sum(j,:);
    end 
end

end

function [U] = temperature_model2(t,N,mat,fit)
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

Nvec=1:N;   %array of values of N
b_n=bncoeff2(Nvec,const.(mat).H.exp,const.L,fit.(mat)(1));  %array of values of coefficient b_n
lambda_n=(((2.*Nvec)-1).*pi)./(2*const.L);    %array of values of coefficient lambda_n

%initialize arrays for values of U(x,t)
%these are 2D arrays with rows corresponding to x and columns to t
%U is a structure with such a 2D array for each alpha value

for ii = 1:length(const.(mat).alpha2)
    U.exp.model2.("alpha"+string(ii))=zeros(length(const.x),length(t));
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
        U.exp.model2.("alpha"+string(jj))(j,:) = fit.(mat)(2) + (const.(mat).H.exp*X) + Z_sum(j,:);
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
%getConst defines material and dimensional constants and reads in data
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
a_delta=const.Alum28.alpha1*0.5;
const.Alum28.alpha2 = linspace((const.Alum28.alpha1-a_delta),(const.Alum28.alpha1+a_delta),11);
const.Alum25.alpha2 = const.Alum28.alpha2;
a_delta=const.Brass26.alpha1*0.7;
const.Brass26.alpha2 = linspace((const.Brass26.alpha1-a_delta),(const.Brass26.alpha1+a_delta),11);
const.Brass29.alpha2 = const.Brass26.alpha2;
a_delta=const.steel.alpha1*0.5;
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