%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

load('lab06_section04_signal7in.mat');
t7 = t;
sig7 = ch1;
Fs7 = Fs;

load('lab06_section04_signal9in.mat');
t9 = t;
sig9 = ch1;
Fs9 = Fs;

load('lab06_section04_signal11in.mat');
t11 = t;
sig11 = ch1;
Fs11 = Fs;


%% I can't stand Matlab's default tiny font size, so make the default
% bigger:

set(0,'DefaultAxesFontSize',16);


%% to save space and limit the number of figures to submit, try plotting all
% of the signals on the same plot

% start by creating figure and axis handles

h1 = figure(1);
ax1 = axes;

plot(ax1,t7,sig7,'b');

% 'hold' ensures the next plot doesn't overwrite the first

hold(ax1,'on');

plot(ax1,t9,sig9,'r');

plot(ax1,t11,sig11,'k');

% don't forget to label the axes, and include units

xlabel(ax1,'Time (sec)');
ylabel(ax1,'Signal (volts)');
title(ax1,'Beam vibration at different lengths');

% A legend is also helpful:

legend(ax1,'length 1','length 2','length 3');


%% Alternatively, if the signals are too "busy" or overlapping too much, put them on different subplots

h2 = figure(2);
ax2(1) = subplot(311);
ax2(2) = subplot(312);
ax2(3) = subplot(313);


plot(ax2(1),t7,sig7,'b');

% we don't need to "hold on" because they'll be on different subplots

plot(ax2(2),t9,sig9,'r');

plot(ax2(3),t11,sig11,'k');

% don't forget to label the axes, and include units

xlabel(ax2(3),'Time (sec)');
ylabel(ax2(1),'Signal (volts)');
ylabel(ax2(2),'Signal (volts)');
ylabel(ax2(3),'Signal (volts)');

title(ax2(1),'Beam vibration at length 1');
title(ax2(2),'Beam vibration at length 2');
title(ax2(3),'Beam vibration at length 3');