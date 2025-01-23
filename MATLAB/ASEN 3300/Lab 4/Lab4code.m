%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

clear all; close all;

%% dataset 1 - 1st 100 points
%load data and plot in time domain
load("lab04_analysis_signal1.mat");
t1 = t; x1 = x; Fs1 = Fs;   N = 100;
figure(); hold on; grid on;
plot(t1(1:N),x1(1:N));
title("Signal 1 - Time Domain - 1st 100 points");
xlabel("t");
ylabel("x");
print('TimePlot1','-dpng');
hold off;

%fourier transform data
x1dft = fft(x1(1:100));    %fourier transform
x1dft = x1dft(1:N/2+1); %discard negative frequencies
psdx1 = (1/(Fs1*N))*abs(x1dft).^2;  %calculate power spectral density
psdx1(2:end-1) = 2*psdx1(2:end-1); %double to account for neg freqs
freq1 = linspace(0,Fs1/2,length(psdx1));

%plot power spectral density
figure();   hold on;    grid on;
plot(freq1,10*log10(psdx1));
title("Signal 1 - Periodogram using FFT - 1st 100 points");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency dB(Vrms^2/Hz");
print('FreqPlot1','-dpng');
hold off;

%% dataset 1 - entire vector
%fourier transform data
x1dft = fft(x1);    N = length(x1);
x1dft = x1dft(1:floor(N/2)+1); %discard negative frequencies
psdx1 = (1/(Fs1*N))*abs(x1dft).^2;  %calculate power spectral density
psdx1 (2:end-1) = 2*psdx1(2:end-1); %double to account for neg freqs
freq1 = linspace(0,Fs1/2,length(psdx1));

%plot power spectral density
figure();   hold on;    grid on;
plot(freq1,10*log10(psdx1));
title("Signal 1 - Periodogram using FFT - Entire Vector");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency dB(Vrms^2/Hz");
print('FreqPlot1_2','-dpng');
hold off;

%% dataset 2
%load data and plot in time domain
load("lab04_analysis_signal2.mat");
t2 = t; x2 = x; Fs2 = Fs;   N = length(x2);
figure(); hold on; grid on;
plot(t2(1:100),x2(1:100));
title("Signal 2 - Time Domain");
xlabel("t");
ylabel("x");
print('TimePlot2','-dpng');
hold off;

%fourier transform data
x2dft = fft(x2);    %fourier transform
x2dft = x2dft(1:floor(N/2)+1); %discard negative frequencies
psdx2 = (1/(Fs2*N))*abs(x2dft).^2;  %calculate power spectral density
psdx2 (2:end-1) = 2*psdx2(2:end-1); %double to account for neg freqs
freq2 = linspace(0,Fs2/2,length(psdx2));

%plot power spectral density
figure();   hold on;    grid on;
plot(freq2,10*log10(psdx2));
title("Signal 2 - Periodogram using FFT");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency dB(Vrms^2/Hz");
print('FreqPlot2','-dpng');
hold off;


%% dataset 3
load("lab04_analysis_signal3.mat");
t3 = t; x3 = x; Fs3 = Fs;   N = length(x3);
figure(); hold on; grid on;
plot(t3(1:100),x3(1:100));
title("Signal 3 - Time Domain");
xlabel("t");
ylabel("x");
print('TimePlot3','-dpng');
hold off;

%fourier transform data
x3dft = fft(x3);    %fourier transform
x3dft = x3dft(1:floor(N/2)+1); %discard negative frequencies
psdx3 = (1/(Fs3*N))*abs(x3dft).^2;  %calculate power spectral density
psdx3 (2:end-1) = 2*psdx3(2:end-1); %double to account for neg freqs
freq3 = linspace(0,Fs3/2,length(psdx3));

%plot power spectral density
figure();   hold on;    grid on;
plot(freq3,10*log10(psdx3));
title("Signal 3 - Periodogram using FFT");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency dB(Vrms^2/Hz");
print('FreqPlot3','-dpng');
hold off;
