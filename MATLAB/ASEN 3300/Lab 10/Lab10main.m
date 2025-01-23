%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

%% Lab 10 signal AM50

%read in data
data = readmatrix("lab10_section03_signalAM50.csv");
time = data(:,1);
amplitude = data(:,2);
Fs = 1/mean(diff(time));  %sampling frequency
L = length(amplitude);

%plot time series
% scatter(time,amplitude); hold on; grid on;
% title("Time Series Signal AM 50");
% xlabel ("Time (s)");
% ylabel ("Amplitude (V)");
% hold off;

 N=length(amplitude);
 data_fft = fft(amplitude);   %run Fast Fourier Transform

 %data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs

 psdx1 = (1/(1000*N))*abs(data_fft).^2;
 psdx1 (2:end-1) = 2*psdx1(2:end-1);

 freq1 = Fs/L*(0:L-1);
 signal = find(max(psdx1));

 %plot power spectral density
 figure();   hold on;    grid on;
 plot(freq1,10*log10(psdx1));
 title("Signal AM 50 - Periodogram using FFT");
 xlabel("Frequency (Hz)");
 ylabel("Power/Frequency dB(Vrms^2/Hz");
%print('SignalAM50','-dpng');
 hold off;

%%  Lab 10 signal AM100

%read in data
data = readmatrix("lab10_section03_signalAM100.csv");
time = data(:,1);
amplitude = data(:,2);
Fs = 1/mean(diff(time));  %sampling frequency
L = length(amplitude);

%plot time series
% scatter(time,amplitude); hold on; grid on;
% title("Time Series Signal AM 100");
% xlabel ("Time (s)");
% ylabel ("Amplitude (V)");
% hold off;

 N=length(amplitude);
 data_fft = fft(amplitude);   %run Fast Fourier Transform

 %data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs

 psdx1 = (1/(1000*N))*abs(data_fft).^2;
 psdx1 (2:end-1) = 2*psdx1(2:end-1);

 freq1 = Fs/L*(0:L-1);
 signal = find(max(psdx1));

 %plot power spectral density
 figure();   hold on;    grid on;
 plot(freq1,10*log10(psdx1));
 title("Signal AM 100 - Periodogram using FFT");
 xlabel("Frequency (Hz)");
 ylabel("Power/Frequency dB(Vrms^2/Hz");
%print('FreqPlot1_2','-dpng');
 hold off;

%% Lab 10 signal AM1002kHz
% 
% %read in data
% data = readmatrix("lab10_section03_signalAM100_2kHz.csv");
% time = data(:,1);
% amplitude = data(:,2);
% Fs = 1/mean(diff(time));  %sampling frequency
% L = length(amplitude);
% 
% %plot time series
% scatter(time,amplitude); hold on; grid on;
% title("Time Series Signal AM 100");
% xlabel ("Time (s)");
% ylabel ("Amplitude (V)");
% hold off;
% 
%  N=length(amplitude);
%  data_fft = fft(amplitude);   %run Fast Fourier Transform
% 
%  %data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs
% 
%  psdx1 = (1/(1000*N))*abs(data_fft).^2;
%  psdx1 (2:end-1) = 2*psdx1(2:end-1);
% 
%  freq1 = Fs/L*(0:L-1);
%  signal = find(max(psdx1));
% 
%  %plot power spectral density
%  figure();   hold on;    grid on;
%  plot(freq1,10*log10(psdx1));
%  title("Signal AM 100-2kHz - Periodogram using FFT");
%  xlabel("Frequency (Hz)");
%  ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
%  hold off;

%  %% Lab 10 Signal FM
% 
%  %read in data
% data = readmatrix("lab10_section04_signalFM.csv");
% time = data(:,1);
% amplitude = data(:,2);
% Fs = 1/mean(diff(time));  %sampling frequency
% L = length(amplitude);
% 
% %plot time series
% % scatter(time,amplitude); hold on; grid on;
% % title("Time Series Signal AM 100");
% % xlabel ("Time (s)");
% % ylabel ("Amplitude (V)");
% % hold off;
% 
%  N=length(amplitude);
%  data_fft = fft(amplitude);   %run Fast Fourier Transform
% 
%  %data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs
% 
%  psdx1 = (1/(1000*N))*abs(data_fft).^2;
%  psdx1 (2:end-1) = 2*psdx1(2:end-1);
% 
%  freq1 = Fs/L*(0:L-1);
%  %signal = find(max(psdx1));
% 
%  %plot power spectral density
%  figure();   hold on;    grid on;
%  plot(freq1,10*log10(psdx1));
%  title("Signal FM - Periodogram using FFT");
%  xlabel("Frequency (Hz)");
%  ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
%  hold off;

%% Modulation
load('asen3300mod.mat')
L = length(signal);
T = 1/fs;
time = 0:T:(L-1)*T;

%plot time series
plot(time(48000:96001),signal(48000:96001))
title("Time Series");
xlabel("Time (s)");
ylabel("Amplitude (V)");

%find fft
data_fft = fft(signal);   %run Fast Fourier Transform
psdx1 = (1/(1000*L))*abs(data_fft).^2;
psdx1 (2:end-1) = 2*psdx1(2:end-1);
freq1 = fs/L*(0:L-1);

%plot fft
 %plot power spectral density
%  figure();   hold on;    grid on;
%  plot(freq1,10*log10(psdx1));
%  title("Modulated Signal - Periodogram using FFT");
%  xlabel("Frequency (Hz)");
%  ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
%  hold off;

 signal_AM = demod(signal,fc,fs,'am');
 signal_FM = demod(signal,fc,fs,'fm');

 % figure();
 % plot(time,signal_AM);
 % title("Signal Demodulated using AM");

figure();
plot(time(48000:96001),signal_FM(48000:96001));
title("Signal Demodulated using FM");

%find fft
data_fft = fft(signal_FM);   %run Fast Fourier Transform
psdx1 = (1/(1000*L))*abs(data_fft).^2;
psdx1 (2:end-1) = 2*psdx1(2:end-1);
freq1 = fs/L*(0:L-1);

%plot fft
 %plot power spectral density
 figure();   hold on;    grid on;
 plot(freq1(1:(length(freq1)/2)),10*log10(psdx1(1:(length(psdx1)/2))));
 title("Demodulated Signal (FM) - Periodogram using FFT");
 xlabel("Frequency (Hz)");
 ylabel("Power/Frequency dB(Vrms^2/Hz");
%print('FreqPlot1_2','-dpng');
 hold off;

%% Noisy signal
%% Modulation
load('asen3300mod_noisy.mat')
L = length(signalnoisy);
T = 1/fs;
time = 0:T:(L-1)*T;

%plot time series
plot(time(48000:96001),signalnoisy(48000:96001))
title("Time Series");
xlabel("Time (s)");
ylabel("Amplitude (V)");

%find fft
data_fft = fft(signalnoisy);   %run Fast Fourier Transform
psdx1 = (1/(1000*L))*abs(data_fft).^2;
psdx1 (2:end-1) = 2*psdx1(2:end-1);
freq1 = fs/L*(0:L-1);

%plot fft
 %plot power spectral density
%  figure();   hold on;    grid on;
%  plot(freq1,10*log10(psdx1));
%  title("Modulated Signal - Periodogram using FFT");
%  xlabel("Frequency (Hz)");
%  ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
%  hold off;

 signal_AM = demod(signalnoisy,fc,fs,'am');
 signal_FM = demod(signalnoisy,fc,fs,'fm');

 % figure();
 % plot(time,signal_AM);
 % title("Signal Demodulated using AM");

figure();
plot(time(48000:96001),signal_FM(48000:96001));
title("Signal Demodulated using FM");

%find fft
data_fft = fft(signal_FM);   %run Fast Fourier Transform
psdx1 = (1/(1000*L))*abs(data_fft).^2;
psdx1 (2:end-1) = 2*psdx1(2:end-1);
freq1 = fs/L*(0:L-1);

%plot fft
 %plot power spectral density
 figure();   hold on;    grid on;
 plot(freq1(1:(length(freq1)/2)),10*log10(psdx1(1:(length(psdx1)/2))));
 title("Demodulated Signal (FM) - Periodogram using FFT");
 xlabel("Frequency (Hz)");
 ylabel("Power/Frequency dB(Vrms^2/Hz");
%print('FreqPlot1_2','-dpng');
 hold off;