%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

clear; close all;

%%10 Hz signal

% data = readmatrix("teraterm.log");
% data = data(2:(end-1));
% 
% N=length(data);
% data_fft = fft(data);   %run Fast Fourier Transform
% 
% data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs
% 
% psdx1 = (1/(1000*N))*abs(data_fft).^2;
% psdx1 (2:end-1) = 2*psdx1(2:end-1);
% 
% freq1 = linspace(0,1000/2,length(psdx1));
% signal = find(max(psdx1));
% 
% %plot power spectral density
% figure();   hold on;    grid on;
% plot(freq1,10*log10(psdx1));
% title("Signal 1 - Periodogram using FFT - Entire Vector");
% xlabel("Frequency (Hz)");
% ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
% hold off;

%%100 Hz signal
% data = readmatrix("teraterm100Hz");
% data = data(2:(end-1));
% 
% N=length(data);
% data_fft = fft(data);   %run Fast Fourier Transform
% 
% data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs
% 
% psdx1 = (1/(1000*N))*abs(data_fft).^2;
% psdx1 (2:end-1) = 2*psdx1(2:end-1);
% 
% freq1 = linspace(0,1000/2,length(psdx1));
% signal = find(max(psdx1));
% 
% %plot power spectral density
% figure();   hold on;    grid on;
% plot(freq1,10*log10(psdx1));
% title("Signal 2 - Periodogram using FFT - Entire Vector");
% xlabel("Frequency (Hz)");
% ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
% hold off;

%%600 Hz signal
% data = readmatrix("teraterm600Hz");
% data = data(2:(end-1));
% 
% N=length(data);
% data_fft = fft(data);   %run Fast Fourier Transform
% 
% data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs
% 
% psdx1 = (1/(1000*N))*abs(data_fft).^2;
% psdx1 (2:end-1) = 2*psdx1(2:end-1);
% 
% freq1 = linspace(0,1000/2,length(psdx1));
% signal = find(max(psdx1));
% 
% %plot power spectral density
% figure();   hold on;    grid on;
% plot(freq1,10*log10(psdx1));
% title("Signal 600 Hz - Periodogram using FFT - Entire Vector");
% xlabel("Frequency (Hz)");
% ylabel("Power/Frequency dB(Vrms^2/Hz");
% %print('FreqPlot1_2','-dpng');
% hold off;

%%1100 Hz signal
data = readmatrix("teraterm1100Hz");
data = data(2:(end-1));

N=length(data);
data_fft = fft(data);   %run Fast Fourier Transform

data_fft = data_fft(1:floor(N/2)+1);    %discard negative freqs

psdx1 = (1/(1000*N))*abs(data_fft).^2;
psdx1 (2:end-1) = 2*psdx1(2:end-1);

freq1 = linspace(0,1000/2,length(psdx1));
signal = find(max(psdx1));

%plot power spectral density
figure();   hold on;    grid on;
plot(freq1,10*log10(psdx1));
title("Signal 1100 Hz - Periodogram using FFT - Entire Vector");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency dB(Vrms^2/Hz");
%print('FreqPlot1_2','-dpng');
hold off;