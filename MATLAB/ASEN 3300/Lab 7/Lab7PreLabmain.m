%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

%%Question 1
% bits = 4;
% min_v = 0;
% max_v = 3.3;
% voltage_range = 0:0.25:3.25;
% bin_range = Voltage2Bin(min_v,max_v,bits,voltage_range);
% binary_bin_range_01 = dec2bin(bin_range);
% 
% figure();   hold on;    grid on;
% scatter(voltage_range,bin_range);
% xlabel("Voltage (V)");
% ylabel("Bin Number");
% title("Voltage vs. Bin Number (4-bit)");
% hold off;
% 
% bits = 8;
% bin_range = Voltage2Bin(min_v,max_v,bits,voltage_range);
% binary_bin_range_02 = dec2bin(bin_range);
% figure();   hold on;    grid on;
% scatter(voltage_range,bin_range);
% xlabel("Voltage (V)");
% ylabel("Bin Number");
% title("Voltage vs. Bin Number (8-bit)");
% hold off;
% 
% bits = 12;
% bin_range = Voltage2Bin(min_v,max_v,bits,voltage_range);
% binary_bin_range_03 = dec2bin(bin_range);
% figure();   hold on;    grid on;
% scatter(voltage_range,bin_range);
% xlabel("Voltage (V)");
% ylabel("Bin Number");
% title("Voltage vs. Bin Number (12-bit)");
% hold off;

%%Question 2
amplitude = 3.3;    %min voltage is 0, max is 3.3, b/c (3.3)/2 = 1.65
offset = 1.65;
min_v = (amplitude/2)-offset;
max_v = (amplitude/2)+offset;
bits = 12;
voltage_range = (amplitude/2)*sin(0:(pi/30):(2*pi))+offset;
bin_range = Voltage2Bin(min_v,max_v,bits,voltage_range);
figure();   hold on; grid on;
plot(bin_range);
title("Sinusoidal Bin Number vs. Array Number");
xlabel("Index");
ylabel("Bin Number");
hold off;

function bin = Voltage2Bin(min_v, max_v, bits, voltage)
%INPUTS     min_v   scalar value of minimum voltage range
%           max_v   scalar value of maximum voltage range
%           bits    scalar value of number of bits in range
%           voltage scalar or vector value of voltage to convert to bins
%
%OUTPUTS    vector of bin numbers that voltage fits into
%
%METHODOLOGY    This function takes in a voltage range and number of bits,
%divides that range into bins, and finds which bin the input voltage fits
%into, and outputs that number.
num_bins = 2^bits;
slope = (max_v - min_v)/num_bins;
bin = floor(voltage./slope);
end
