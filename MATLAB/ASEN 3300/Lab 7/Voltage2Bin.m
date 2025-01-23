%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

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

