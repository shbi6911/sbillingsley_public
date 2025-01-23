%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

%This script generates a vector of temperature steps, calculates speed of
%sound for those temperatures, generates a random speed and selects a
%random speed of sound, then calculates a Mach number based on these chosen
%values.  It then outputs flow regime and Mach number.

TvectorA = -70:10:100;      %creating temperature vectors and concatenating
TvectorB = 150:50:1500;
Tvector = [TvectorA TvectorB];
SoundSpeedvector = sqrt (1.4 * 287.058 * (Tvector + 273.16));   %calculating speed of sound in air at each temperature
rng(uint64(now*1000));      %seed random number generator
x=randi(46);                %generate random integer and use it to choose random speed of sound value
SoundSpeed = SoundSpeedvector(x);
VehicleSpeed = 10 + (1500-10)*rand;     %generate a random vehicle speed between 10 and 1500
MachNumber = VehicleSpeed / SoundSpeed; %assign a Mach number to the speed
Regimes = ["Incompressible", "Subsonic", "Transonic", "Sonic", "Supersonic", "Hypersonic"];
if MachNumber<0.3                       %determine flow regime based on Mach number
    y = Regimes (1);
elseif MachNumber < 0.8
    y = Regimes (2);
elseif MachNumber < 1
    y = Regimes (3);
elseif MachNumber == 1
    y = Regimes (4);
elseif MachNumber < 5
    y = Regimes (5);
elseif MachNumber > 5
    y = Regimes (6);
end
fprintf ('%s flight regime and MachNumber is %7.6f\n', y, MachNumber);  %output flow regime and Mach number
