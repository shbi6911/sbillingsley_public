%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

%% Define constants
function const = getConst()
    const.mu.Earth = 3.986004415e5;
    const.mu.Moon = 4902.799;
    const.mu.Mars = 4.305e4;        %km^3/s^2
    const.mu.Sun = 1.32712428e11;
    const.mu.Saturn = 3.794e7;
    const.mu.Jupiter = 1.268e8;

    const.radius.Earth = 6378.1363;
    const.radius.Moon = 1738;
    const.radius.Mars = 3397.2;     %km
    const.radius.Saturn = 60268;
    const.radius.Jupiter = 71492;

    const.a.Earth = 1.0000010178;   %AU
    const.a.Mars = 1.52;
    const.a.Saturn = 9.554909595;
    const.a.Jupiter = 5.202603191;
    
    const.AU = 149597870.7;         %km
    const.G = 6.6743e-20;           %km^3/(kg*s^2)
end