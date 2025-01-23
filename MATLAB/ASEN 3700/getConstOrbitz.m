%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024

function const = getConstOrbitz()
%getConstOrbitz outputs a structure of planetary parameters relevant to
%calculation in ASEN3700, referenced from Orbital Mechanics by Curtis
%
%INPUTS     none 
% 
%OUTPUTS    const   a struct of constants for orbital parameters
%           includes Sun, Mercury, Venus, Earth, Luna, Mars, Jupiter,
%           Saturn, Uranus, Neptune, and Pluto
%           (fields)
%           .mu     gravitational parameter (km^3/s^2)
%           .SOI    radius of gravitational sphere of influence (km)
%           .radius average radius (km)
%           .mass   mass (kg)
%           .rotation   sidereal rotation period (s)
%           .inclination    inclination of the body's equator to its orbit
%           .a      semi-major axis of orbit (km) (N/A for Sun)
%           .e      eccentricity of orbit (unitless) (N/A for Sun)
%           .ecliptic   inclination of body's orbit to Solar System
%               ecliptic plane (radians)    (N/A for Sun)
%           .period orbital period around the Sun (s)   (N/A for Sun)
%METHODOLOGY    this function outputs a struct of constants from Table A-1
%and Table A-2 in Orbital Mechanics by Curtis
%
%   mu values in km^3/s^2
const.Sun.mu = 132712440018;
const.Mercury.mu = 22032;
const.Venus.mu = 324859;
const.Earth.mu = 398600;
const.Luna.mu = 4905;
const.Mars.mu = 42828;
const.Jupiter.mu = 126686534;
const.Saturn.mu = 37931187;
const.Uranus.mu = 5793939;
const.Neptune.mu = 6836529;
const.Pluto.mu = 871;
%   SOI radius in km
const.Sun.SOI = 112000;
const.Mercury.SOI = 616000;
const.Venus.SOI = 925000;
const.Earth.SOI = 925000;
const.Luna.SOI = 66100;
const.Mars.SOI = 577000;
const.Jupiter.SOI = 48200000;
const.Saturn.SOI = 54800000;
const.Uranus.SOI = 51800000;
const.Neptune.SOI = 86600000;
const.Pluto.SOI = 3080000;
%   average radii (km)
const.Sun.radius = 696000;
const.Mercury.radius = 2440;
const.Venus.radius = 6052;
const.Earth.radius = 6378;
const.Luna.radius = 1737;
const.Mars.radius = 3396;
const.Jupiter.radius = 71490;
const.Saturn.radius = 60270;
const.Uranus.radius = 25560;
const.Neptune.radius = 24764;
const.Pluto.radius = 1187;
% masses (kg)
const.Sun.mass = 1.989e30;
const.Mercury.mass = 330.2e21;
const.Venus.mass = 4.869e24;
const.Earth.mass = 5.974e24;
const.Luna.mass = 73.48e21;
const.Mars.mass = 641.9e21;
const.Jupiter.mass = 1.899e27;
const.Saturn.mass = 568.5e24;
const.Uranus.mass = 86.83e24;
const.Neptune.mass = 102.4e24;
const.Pluto.mass = 13.03e21;
%sidereal period (s)
const.Sun.rotation = 25.38*24*3600;
const.Mercury.rotation = 58.65*24*3600;
const.Venus.rotation = 243*24*3600;
const.Earth.rotation = 23.9345*3600;
const.Luna.rotation = 27.32*24*3600;
const.Mars.rotation = 24.62*3600;
const.Jupiter.rotation = 9.925*3600;
const.Saturn.rotation = 10.66*3600;
const.Uranus.rotation = 17.24*3600;
const.Neptune.rotation = 16.11*3600;
const.Pluto.rotation = 6.387*24*3600;
%inclination of equator to planetary orbit plane (radians)
const.Sun.inclination = deg2rad(7.25);
const.Mercury.inclination = deg2rad(0.01);
const.Venus.inclination = deg2rad(177.4);
const.Earth.inclination = deg2rad(23.45);
const.Luna.inclination = deg2rad(6.68);
const.Mars.inclination = deg2rad(25.19);
const.Jupiter.inclination = deg2rad(3.13);
const.Saturn.inclination = deg2rad(26.73);
const.Uranus.inclination = deg2rad(97.77);
const.Neptune.inclination = deg2rad(28.32);
const.Pluto.inclination = deg2rad(122.5);
%semimajor axis of orbit (km) 
const.Mercury.a = 57.91e6;
const.Venus.a = 108.2e6;
const.Earth.a = 149.6e6;
const.Luna.a = 384.4e3;
const.Mars.a = 227.9e6;
const.Jupiter.a = 778.6e6;
const.Saturn.a = 1.433e9;
const.Uranus.a = 2.872e9;
const.Neptune.a = 4.495e9;
const.Pluto.a = 5.906e9;
%eccentricity
const.Mercury.e = 0.2056;
const.Venus.e = 0.0067;
const.Earth.e = 0.0167;
const.Luna.e = 0.0549;
const.Mars.e = 0.0935;
const.Jupiter.e = 0.0489;
const.Saturn.e = 0.0565;
const.Uranus.e = 0.0457;
const.Neptune.e = 0.0113;
const.Pluto.e = 0.2488;
%inclination to ecliptic plane (radians)
const.Mercury.ecliptic = deg2rad(7.00);
const.Venus.ecliptic = deg2rad(3.39);
const.Earth.ecliptic = deg2rad(0.00);
const.Luna.ecliptic = deg2rad(5.145);
const.Mars.ecliptic = deg2rad(1.850);
const.Jupiter.ecliptic = deg2rad(1.304);
const.Saturn.ecliptic = deg2rad(2.485);
const.Uranus.ecliptic = deg2rad(0.772);
const.Neptune.ecliptic = deg2rad(1.769);
const.Pluto.ecliptic = deg2rad(17.16);
%orbital period around sun (s)
const.Mercury.period = 87.97*24*3600;
const.Venus.period = 224.7*24*3600;
const.Earth.period = 365.256*24*3600;
const.Luna.period = 27.322*24*3600;
const.Mars.period = 1.881*365*24*3600;
const.Jupiter.period = 11.86*365*24*3600;
const.Saturn.period = 29.46*365*24*3600;
const.Uranus.period = 84.01*365*24*3600;
const.Neptune.period = 164.8*365*24*3600;
const.Pluto.period = 247.9*365*24*3600;
end