%--------------------------------------------------------------------------
%
%     TEST CASE FOR RUNNING JB2006
% 
% Last modified:   2022/12/14   Meysam Mahooti
%
%--------------------------------------------------------------------------
clc
clear all
format long g

global const PC % astronomical constants & planetary coefficients
SAT_Const
constants
load DE440Coeff.mat
PC = DE440Coeff;

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
for i=1:33
    tline = fgetl(fid);
end
numrecsobs = str2num(tline(21:end));
tline = fgetl(fid);
for i=1:numrecsobs
    eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
end
for i=1:4
    tline = fgetl(fid);
end
numrecspred = str2num(tline(22:end));
tline = fgetl(fid);
for i=numrecsobs+1:numrecsobs+numrecspred
    eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
end
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ---------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c Ssrc
%  ---------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read Ap data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d',[12 inf]);
fclose(fid);

year = 2001;
doy = 200;
[month,day,hour,minute,sec] = days2mdh(year,doy);
MJD_UTC = Mjday(year,month,day,hour,minute,sec);

% SET SOLAR INDICES
JD = floor(MJD_UTC+2400000.5);
i = find(JD==SOLdata(3,:),1,'first');
SOL = SOLdata(:,i);
F10B = SOL(5);
S10B = SOL(7);
XM10 = SOL(8);

% USE 1 DAY LAG FOR EUV
SOL = SOLdata(:,i-1);
F10 = SOL(4);
S10 = SOL(6);

% USE 5 DAY LAG FOR MG FUV INFLUENCE
SOL = SOLdata(:,i-5);
XM10B = SOL(9);

% USE 6.7 HR LAG FOR Ap INFLUENCE
[year,month,day,hour,minute,sec] = invjday(MJD_UTC-6.7/24);
doy = finddays(year,month,day,hour,minute,sec);
i = find(year==APdata(1,:) & floor(doy)==APdata(2,:),1,'first');
AP = APdata(5:12,i);

if hour < 3
    Ap = AP(1);
elseif hour < 6
    Ap = AP(2);
elseif hour < 9
    Ap = AP(3);
elseif hour < 12
    Ap = AP(4);
elseif hour < 15
    Ap = AP(5);
elseif hour < 18
    Ap = AP(6);
elseif hour < 21
    Ap = AP(7);
else
    Ap = AP(8);
end

GEO(1) = F10;
GEO(2) = F10B;
GEO(3) = Ap;

% CONVERT POINT OF INTEREST LOCATION (RADIANS AND KM)
% CONVERT LONGITUDE TO RA
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT = MJD_UTC + TT_UTC/86400;
GWRAS = iauGmst06(const.DJM0, MJD_UT1, const.DJM0, MJD_TT);
XLON = 60*const.Rad;
SAT(1) = mod(GWRAS + XLON, const.pi2);
SAT(2) = -70*const.Rad; 
SAT(3) = 400;

% Difference between ephemeris time and universal time
% [year, month, day, hour, minute, sec] = invjday(MJD_UTC);
% days = finddays(year, month, day, hour, minute, sec);
% ET_UT = ETminUT(year+days/365.25);
% MJD_ET = MJD_UTC+ET_UT/86400;
% [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,...
%  r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_ET);

MJD_TDB = Mjday_TDB(MJD_TT);
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB);

% SET Sun's right ascension and declination (RADIANS)
ra_Sun  = atan2(r_Sun(2), r_Sun(1));
dec_Sun = atan2(r_Sun(3), sqrt(r_Sun(1)^2+r_Sun(2)^2));
SUN(1)  = ra_Sun;
SUN(2)  = dec_Sun;

% COMPUTE DENSITY KG/M3 RHO
[TEMP,RHO] = JB2006(MJD_UTC,SUN,SAT,GEO,S10,S10B,XM10,XM10B)

