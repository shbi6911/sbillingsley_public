%--------------------------------------------------------------------------
%
%     TEST CASE FOR RUNNING JB2008
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
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

year = 2001;
doy = 200;
[month,day,hour,minute,sec] = days2mdh(year,doy);
MJD_UTC = Mjday(year,month,day,hour,minute,sec);

% USE 1 DAY LAG FOR F10 AND S10 FOR JB2008
JD = floor(MJD_UTC-1+2400000.5);
i = find(JD==SOLdata(3,:),1,'first');
SOL = SOLdata(:,i);
F10 = SOL(4);
F10B = SOL(5);
S10 = SOL(6);
S10B = SOL(7);

% USE 2 DAY LAG FOR M10 FOR JB2008
SOL = SOLdata(:,i-1);
XM10 = SOL(8);
XM10B = SOL(9);

% USE 5 DAY LAG FOR Y10 FOR JB2008
SOL = SOLdata(:,i-4);
Y10 = SOL(10);
Y10B = SOL(11);

% doy = finddays(year,month,day,hour,minute,sec);
i = find(year==DTCdata(1,:) & floor(doy)==DTCdata(2,:),1,'first');
DTC = DTCdata(:,i);
ii = floor(hour)+3;
DSTDTC = DTC(ii);

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
SAT(3) = 150;

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
[TEMP,RHO] = JB2008(MJD_UTC,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC)

