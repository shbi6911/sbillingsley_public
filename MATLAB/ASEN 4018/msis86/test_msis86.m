clc
clear
format long g

ap=zeros(7,1);
an=zeros(6,1);
sa=zeros(3,1);
su=zeros(2,1);
sf=zeros(3,1);
te=zeros(2,1);

rl = [-60, -30, 0, 30, 60];

iday  = 172;	% day of the year
sec   = 29000.;	% utc in seconds
alt   = 550.;	% altitude in km
glat  = 60.;	% latitude geodetic in degrees
glong = -70.;	% longitude in degrees
stl	  = 16.;	% local solar time in hours
f107a = 150.;	% mean solar flux over 90 days
f107  = 150.;	% solar flux from previous day
ap(1) = 4;		% dayly Ap
ap(2) = 4;		% 3 hr ap index for current time  
ap(3) = 4;		% 3 hr ap index for 3 hrs before current time
ap(4) = 4;		% 3 hr ap index for 6 hrs before current time
ap(5) = 4;		% 3 hr ap index for 9 hrs before current time
ap(6) = 4;		% average of eight 3 hr ap indicies from 12 to 33 hrs prior to current time
ap(7) = 4;		% average of eight 3 hr ap indicies from 36 to 59 hrs prior to current time

[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

iday  = 81;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

iday  = 172;
sec   = 75000.;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

sec   = 29000.;
alt   = 100.;
ap(7) = 40;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

alt   = 400.;
glat  = 0.;
ap(6) = 40;
ap(7) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

glat  = 60.;
glong = 0.;
ap(5) = 40;
ap(6) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

glong = -70.;
stl	  = 4.;
ap(4) = 40;
ap(5) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

stl	  = 16.;
f107a = 70.;
ap(3) = 40;
ap(4) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

f107a = 150.;
f107  = 180.;
ap(2) = 40;
ap(3) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

f107  = 150.;
ap(1) = 40;
ap(2) = 0;
[d, t] = msis86(iday, sec, alt, glat, glong, stl, f107a, f107, ap);
fprintf (' %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %10.4g %10.4g \n',...
	d(1), d(2), d(3), d(4), d(5), d(6), d(7), d(8), t(1), t(2));

