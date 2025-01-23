%By:        Shane Billingsley
%Class:     ASEN 3300 Aerospace Electronics & Communications
%Date:      Spring 2024

clear; close all;

fopen("LAB12_EVEN_STATIC.txt");
gps_even = GPS_parser("GPSLOG05_EVEN_Static.TXT","LAB12_EVEN_STATIC.txt");
fopen("LAB12_ODD_STATIC.txt");
gps_odd = GPS_parser("GPSLOG03_ODD_Static.TXT","LAB12_ODD_STATIC.txt");
fopen("LAB12_DYNAMIC.txt");
gps_dynamic = GPS_parser("GPSLOG01_Dynamic.TXT","LAB12_DYNAMIC.txt");

gps_even_lat_mean = mean(gps_even.lat);
gps_even_lat_stddev = std(gps_even.lat);

gps_even_long_mean = mean(gps_even.long);
gps_even_long_stddev = std(gps_even.long);

gps_even_alt_mean = mean(gps_even.alt);
gps_even_alt_stddev = std(gps_even.alt);

gps_odd_lat_mean = mean(gps_odd.lat);
gps_odd_lat_stddev = std(gps_odd.lat);

gps_odd_long_mean = mean(gps_odd.long);
gps_odd_long_stddev = std(gps_odd.long);

gps_odd_alt_mean = mean(gps_odd.alt);
gps_odd_alt_stddev = std(gps_odd.alt);
