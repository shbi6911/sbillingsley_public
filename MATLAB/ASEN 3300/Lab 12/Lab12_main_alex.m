clc; clear; close all;

oddIn = 'GPSLOG03_ODD_Static.txt';
oddOut = 'OddGPS';

evenIn = 'GPSLOG05_EVEN_Static.txt';
evenOut = 'EvenGPS';

dynamicIn = 'GPSLOG01_Dynamic.txt';
dynamicOut = 'DynamicGPS';

oddGPS = GPS_parser(oddIn,oddOut);
evenGPS = GPS_parser(evenIn,evenOut);
dynamicGPS = GPS_parser(dynamicIn,dynamicOut);

oddLLA = [oddGPS.lat' oddGPS.long' oddGPS.alt'];
evenLLA = [evenGPS.lat' evenGPS.long' evenGPS.alt'];
dynamicLLA = [dynamicGPS.lat' dynamicGPS.long' dynamicGPS.alt'];

oddXYZ = lla2ecef(oddLLA);
evenXYZ = lla2ecef(evenLLA);
dynamicXYZ = lla2ecef(dynamicLLA);

oddMn = mean(oddXYZ)
evenMn = mean(evenXYZ)
dynamicMn = mean(dynamicXYZ)

oddDev = std(oddXYZ)
evenDev = std(evenXYZ)
dynamicXYZ = std(dynamicXYZ)

figure(1)
hold on
grid on
scatter3(oddXYZ(:,1),oddXYZ(:,2),oddXYZ(:,3),[],oddXYZ(:,3))
colormap("turbo")
colorbar
xlabel('X-Coordinate [m]')
ylabel('Y-Coordinate [m]')
zlabel('Z_coordinate [m]')
title('Plot of Odd GPS Coordinates')

pos_acc_even = norm(evenDev);
pos_acc_odd = norm(oddDev);




