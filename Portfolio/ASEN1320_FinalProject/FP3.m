%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

%% loading and initializing data
load Sat1Position.csv;
load Sat2Position.csv;
load Sat1Visibility.csv;
load Sat2Visibility.csv;
load CBPosition.csv;
Sat1Visibility = logical(Sat1Visibility);
Sat2Visibility = logical(Sat2Visibility);
Sat1Comm = Sat1Position(Sat1Visibility,:);
Sat2Comm = Sat2Position(Sat2Visibility,:);
Sat1GrndComm = CBPosition(Sat1Visibility,:);
Sat2GrndComm = CBPosition(Sat2Visibility,:);

%% create plot

[X,Y,Z] = sphere(100);
figure;
hold on; grid on;
axis([-7000 7000 -7000 7000 -7000 7000],'manual','equal');
title('Plot of ISS and Hubble orbits over one day, showing Communication Points with Canberra Ground Station');
xlabel('x-distance from Earth Center (km)');
ylabel('y-distance from Earth Center (km)');
zlabel('z-distance from Earth Center (km)');
surf((X*6378),(Y*6378),(Z*6378),'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.35);
plot3(0,0,0,'r.','MarkerSize',30);

%% plot static image
%  
% view(-45,20);
% plot3(CBPosition(:,1), CBPosition(:,2),CBPosition(:,3),'g','LineWidth',2);
% plot3(Sat1GrndComm(:,1), Sat1GrndComm(:,2), Sat1GrndComm(:,3), 'k.', 'MarkerSize',20);
% plot3(Sat2GrndComm(:,1), Sat2GrndComm(:,2), Sat2GrndComm(:,3), 'r.', 'MarkerSize',20);
% plot3(Sat1Position(:,1),Sat1Position(:,2),Sat1Position(:,3),'r','LineWidth',2);
% plot3(Sat1Comm(:,1),Sat1Comm(:,2),Sat1Comm(:,3),'k.','MarkerSize',20);
% plot3(Sat2Position(:,1),Sat2Position(:,2),Sat2Position(:,3),'k','LineWidth',2);
% plot3(Sat2Comm(:,1),Sat2Comm(:,2),Sat2Comm(:,3),'r.','MarkerSize',20);
% legend ('Earth Surface','Earth Center','Canberra','','','ISS','ISS Comm Points','Hubble','Hubble Comm Points');

%% create animated lines and initial objects

view(-85,20);
Sat1An = animatedline(Sat1Position(1,1), Sat1Position(1,2),Sat1Position(1,3),...
    'LineWidth',2, 'Color','r','MaximumNumPoints',45);
Sat2An = animatedline(Sat2Position(1,1), Sat2Position(1,2),Sat2Position(1,3),...
    'LineWidth',2,'Color','k','MaximumNumPoints',45);
GrndAn = animatedline(CBPosition(1,1), CBPosition(1,2), CBPosition(1,3),...
    'LineWidth',2,'Color','y');

%% update the animated lines in a loop

for ii = 2:1441
    addpoints(Sat1An,Sat1Position(ii,1),Sat1Position(ii,2),Sat1Position(ii,3));
    addpoints(Sat2An,Sat2Position(ii,1),Sat2Position(ii,2),Sat2Position(ii,3));
    addpoints(GrndAn,CBPosition(ii,1),CBPosition(ii,2),CBPosition(ii,3));
    if Sat1Visibility(ii)
        plot3(Sat1Position(ii,1),Sat1Position(ii,2),Sat1Position(ii,3),'y.','MarkerSize',20);
        plot3(CBPosition(ii,1), CBPosition(ii,2), CBPosition(ii,3), 'r.', 'MarkerSize',20);
    end
    if Sat2Visibility(ii)
        plot3(Sat2Position(ii,1),Sat2Position(ii,2),Sat2Position(ii,3),'y.','MarkerSize',20);
        plot3(CBPosition(ii,1), CBPosition(ii,2), CBPosition(ii,3), 'k.', 'MarkerSize',20);
    end
    view((-115+(ii/4)), 20);
    %saveas(gcf,sprintf('image%d.jpg',ii));
    drawnow
end

