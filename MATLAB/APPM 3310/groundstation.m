%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

initialpos = [-4460.49;2682.08;-3674.26]; %[km]
tspan= 0:86400; %seconds in a day
rotangle = tspan*(360/86400); %rotation angle in degrees per second
rotanglerad = rotangle.*(pi/180); %rotation angle in radians per second
axis = [0;0;1]; %rotate around the z-axis

rotations = zeros(3,length(rotanglerad));
for i=1:length(rotanglerad)
    newvec = rotate(rotanglerad(i),axis,initialpos);
    rotations(:,i)=newvec;
end

[X,Y,Z] = sphere(100);
figure (1); hold on; grid on;
plot3(-4460.49,2682.08,-3674.26,'.','MarkerSize',30,'Color',[0.9290 0.6940 0.1250]);
plot3(rotations(1,:),rotations(2,:),rotations(3,:),'r','LineWidth',2);
surf((X*6378),(Y*6378),(Z*6378),'FaceColor',[0 0.4470 0.7410],...
    'EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.35);
plot3(0,0,0,'r.','MarkerSize',30);
title ("Rotation of a Ground Station Over One Day");
xlabel ("X-Position [km]");
ylabel ("Y-Position [km]");
zlabel ("Z-Position [km]");
legend ("Canberra Ground Station Initial Location","Ground Station Path");