x_vec = -360:1:360;
y_vec = -180:1:180;
[X,Y] = meshgrid(x_vec,y_vec);

Z = sind(X) + cos(Y).^2;

figure(1);
surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none');