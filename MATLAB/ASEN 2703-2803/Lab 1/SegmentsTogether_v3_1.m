%By:        Shane Billingsley
%Class:     ASEN 2803 Dynamics & Controls Lab
%Date:      Spring 2023

clc; close all; clear;

%% Creating 3D Plot

m = 200;
g = 9.81;
xDistance = 0;

startPoint = [0 0 125];
splitDrop2Trans = [60 0 45];
splitTransFlat = [100 0 25];
splitTrans2Parabola = [140 0 45];
splitParabola2Level = [293.6 0 45];
splitLevel2Loop = [333.6 0 25];
splitLoop2Bank = [413.6 0 25];
splitBankPeak = [453.6 40 25];
splitBank2Brake = [413.6 80 25];
endPoint = [52 80 25];

% 3D Position Plot
f1 = figure(1);


%% Points of interest
figure(1),plot3(startPoint(1), startPoint(2), startPoint(3), '|', MarkerSize=10); hold on
% figure(1),plot3(split0_6(1), split0_6(2), split0_6(3), '.', MarkerSize=15); hold on
figure(1),plot3(splitTransFlat(1), splitTransFlat(2), splitTransFlat(3), '|', MarkerSize=10); hold on
% figure(1),plot3(split1_3(1), split1_3(2), split1_3(3), '.', MarkerSize=15); hold on
% figure(1),plot3(split1_6(1), split1_6(2), split1_6(3), '.', MarkerSize=15); hold on
figure(1),plot3(splitLevel2Loop(1), splitLevel2Loop(2), splitLevel2Loop(3), '|', MarkerSize=10); hold on
figure(1),plot3(splitLoop2Bank(1), splitLoop2Bank(2), splitLoop2Bank(3), '|', MarkerSize=10); hold on
% figure(1),plot3(split3_5(1), split3_5(2), split3_5(3), '.', MarkerSize=15); hold on
figure(1),plot3(splitBank2Brake(1), splitBank2Brake(2), splitBank2Brake(3), '|', MarkerSize=10); hold on
figure(1),plot3(endPoint(1), endPoint(2), endPoint(3), '|', MarkerSize=10); hold on

%% Section 1 - Drop

% Put in g force calcutions

% Straight line
xDrop = 0:60;
yDrop = zeros([1 length(xDrop)]);
zDrop = -4/3*xDrop + 125;
figure(1),plot3(xDrop, yDrop, zDrop);

% Curved Part
xDrop2Parabola = 60:140;
yDrop2Parabola = zeros([1 length(xDrop2Parabola)]);
zDrop2Parabola = -sqrt(2500-(xDrop2Parabola-100).^2) + 75;
figure(1),plot3(xDrop2Parabola, yDrop2Parabola, zDrop2Parabola);

% G-Forces

% velocity at the end of the straight drop portion
v_anytime = sqrt(2*g*(startPoint(3)-splitDrop2Trans(3)));
% delta_x_Drop = sqrt((startPoint(3)-split0_6(3))^2 + (startPoint(1)-split0_6(1))^2);
% a_Drop = (v_flat^2)/(2*delta_x_Drop);

gForcesDropUp = ones([1 length(xDrop)]) * cos(atan(4/3));
gForces_Drop2ParabolaUp = cos(atan((xDrop2Parabola-100)./sqrt(-xDrop2Parabola.^2+200*xDrop2Parabola-7500))) + (50+sqrt(2500-(xDrop2Parabola-100).^2))/25;


%% Section 2 - Parabola

% Parabola
xParabola = 140:293.6;
yParabola = zeros([1 length(xParabola)]);
zParabola = 45 + (xParabola-140)*(4/3) - (9.81*(xParabola-140).^2)/(3139.2*(cos(atan(4/3)))^2);
figure(1),plot3(xParabola, yParabola, zParabola);

% Curved Part
xParabola2Level = 293.6:0.5:333.6;
yParabola2Level = zeros([1 length(xParabola2Level)]);
zParabola2Level = -sqrt(2500-(xParabola2Level-333.6).^2) + 75;
figure(1),plot3(xParabola2Level, yParabola2Level, zParabola2Level);

h_start = 125; %starting height of coaster(m)
h_incoming = 45; %starting height of parabola (m)
x_0 = 0; %starting horizontal position for parabola (m)

v_incoming = sqrt(2*g*(h_start-h_incoming)); %incoming velocity (m/s)

theta_para = 4/3; %incoming slope
angle_launch = (theta_para); %incoming angle (degrees)

v_x0_para = v_incoming*cos(angle_launch); %initial x component velocity (m/s)
v_y0_para = v_incoming*sin(angle_launch); %initial y component velocity (m/s)

%calculate max height using conservation of energy
y_max_para = (v_y0_para^2)/(2*g) + h_incoming; %max height of parabola (m)

%calculate time using quadratic formula
t1 = (v_y0_para+sqrt(v_y0_para^2+2*g*(y_max_para-h_incoming)))/g;
t2 = (v_y0_para-sqrt(v_y0_para^2+2*g*(y_max_para-h_incoming)))/g;

t_para = 0:0.01:t1*2; %time interval (s)

x_para = x_0 + v_x0_para.*t_para; %x position as a function of time

y_para = h_incoming+tan(angle_launch).*(x_para-x_0)-0.5*g.*((x_para-x_0)./v_x0_para).^2; %y position paramaterized with x position

%clear values outside of parabolic segement
for i = 1:length(y_para)
    if y_para(i) < h_incoming
        y_para(i) = 0;
        t_para(i) = 0;
    end
end

t_para = [0,nonzeros(t_para)'];
y_para = nonzeros(y_para)';

%find x position corresponding to max height
for i = 1:length(y_para)
    if y_para(i) == max(y_para)
        x_mid = x_para(i);
    end
end

x_para = x_para(1,1:length(y_para));
x_end = x_para(length(x_para));

%calculate velocity
v_y_para = v_y0_para - g.*t_para;
v_x_para = ones(1,length(v_y_para)).*v_x0_para;
v_para = sqrt(v_y_para.^2+v_x_para.^2);

%derivatives for G force calculation
dydx = tan(angle_launch)-g.*((x_para-x_0)./v_x0_para^2);
d2ydx2 = -g/v_x0_para^2;
rho = ((1+dydx.^2).^(3/2))./abs(d2ydx2);

%calculate inner angle
phi = atan(dydx);

%G force
Gs_para = cos(phi)-((v_para.^2)./(g.*rho));

%g level out

gForces_Para2LevelUp = gForces_Drop2ParabolaUp(1:41);


%% Section 3 - Loop

% radius: 40m

thetaLoop = 0:pi/50:2*pi;
x_thetaLoop = 40*cos(thetaLoop) + (40 + 333.6);
y_thetaLoop = zeros([1 length(thetaLoop)]);
z_thetaLoop = 40*sin(thetaLoop) + (40+25);
figure(1),plot3(x_thetaLoop, y_thetaLoop, z_thetaLoop);

figure(1),plot3([splitLevel2Loop(1) splitLoop2Bank(1)], [splitLevel2Loop(2) splitLoop2Bank(2)], [splitLevel2Loop(3) splitLoop2Bank(3)]);

%G-forces
%create height vector
Vi = sqrt(9.81*2*(125-25));
Radius = 40;
hmax = 2*Radius;
hvec = 0:0.0001:hmax;

AbsVelvec = zeros(2*length(hvec), 1);
GLoopVec = zeros(2*length(hvec), 3);
dispvec = zeros(2*length(hvec), 1);

for i = 1:length(hvec)

AbsVelvec(i) = Velabs(Vi,g,hvec(i));
GLoopVec(i,:) = Gnorm(AbsVelvec(i),Radius);
GLoopVec(i,:) = GLoopVec(i,:) + Ggrav(g,hvec(i),Radius);
GLoopVec(i,:) = GLoopVec(i,:)./g;

GLoopVec(length(GLoopVec)+(1-i),:) = [-GLoopVec(i,1), GLoopVec(i,2), GLoopVec(i,3)];
AbsVelvec(length(GLoopVec)+(1-i)) = AbsVelvec(i);
dispvec(i) = hvec(i);
dispvec(length(AbsVelvec)+(1-i)) = dispvec(i);

end

GLoopVec(:,1)=0;
ArcLengthVec = linspace(0,Arcl(Radius),length(dispvec));


%% Section 4 - Bank

% radius: 40m

thetaBank = -pi/2:pi/50:pi/2;
z_thetaBank = zeros([1 length(thetaBank)]) + 25;
x_thetaBank = 40*cos(thetaBank) + 413.6;
y_thetaBank = 40*sin(thetaBank) + 40;
figure(1),plot3(x_thetaBank, y_thetaBank, z_thetaBank);

%G-forces

% define banked turn
theta_deg = 60;
theta_bnk = theta_deg*(pi/180);
rho = 40;
phi = linspace(0,pi,1000);
X = rho*sin(phi); Y = rho*cos(phi);
Z = zeros(1,1000);

% define other constants
g = 9.81;
h_0 = 125;
h = 25;
v_para = sqrt(2*g*(h_0-h));

%find g-forces
centri = v_para^2/(g*rho);
G_Bank = abs([0 (sin(theta_bnk)-(centri*cos(theta_bnk))) (cos(theta_bnk)+(centri*sin(theta_bnk)))]);

%% Section 5 - Braking

xBrake = 413.6:-1:52;
yBrake = ones([1 length(xBrake)]) * 80;
zBrake = ones([1 length(xBrake)]) * 25;
figure(1),plot3(xBrake, yBrake, zBrake);

xlim([0 500]);
ylim([0 500]);
zlim([0 500]);
title('Roller Coaster');

% G-Forces

% Incoming velocity I picked arbitrarily, will depend on previous values
v_intoBrake = velocity(zBrake(1));
v_Drop_Bottom = 0;
delta_x5 = splitBank2Brake(1) - endPoint(1);
a_Drop = (v_Drop_Bottom^2-v_intoBrake^2)/(2*delta_x5);

% 600 picked arbitrarily, used as place holder while other g-forces calced
gBrakeBack = ones([1 length(xBrake)]) * abs(a_Drop/g);

xDistance = 600 + flip(xBrake);

brakingForce = m * a_Drop;
fprintf('Section 5 (constant) G Forces: %s\n', num2str(gBrakeBack(1)))

%% Arc length
%get array of xyz points
position(1,1:61) = 0:60;
position(2,1:61) = 0;
position(3,1:61) = -4/3*(0:60) + 125;

position(1,62:141) = 61:140;
position(2,62:141) = 0;
position(3,62:141) = -sqrt(2500-((61:140)-100).^2) + 75;

position(1,142:294) = 141:293;
position(2,142:294) = 0;
position(3,142:294) = 45 + ((141:293)-140)*(4/3) - (9.81*((141:293)-140).^2)/(3139.2*(cos(atan(4/3)))^2);

position(1,295:334) = 294:333;
position(2,295:334) = 0;
position(3,295:334) = -sqrt(2500-((294:333)-333.6).^2) + 75;

position(1,335:374) = 334:373;
position(2,335:374) = 0;
position(3,335:374) = 25;

thetaArc = -pi/2 + pi/50 :pi/50: 3*pi/2 ;
position(1,375:474) = 40*cos(thetaArc) + (40 + 333.6);
position(2,375:474) = 0;
position(3,375:474) = 40*sin(thetaArc) + (40+25);

position(1,475:514) = 374:413;
position(2,475:514) = 0;
position(3,475:514) = 25;

thetaArc = -pi/2:pi/50:pi/2;
position(1,515:565) = 40*cos(thetaArc) + 413.6;
position(2,515:565) = 40*sin(thetaArc) + 40;
position(3,515:565) = 25;

position(1,566:927) = 413:-1:52;
position(2,566:927) = 80;
position(3,566:927) = 25;

% figure(3)
% plot3(position(1,:),position(2,:),position(3,:))
% xlim([0 500]);
% ylim([0 500]);
% zlim([0 500]);

trackLength(1) = 0;
for i = 2:927
    trackLength(i) = trackLength(i-1) + sqrt((position(1,i)-position(1,i-1))^2 + (position(2,i)-position(2,i-1))^2 + (position(3,i)-position(3,i-1))^2);
end

%% Plot G's

G_upDown = [gForcesDropUp, gForces_Drop2ParabolaUp(2:81), zeros(1,152), gForces_Para2LevelUp, ones(1,40), GLoopVec(1:16000:end,2).', ones(1,40), ones(1,50).*G_Bank(3), ones(1,362)];
figure(2)
subplot(2,2,1)
plot(trackLength,G_upDown)
title('Vertical Gs vs. Track Length');

G_side = [zeros(1,514), ones(1,51).*G_Bank(2), zeros(1,362)];
subplot(2,2,2)
plot(trackLength,G_side)
title('Horizontal Gs vs. Track Length');

G_frontBack = [zeros(1,565), gBrakeBack];
subplot(2,2,3)
hold on
plot(trackLength,G_frontBack)
plot([1250 1250],[0,0.2765])
title('Front/Back Gs vs. Track Length');
hold off

G_total_mag = zeros(1,927);
for i = 1:927
    G_total_mag(i) = sqrt(G_upDown(i).^2 + G_side(i).^2 + G_frontBack(i).^2);
end
subplot(2,2,4)
plot(trackLength,G_total_mag)
title('Total Magnitude of Gs vs. Track Length');

figure(5);
hold on; grid on;
c = colorbar;
c.Label.String = "Total Magnitude of Gs";
color_line3(position(1,:),position(2,:),position(3,:),G_total_mag,'LineWidth',3);
title ("Total Magnitude of Gs on Track");
xlabel ("X-Position (m)");
ylabel ("Y-Position (m)");
zlabel ("Z-Position (m)");
ylim([-10 90]);
xlim([-25 500]);
zlim([0 130]);
view(45,25);

hold off;

%% Functions

function v_anytime = velocity(height)
v_anytime = sqrt(2*9.81*(125-height));
end

function [Vabs] = Velabs(Vi,g,h)

    Vabs = sqrt((Vi^2)-(2*g*h));
end

function [Arcl] = Arcl(r)

    Arcl= 2*((pi*r));
end 

function [Gnorm] = Gnorm(Vabs,r)

    Gacc = (Vabs^2)/r;

    Gnorm = zeros(1,3);
    Gnorm(2) = Gacc;
end 

function [Ggrav] = Ggrav(g,h,r)

    Ggrav = zeros(1,3);
    Ggrav(2) = cos((h/(2*r))*pi)*g;
    Ggrav(1) = sin((h/(2*r))*pi)*g;
end

function h = color_line3(x, y, z, c, varargin)
% color_line3 plots a 3-D "line" with c-data as color
%
%       h = color_line(x, y, z, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line(x, y, z, c, mark) 
%          or
%       h = color_line(x, y, z, c, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x      x-data
%       y      y-data
%       z      z-data
%       c      4th dimension for colouring
%       mark   for scatter plots with no connecting line
%
% out:  h   handle of the surface object

% Copyright (c) 2009, Georg Stillfried
% Copyright (c) 2009, Pekka Kumpulainen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none');
  
if nargin ==5
    switch varargin{1}
        case {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'}
            set(h,'LineStyle','none','Marker',varargin{1})
        otherwise
            error(['Invalid marker: ' varargin{1}])
    end

elseif nargin > 5
    set(h,varargin{:})
end
end
