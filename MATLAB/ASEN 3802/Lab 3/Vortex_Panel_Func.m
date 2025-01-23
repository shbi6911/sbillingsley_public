%% HOUSEKEEPING
clc
clear
close all

%% DECLARING VARIABLES

VINF = 30;
ALPHA = linspace(-20,20,41);
ALPHA_rad = ALPHA.*(pi/180);

%% GATHERING NACA DATA FOR Q1

A = readmatrix("NACA 0012 5.txt");
XB_0012_5 = [A(:,1)]';
YB_0012_5 = [A(:,2)]';
A = readmatrix("NACA 0012 10.txt");
XB_0012_10 = [A(:,1)]';
YB_0012_10 = [A(:,2)]';
A = readmatrix("NACA 0012 20.txt");
XB_0012_20 = [A(:,1)]';
YB_0012_20 = [A(:,2)]';
A = readmatrix("NACA 0012 30.txt");
XB_0012_30 = [A(:,1)]';
YB_0012_30 = [A(:,2)]';
A = readmatrix("NACA 0012 40.txt");
XB_0012_40 = [A(:,1)]';
YB_0012_40 = [A(:,2)]';
A = readmatrix("NACA 0012 50.txt");
XB_0012_50 = [A(:,1)]';
YB_0012_50 = [A(:,2)]';
A = readmatrix("NACA 0012 60.txt");
XB_0012_60 = [A(:,1)]';
YB_0012_60 = [A(:,2)]';
A = readmatrix("NACA 0012 1000.txt");
XB_0012_1000 = [A(:,1)]';
YB_0012_1000 = [A(:,2)]';

%% CALCULATING CL FOR Q1
g = 1;
for a = 12
CL_0012_5(g) = Vortex_Panel(XB_0012_5,YB_0012_5,VINF,a);
CL_0012_10(g) = Vortex_Panel(XB_0012_10,YB_0012_10,VINF,a);
CL_0012_20(g) = Vortex_Panel(XB_0012_20,YB_0012_20,VINF,a);
CL_0012_30(g) = Vortex_Panel(XB_0012_30,YB_0012_30,VINF,a);
CL_0012_40(g) = Vortex_Panel(XB_0012_40,YB_0012_40,VINF,a);
CL_0012_50(g) = Vortex_Panel(XB_0012_50,YB_0012_50,VINF,a);
CL_0012_60(g) = Vortex_Panel(XB_0012_60,YB_0012_60,VINF,a);
CL_0012_1000(g) = Vortex_Panel(XB_0012_1000,YB_0012_1000,VINF,a);
fprintf("Iteration " + num2str(g)+"\n")
g = g+1;
end

fprintf("NACA 0012 CL @ AOA = 12 is " + num2str(Vortex_Panel(XB_0012_50,YB_0012_50,VINF,12)) + "\n")




%% PLOTTING Q1

figure(1);
plot(5,CL_0012_5,"*",'MarkerSize',7)
hold on
plot(10,CL_0012_10,"*",'MarkerSize',7)
grid on
plot(20,CL_0012_20,"*",'MarkerSize',7)
plot(30,CL_0012_30,"*",'MarkerSize',7)
plot(40,CL_0012_40,"*",'MarkerSize',7)
plot(50,CL_0012_50,"*",'MarkerSize',7)
plot(60,CL_0012_60,"*",'MarkerSize',7)
yline(CL_0012_1000)
yline(CL_0012_1000*.99,'--r')
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 0012 with Varying Amount of Panels")
legend("5 Panels","10 Panels","20 Panels","30 Panels","40 Panels","50 Panels","60 Panels","1000 Panels","1000 Panels 1 Percent Threshold",'location','southeastoutside')
saveas(gcf,"NACA 0012 with Varying Amount of Panels")

%% GATHERING NACA DATA FOR Q2

A = readmatrix("NACA 0006 50.txt");
XB_0006 = [A(:,1)]';
YB_0006 = [A(:,2)]';
A = readmatrix("NACA 0012 50.txt");
XB_0012 = [A(:,1)]';
YB_0012 = [A(:,2)]';
A = readmatrix("NACA 0018 50.txt");
XB_0018 = [A(:,1)]';
YB_0018 = [A(:,2)]';



%% CALCULATING CL

g = 1;
for a = ALPHA
CL_0006(g) = Vortex_Panel(XB_0006,YB_0006,VINF,a);
CL_0012(g) = Vortex_Panel(XB_0012,YB_0012,VINF,a);
CL_0018(g) = Vortex_Panel(XB_0018,YB_0018,VINF,a);
CL_TAT_0006(g) = ThinAirfoilTheory_symmetric(a);
CL_TAT_0012(g) = ThinAirfoilTheory_symmetric(a);
CL_TAT_0018(g) = ThinAirfoilTheory_symmetric(a);
fprintf("Iteration " + num2str(g)+"\n")
g = g+1;
end

B = polyfit(ALPHA_rad,CL_0006,1);
dcl_0006 = B(1);
a_0_0006 = B(2);
B = polyfit(ALPHA_rad,CL_0012,1);
dcl_0012 = B(1);
a_0_0012 = B(2);
B = polyfit(ALPHA_rad,CL_0018,1);
dcl_0018 = B(1);
a_0_0018 = B(2);

figure(10)
table([a_0_0006;a_0_0012;a_0_0018],[dcl_0006;dcl_0012;dcl_0018],'VariableNames',["Zero Lift Angle of Attack","Lift Slope"],'RowNames',["NACA 0006";"NACA 0012";"NACA 0018"])
saveas(gcf,"A0 vs dcl.png")
%% GATHERING WIND TUNNEL DATA

NACA0006_dat = readmatrix("NACA 0006 DATA.csv");
NACA0012_dat = readmatrix("0012_cl_alpha_data.txt");

%% PLOTTING Q2

figure(2)
plot(ALPHA,CL_0006)
hold on
grid on
plot(ALPHA,CL_TAT_0006)
plot(NACA0006_dat(:,1),NACA0006_dat(:,2))
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 0006")
legend("Vortex Panel Method","Thin Aifroil Theory","Wind Tunnel Data","Location","southeastoutside")
saveas(gcf,"CL_0006 Q2.png")

figure(3)
plot(ALPHA,CL_0012)
hold on
grid on
plot(ALPHA,CL_TAT_0012)
plot(NACA0012_dat(:,1),NACA0012_dat(:,2))
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 0012")
legend("Vortex Panel Method","Thin Aifroil Theory","Wind Tunnel Data","Location","southeastoutside")
saveas(gcf,"CL_0012 Q2.png")

figure(4)
plot(ALPHA,CL_0018)
hold on
grid on
plot(ALPHA,CL_TAT_0018)
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 0018")
legend("Vortex Panel Method","Thin Aifroil Theory","Location","southeastoutside")
saveas(gcf,"CL_0018 Q2.png")

%% GETTING DATA FOR Q3

A = readmatrix("NACA 2412 50.txt");
XB_2412 = [A(:,1)]';
YB_2412 = [A(:,2)]';
A = readmatrix("NACA 4412 50.txt");
XB_4412 = [A(:,1)]';
YB_4412 = [A(:,2)]';

%% CALCULATING CL

g = 1;
for a = ALPHA
CL_2412(g) = Vortex_Panel(XB_2412,YB_2412,VINF,a);
CL_4412(g) = Vortex_Panel(XB_4412,YB_4412,VINF,a);

fprintf("Iteration " + num2str(g)+"\n")
g = g+1;
end

B = polyfit(ALPHA_rad,CL_2412,1);
dcl_2412 = B(1);
a_0_2412 = Cambered_A_0(CL_2412,ALPHA);
B = polyfit(ALPHA_rad,CL_4412,1);
dcl_4412 = B(1);
a_0_4412 = Cambered_A_0(CL_4412,ALPHA)

table([a_0_2412;a_0_4412],[dcl_2412;dcl_4412],'VariableNames',["Zero Lift Angle of Attack","Lift Slope"],'RowNames',["NACA 2412";"NACA 4412"])

g = 1;
for a = ALPHA

    CL_TAT_2412(g) = TAT_cambered(a,a_0_2412);
    CL_TAT_4412(g) = TAT_cambered(a,a_0_4412);
    fprintf("Iteration " + num2str(g)+"\n")
    g = g+1;

end

%% GATHERING DATA FROM BOOK 

NACA4412_dat = readmatrix("NACA 4412 DATA.csv");
NACA2412_dat = readmatrix("2412_cl_alpha_data.txt");


%% PLOTTING Q3

figure(5)
plot(ALPHA,CL_0012)
hold on
grid on
plot(ALPHA,CL_TAT_0012)
plot(NACA0012_dat(:,1),NACA0012_dat(:,2))
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 0012")
legend("Vortex Panel Method","Thin Aifroil Theory","Wind Tunnel Data","Location","southeastoutside")
saveas(gcf,"CL_0012 Q3.png")

figure(6)
plot(ALPHA,CL_2412)
hold on
grid on
plot(ALPHA,CL_TAT_2412)
plot(NACA2412_dat(:,1),NACA2412_dat(:,2))
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 2412")
legend("Vortex Panel Method","Thin Aifroil Theory","Wind Tunnel Data","Location","southeastoutside")
saveas(gcf,"CL_2412 Q3.png")

figure(7)
plot(ALPHA,CL_4412)
hold on
grid on
plot(ALPHA,CL_TAT_4412)
plot(NACA4412_dat(:,1),NACA4412_dat(:,2))
xlabel("AOA (Degrees)")
ylabel("Cl")
title("NACA 4412")
legend("Vortex Panel Method","Thin Aifroil Theory","Wind Tunnel Data","Location","southeastoutside")
saveas(gcf,"CL_4412 Q3.png")

%% FUNCTIONS

function [CL] = ThinAirfoilTheory_symmetric(ALPHA)

ALPHA_rad = ALPHA * (pi/180);
CL = 2*pi*ALPHA_rad;

end

function [A_0] = Cambered_A_0(CL,ALPHA)

ALPHA_rad = ALPHA*(pi/180);

A_0 = ALPHA_rad - (CL/(2*pi));
A_0 = (sum(A_0)/length(A_0))*(180/pi);

end

function [CL] = TAT_cambered(ALPHA,A0)

ALPHA_rad = ALPHA * (pi/180);
A0_rad = A0 * (pi/180);

CL = 2*pi*(ALPHA_rad-A0_rad);

end

function [CL] = Vortex_Panel(XB,YB,VINF,ALPHA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                           %
%                                  %
% XB  = Boundary Points x-location %
% YB  = Boundary Points y-location %
% VINF  = Free-stream Flow Speed   %
% ALPHA = AOA                      %
%                                  %
% Output:                          %
%                                  %
% CL = Sectional Lift Coefficient  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%

ALPHA = ALPHA*pi/180;

%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%

CHORD = max(XB)-min(XB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships:                                  %
%                                                             %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I = 1:M
    IP1 = I+1;
    X(I) = 0.5*(XB(I)+XB(IP1));
    Y(I) = 0.5*(YB(I)+YB(IP1));
    S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
    THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
    SINE(I) = sin( THETA(I) );
    COSINE(I) = cos( THETA(I) );
    RHS(I) = sin( THETA(I)-ALPHA );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:             %
%                                        %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I = 1:M
    for J = 1:M
        if I == J
            CN1(I,J) = -1.0;
            CN2(I,J) = 1.0;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
        else
            A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
            B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
            C = sin( THETA(I)-THETA(J) );
            D = cos( THETA(I)-THETA(J) );
            E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
            F = log( 1.0 + S(J)*(S(J)+2*A)/B );
            G = atan2( E*S(J), B+A*S(J) );
            P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
              + (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
            Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
              - (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
            CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
            CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
            CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
            CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:           %
%                                      %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I = 1:M
    AN(I,1) = CN1(I,1);
    AN(I,MP1) = CN2(I,M);
    AT(I,1) = CT1(I,1);
    AT(I,MP1) = CT2(I,M);
    for J = 2:M
        AN(I,J) = CN1(I,J) + CN2(I,J-1);
        AT(I,J) = CT1(I,J) + CT2(I,J-1);
    end
end
AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;
for J = 2:M
    AN(MP1,J) = 0.0;
end
RHS(MP1) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%

GAMA = AN\RHS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for I = 1:M
    V(I) = cos( THETA(I)-ALPHA );
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
    end
    CP(I) = 1.0 - V(I)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION/CHORD;
end