%calculates all forms of stability and creates a trim diagram
clc;close all;
%constants
x_cg = 250.832*10^-3; %m
x_ac_w = 224.1*10^-3; %m
x_ac_vt = 388*10^-3; %m
span = 500.6*10^-3; %m
taper_ratio = 1;
chord_wing = 60*10^-3; %m
wing_span = 425*10^-3; %m
planform_wing = (18490.5*2)*10^-6; %m^2
planform_Htail = (8327.025*2)*10^-6; %m^2
planform_vtail = (3714.8*2)*10^-6; %m^2
C_mac = 0; %N
weight = 0.255; %kg
rho = 1.2250;  %kg/m^3 %currently S.L
ep_partial = 0.35; %approx
a0 = 0.11; %found using naca 0006
oswalds_eff = 0.665; %using model data
it = 1; %can change
L_D_max = 4.63;
glide_angle = 1/L_D_max;

%const equations
AR = (span^2)/planform_wing;
lt = x_ac_vt - x_cg; %m
a = a0/(1+a0/(pi * oswalds_eff * AR)); %per degree
volume_H = (lt * planform_Htail) / (chord_wing * planform_wing); 
volume_V = (planform_vtail * (x_ac_vt - x_cg)) / (planform_wing * wing_span);
neutral_pt = (x_ac_w/chord_wing) + volume_H * (1-ep_partial); 
h_cg = x_cg/chord_wing;
static_margin = neutral_pt - (x_cg/chord_wing); 

%place holders
angle_of_attack = zeros(20,1);
C_m_cg = zeros(20,1); 

%variable eqs
for i = 1:20
    angle_of_attack(i,1) = i-10; %deg
    C_m_cg(i,1) = a*angle_of_attack(i,1)*(((x_cg/chord_wing) - (x_ac_w/chord_wing)) - volume_H*(1 - ep_partial)) + (volume_H * a * it) ;
end
%graph
plot(C_m_cg,angle_of_attack,'Linewidth',2)
title('moment coefficient vs Alpha');
xlabel('alpha(degrees)');
ylabel('moment coefficient');