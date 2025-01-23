%By:        Shane Billingsley
%Class:     ASEN 3712 Structures
%Date:      Fall 2023

%a script to calculate force necessary to elongate a bar in a frame 
% calculated with the small displacement approximation and without

%set some constants
delta = 0:0.001:1;      %range of horizontal displacement
E = 100 * 10^6;         %elastic modulus
A = 0.001;              %cross-sectional area
L_ref = 5;              %reference (initial) length

%calculate force with small displacement approximation
P_small = ((32*E*A).*delta)/(25*L_ref);

%calculate elongation without approximation (stepwise)
L_new=sqrt(9+(4+delta).^2); %find new length of bar
strain = (L_new - L_ref)./L_ref;    %find strain of bar
angle = (4+delta)./L_new;       %find cosine of new angle of bar
P_large = 2*E*A.*strain.*angle; %clacluate force necessary

%find error of small displacement approximation assuming P_large is the
%"correct" value
error = ((P_small - P_large)./P_large).*100;

%plotting results

figure (1); hold on; grid on;
plot(delta,P_small);
plot(delta,P_large);
title ("Displacement of Two-Member Frame Under Load");
xlabel ("Horizontal Displacement (\delta) (m)");
ylabel ("Horizontal Force (P) (N)");
legend ("Small Disp. Approx.", "Not Approximated",'Location','southeast');

figure(2);
plot(delta,error);
title ("Error Between Predicted Forces Under Small Disp. and True Value");
xlabel ("Horizontal Displacement (\delta) (m)");
ylabel ("Percentage Error Between Predicted Forces (P) (N)");