%By:        Shane Billingsley
%Class:     ASEN 2804 Vehicle Design Lab
%Date:      Spring 2023

%define design constants

%lengths in mm
L = 320;        %length of bottle
nose = 140;     %additional length of nosecone
rear = 80;      %distance rear of wing from nozzle end

%aerospike length
spike = linspace(0,500,501);

%total length front of spike to rear of wing
length = ((L+nose)-rear)+spike;

%Mach angle
mach = 2;
mu = asin(1/mach);

%calculate wingspans and plot
span = tan(mu)*length*2;

plot(spike,span,'LineWidth',2);
title ("Aerospike Length vs. Max Wingspan");
txt = ['Bottle Length = ',int2str(L),' mm, Nosecone Length = ',...
    int2str(nose),' mm, Wing ',int2str(rear),' mm from nozzle'];
subtitle(txt);
xlabel("Length of Aerospike (mm)");
ylabel("Maximum Wingspan At Rear (mm)");





