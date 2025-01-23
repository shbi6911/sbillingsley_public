%% Question 3.3

Data10_5 = table2array(readtable("24_02_21_001_RWHEEL_10_5"));
Data10_5 = Data10_5(122:593,:); % cutoff values for when motor was producing torque
Data20_10 = table2array(readtable("24_02_21_001_RWHEEL_20_10"));
Data20_10 = Data20_10(130:392,:);
Data20_5 = table2array(readtable("24_02_21_001_RWHEEL_20_5"));
Data20_5 = Data20_5(124:391,:);
Data5_10 = table2array(readtable("24_02_21_001_RWHEEL_5_10"));
Data5_10 = Data5_10(125:1085,:);
Data5_5 = table2array(readtable("24_02_21_001_RWHEEL_5_5"));
Data5_5 = Data5_5(133:600,:);

Data10_5(:,3) = Data10_5(:,3) * 0.1047;  %rpm to rad/s
Data20_10(:,3) = Data20_10(:,3) * 0.1047;
Data20_5(:,3) = Data20_5(:,3) * 0.1047;
Data5_10(:,3) = Data5_10(:,3) * 0.1047;
Data5_5(:,3) = Data5_5(:,3) * 0.1047;

Data10_5(:,4) = Data10_5(:,4) * 0.0335;  %Amps to Nm
Data20_10(:,4) = Data20_10(:,4) * 0.0335;
Data20_5(:,4) = Data20_5(:,4) * 0.0335;
Data5_10(:,4) = Data5_10(:,4) * 0.0335;
Data5_5(:,4) = Data5_5(:,4) * 0.0335;

%Time10_5 = Data10_5(:,1);
%Time20_10 = Data20_10(:,1);
%Time20_5 = Data20_5(:,1);
%Time5_10 = Data5_10(:,1);
%Time5_5 = Data5_5(:,1);


%Plotting Angular velocity over time
figure(1)
subplot(3,2,1)
plot(Data10_5(:,1), Data10_5(:,3))
p10_5 = polyfit(Data10_5(:,1),Data10_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("10 Nm over 5 seconds")

subplot(3,2,2)
plot(Data20_10(:,1), Data20_10(:,3))
p20_10 = polyfit(Data20_10(:,1),Data20_10(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("20 Nm over 10 seconds")

subplot(3,2,3)
plot(Data20_5(:,1), Data20_5(:,3))
p20_5 = polyfit(Data20_5(:,1),Data20_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("20 Nm over 5 seconds")

subplot(3,2,4)
plot(Data5_10(:,1), Data5_10(:,3))
p5_10= polyfit(Data5_10(:,1),Data5_10(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("5 Nm over 10 seconds")

subplot(3,2,5)
plot(Data5_5(:,1), Data5_5(:,3))
p5_5 = polyfit(Data5_5(:,1),Data5_5(:,3),1);
xlabel("time(s)")
ylabel("Angular Velocity (rad/s)")
title("5 Nm over 5 seconds")

%Plotting torque over time , for getting cutoff values
%{
figure(1)
plot(Data10_5(:,1), Data10_5(:,4))
p10_5 = polyfit(Data10_5(:,1),Data10_5(:,3),1);

figure(2)
plot(Data20_10(:,1), Data20_10(:,4))
p20_10 = polyfit(Data20_10(:,1),Data20_10(:,3),1);

figure(3)
plot(Data20_5(:,1), Data20_5(:,4))
p20_5 = polyfit(Data20_5(:,1),Data20_5(:,3),1);

figure(4)
plot(Data5_10(:,1), Data5_10(:,4))
p5_10= polyfit(Data5_10(:,1),Data5_10(:,3),1);

figure(5)
plot(Data5_5(:,1), Data5_5(:,4))
p5_5 = polyfit(Data5_5(:,1),Data5_5(:,3),1);
%}

%MOI = Torque / Angular acceleration
I10_5 = mean(Data10_5(:,4))/p10_5(1);
I20_10 = mean(Data20_10(:,4))/p20_10(1);
I20_5 = mean(Data20_5(:,4))/p20_5(1);
I5_10 = mean(Data5_10(:,4))/p5_10(1);
I5_5 = mean(Data5_5(:,4))/p5_5(1);
I = [I10_5,I20_10,I20_5,I5_10,I5_5];

MeanI = mean(I); %kg m^2
StdDev = std(I);


% Angular Momentum Capacity
%H = I*Omega
H = MeanI*4000*0.1047;  %kg m^2 rad/s
