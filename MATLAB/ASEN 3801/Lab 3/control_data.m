%% Control Data

clear; close all;

data = readmatrix("24_03_06_Control_76.47_26.76",'FileType','text',NumHeaderLines=2);
clean_data = data(:,1:4);
clean_data(:,1) = clean_data(:,1) - clean_data(1,1);
clean_data(:,1) = clean_data(:,1)./1000;

figure(); hold on; grid on;
plot(clean_data(:,1),clean_data(:,2));
plot(clean_data(:,1),clean_data(:,3));
legend("Reference Height","Measured Position","Location",'northoutside');
title("Controlled Output with K1 = 76.46, K2 = 26.76");
xlabel("Time in Seconds");
ylabel("Angular Position (radians)");
hold off;

subset = find(clean_data(:,1)>=10 & clean_data(:,1) <= 22);
clean_data_subset = clean_data(subset,:);
figure(); hold on; grid on;
plot(clean_data_subset(:,1),clean_data_subset(:,2));
plot(clean_data_subset(:,1),clean_data_subset(:,3));
legend("Reference Height","Measured Position","Location",'northoutside');
set(gca,'YDir','reverse');
title("Controlled Output with K1 = 76.46, K2 = 26.76");
xlabel("Time in Seconds");
ylabel("Angular Position (radians)");
hold off;
