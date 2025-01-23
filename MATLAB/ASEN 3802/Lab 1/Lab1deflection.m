%this script calculates maximum beam deflection and plots against
%experimental data

%IMPORTANT the script assumes that script Lab1data has already been run

E=69.5*10^9;    %elastic modulus in Pa
r_2=d/2;        %outer radius
r_1=r_2-t;      %inner radius
A=(pi*r_2^2)-(pi*r_1^2);                %area of a bar in m^2
I=4*((pi/4)*(r_2^4-r_1^4)+(A*(l/2)^2));  %area moment of inertia in m^4

F_A=meandata.Case_1_4(:,4); %extract relevant data

%calculate deflection for Case 1 at point L/2
coeff_1=1/(E*I);
coeff_2=F_A./6;
coeff_3=(F_A*L^2)./8;
deflection_1=coeff_1.*((coeff_2.*(L/2).^3)-(coeff_3.*(L/2)));

%plotting
figure (4);
hold on;
scatter(meandata.Case_1_4(:,1),meandata.Case_1_4(:,6));
scatter(meandata.Case_2_6(:,1),meandata.Case_2_6(:,6));
scatter(meandata.Case_3_4(:,1),meandata.Case_3_4(:,6));
scatter()
title ("Load Case 1 (Single Central Loading)");
xlabel ("Total Load Weight in Newtons (N)");
ylabel ("Measured Displacement at Center");
legend("Case 1 Dataset 4","Case 2 Dataset 6","Case 3 Dataset 4",'location','southeast');
hold off;