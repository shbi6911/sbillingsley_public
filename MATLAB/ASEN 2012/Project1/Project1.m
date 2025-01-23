%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

%% read in and clean data
% opts = detectImportOptions("SampleB");
% disp(opts);
data = readmatrix("SampleB");
time = data(:,2);
thermo1 = data(:,3); %cal1
thermo2 = data(:,4); %water temp
thermo3 = data(:,5); %room temp
thermo4 = data(:,6); %cal2

%% average calorimeter thermocouples and determine error
thermo_cal = mean([thermo1 thermo4],2); %cal1 and cal2
thermo_cal_sigma_x = sqrt(((thermo1 - thermo_cal)+(thermo4 - thermo_cal).^2));
thermo_cal_err = thermo_cal_sigma_x./sqrt(2);

%% plot initial data
figure (1); hold on;
% plot (time,thermo1,"red");
%scatter (time,thermo2,"blue","+");
% %plot (time,thermo3,"black");
% plot (time,thermo4,"magenta");
scatter (time,thermo_cal,"green","+");
% legend("Thermo1","Thermo4","Thermo_cal");

%disp(find(time == 442.461));
%disp(find(time == 681.44));
%disp(find(time == 500.739));
%disp(find(time == 439.625));

%% split arrays
time_line1 = time(1:312);
thermo_cal_line1 = thermo_cal(1:312);
time_line2 = time(480:746);
thermo_cal_line2 = thermo_cal(480:746);
time_line3 = time(315:350);
thermo_cal_line3 = thermo_cal(315:350);

%plot(time_line1,thermo_cal_line1,'red');
%plot(time_line2,thermo_cal_line2,'blue');
%plot(time_line3,thermo_cal_line3,'green');

%%find least-squares fit lines
line1_fit = leastSquaresLine(time_line1,thermo_cal_line1);
line2_fit = leastSquaresLine(time_line2,thermo_cal_line2);
line3_fit = leastSquaresLine(time_line3,thermo_cal_line3);

line = @(fit,x) fit(2)*x + fit(1);

plot(0:1200,line(line1_fit,0:1200),'red','LineWidth',1);
plot(0:1200,line(line2_fit,0:1200),'red','LineWidth',1);
plot(440:525,line(line3_fit,440:525),'red','LineWidth',1);
xline(439,'LineStyle',"--");

%% find values for equation temperatures
t0 = 439;  %time when sample is added
T_L = line(line1_fit,439);
T_H = line(line2_fit,439);
T_avg = (T_L+T_H)/2;         %avg temperature
t_avg = (1/line3_fit(2)*T_avg)-(line3_fit(1)/line3_fit(2)); %time @ T_avg

%% define all equation values
T0 = T_L;            %temp of calorimeter at t0
T1 = thermo2(310);    %temp of sample at t0
T2 = line(line2_fit,t_avg);  %temp of calorimeter at equilibrium
m_c = 0.51*1000; % mass of calorimeter in g
c_c = 0.895; %specific heat of calorimeter in J/(g*degC)
m_s = 39.306; %mass of sample in g

%% calculate specific heat of sample
c_s = (m_c*c_c*(T2-T0))/(m_s*(T1-T2));
disp("C_s");
disp(c_s);

yline(T2,'LineStyle',"--");
yline(T_avg,'LineStyle',"--");
xline(t_avg,'LineStyle',"--");
%disp(T_L);
%disp(T_H);
%disp(T_avg);
%disp(t_avg);
%disp(T2);

%% compute partial derivatives
m_c_deriv = (c_c*(T2-T0))/(m_s*(T1-T2));
m_s_deriv = -(m_c*c_c*(T2-T0))/((m_s^2)*(T1-T2));
T0_deriv = -(m_c*c_c)/(m_s*(T1-T2));
T1_deriv = -(m_c*c_c*(T2-T0))/(m_s*(T1-T2)^2);
T2_deriv = (m_c*c_c*(T1-T0))/(m_s*(T1-T2)^2);

%% define error values
m_c_err = 0.05; %error value of calorimter mass in grams
m_s_err = 0.001; %error value of sample mass in grams

%% calculate sigma_y for each fit line
[sigma_y1,W1] = sigma_y(time_line1,thermo_cal_line1,line1_fit);
[sigma_y2,W2] = sigma_y(time_line2,thermo_cal_line2,line2_fit);
[sigma_y3,W3] = sigma_y(time_line3,thermo_cal_line3,line3_fit);
P1 = covariance(time_line1,W1);
P2 = covariance(time_line2,W2);
P3 = covariance(time_line3,W3);

%% find sigma_q for each line
time_line1_extrap = (0:450)';
sigma_q1 = sigma_q(time_line1_extrap,P1);
time_line2_extrap = (400:1200)';
sigma_q2 = sigma_q(time_line2_extrap,P2);
sigma_q3 = sigma_q(time_line3,P3);

%%  find error of T_L and T_H
disp ("2*sigma_y for T_L")
disp(2*sigma_y1)
disp ("sigma_q for T_L")
disp(2*sigma_q1(t0))
disp ("2*sigma_Y for T_H")
disp (2*sigma_y2)
disp ("sigma_q for T_H")
disp(2*sigma_q2(39))
if (2*sigma_y1) >= (2*sigma_q1(t0))
    T_L_err = sigma_y1;
else
    T_L_err = sigma_q1(t0);
end

if (2*sigma_y2) >= (2*sigma_q2(39))
    T_H_err = sigma_y2;
else
    T_H_err = sigma_q2(39);
end

%%  find error of T0 and T2
T0_err = T_L_err;
T_avg_err = norm([T_L_err T_H_err]);
t_avg_err = norm([((1/line3_fit(2))*T_avg_err),... 
   ((-1/line3_fit(2))*sqrt(P3(1,1))),...
   (((T_avg-line3_fit(1))/line3_fit(2)^2)*sqrt(P3(2,2)))]);
T2_err = norm([(line2_fit(2)*t_avg_err),(t_avg*sqrt(P2(2,2))),(sqrt(P2(1,1)))]);

%% find error of T1
thermo2fit = leastSquaresLine(time(1:310),thermo2(1:310));
[sigma_y_thermo2,~] = sigma_y(time(1:310),thermo2(1:310),thermo2fit);
T1_err = 2*sigma_y_thermo2;

%% calculate final error value
c_s_err = norm([(m_c_err*m_c_deriv),(T2_err*T2_deriv),(T0_err*T0_deriv),...
    (m_s_err*m_s_deriv),(T1_err*T1_deriv)]);
disp("error in C_s")
disp(c_s_err);

function fit = leastSquaresLine(x,y)
N = length(x);
H = [ones(N,1) x];
fit = inv(H'*H)*(H'*y);
end

function P = covariance(x,W)
N = length(x);
H = [ones(N,1) x];
P = inv(H'*W*H);
end

function [sigma_y, W] = sigma_y(x,y,fit)
sigma_y = sqrt(1/(length(x)-2) * sum(((fit(2)*x+fit(1))-y).^2));
W = (1/sigma_y^2)*eye(length(x));
end

function sigma_q = sigma_q(x,P)
N = length(x);
H = [ones(N,1) x];
sigma_q = diag(sqrt(H * P * H'));
end