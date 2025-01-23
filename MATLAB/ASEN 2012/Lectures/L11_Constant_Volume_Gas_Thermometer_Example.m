% Housekeeping
clear all;
clc;
close all;

% Example 8.5 - Constant-Volume Gas Thermometer Example using the matrix
% formulation

% Pressure is the x variable,
% Make sure that data are stored as a column vector
x = [65; 75; 85; 95; 105];

% Temperature is the y variable,
% Make sure that data are stored as a column vector
y = [-20; 17; 42; 94; 127];

% Assume that uncertainity in each data point is the same
% Thus we can find least squares paramenters without the weighting matrix.
% We will then correct for uncertainty in the parameter using the computed
% uncertaity for the weights

% For the least squares 
% set up H-matrix 
H = [ones(length(y),1),x];

% Carryout the least squares without the weighting matrix
Par = ((H' * H)^-1) * H' * y;
A = Par(1);
B = Par(2);

% Errors in the data
Deviation = y - (A + B.*x);
Sigma_y = sqrt((1/(length(y)-length(Par))) * sum(Deviation .* Deviation));

% Using the error in the data as the weights in the weighting matrix so
% that we can determine the error covariance matrix to get the uncertainty
% in the parameters
delta_y (1:length(y)) = Sigma_y;
Diagonal = 1 ./ (delta_y .* delta_y);
W = diag(Diagonal);

% Error covariance matrix
Sigma_Par_2 = (H' * W * H)^-1;

% Uncertainty in paramenter A 
Sigma_A = sqrt(Sigma_Par_2(1,1));

% Double check to show that the fit parameters are are the same when using
% the weighting matrix
Par_1 = ((H' * W * H)^-1) * H' * W * y;

% Predicted values at using the fit equation
y_fit = A + B .* x;

% Extrapolation using the least sqaures fit
% Extrapolation to find the temaprature at a pressure of zero 
y_0 = A + B * 0; % we can see this would be the same as the y-intercept A

% Extrapolation using the matrix formulation
y_0_new = [1 0] * Par;

% Extraploation of the uncertainty to pressure of zero
% If we do this correctly we should get the same uncertainty value as the
% uncertainty in parameter A since A is the y-intercept
Sigma_y0 = sqrt ([1 0]* Sigma_Par_2 *[1; 0]);

% Calculate sigma_q for extrapolated and fit data
xextrap = linspace(0,x(5),25)';
Hnew = [ones(length(xextrap),1),xextrap];
Sigma_q = diag(sqrt(Hnew * Sigma_Par_2 * Hnew'));

% Create Figure 8.4
figure(); hold on; grid on;
errorbar(x,y,delta_y,'o','LineWidth',2)
plot(x,A+B*x,'k','LineWidth',2)
plot(xextrap,A+B*xextrap,'k--','LineWidth',2)
errorbar(0,A,Sigma_A,'o','LineWidth',2,'Color',[0.85	0.325	0.098])
plot([xextrap],A+B*[xextrap]+delta_y(1),'--','Color',[0	0.447 0.741],'LineWidth',1.2)
plot([xextrap],A+B*[xextrap]-delta_y(1),'--','Color',[0	0.447 0.741],'LineWidth',1.2, 'HandleVisibility','off')
plot([xextrap],A+B*[xextrap]+Sigma_q,'--','Color',[0.85	0.325	0.098],'LineWidth',1.2)
plot([xextrap],A+B*[xextrap]-Sigma_q,'--','Color',[0.85	0.325	0.098],'LineWidth',1.2, 'HandleVisibility','off')
legend('data','fit','extrapolation','Abs. Zero Est.','\sigma_T','\sigma_q','Location','SE')
xlabel('P (mmHg)')
ylabel('T (deg. C)')