clear
clc
close all
cmap=colormap(lines(10));

%% Initialize random spread of data around line of y=4x+2
x=linspace(-10,10,101)';
spread = 25;
rng('default'); rng(1); %comment out this line to repeat with new data each time
y=(-x.^2 + 2*x+4) + spread*randn(size(x)); %create new random data
figure(1);set(1,'Units','inches','Position',[ 1,1 ,10,6],'Color','w');hold on;
plot(x,y,'o','Color',cmap(1,:))
xlabel('x');ylabel('y');
legend('data','Location','EastOutside');

%% Calcualte A and B using matrix method
H = [ones(size(x)) x x.^2]; %create H for linear fit
xHat = inv(H'*H)*H'*y; %solve for xHat (weights are not necessary)
A=xHat(1); %intercept
B=xHat(2); %x term
C=xHat(3); %square term

%% Plot the line of best fit on top of the data
yfit=C*x.^2+B*x+A; %line of best fit
plot(x,yfit,'Color',cmap(2,:),'LineWidth',2)
legend('data','best fit line','Location','EastOutside');

%% Calculate sigma_y (uncertainty in y data) and plot 95% confidence interval
sigma_y = sqrt(1/(length(x)-2) * sum((yfit-y).^2)); %uncertainty in y data
plot(x,yfit+2*sigma_y,'--','Color',cmap(5,:),'LineWidth',2) %plotting 95% confidence interval of y data
plot(x,yfit-2*sigma_y,'--','Color',cmap(5,:),'LineWidth',2, 'HandleVisibility','off') %plotting 95% confidence interval of y data
legend('data','best fit line','95% confidence \sigma_y','Location','EastOutside');

%% Calculate sigma_A, sigma_B and sigma_q (uncertainty in A and B terms and of the line of bets fit) using matrices and plot 95% confidence interval of sigma_q
W=1/(sigma_y^2)*eye(length(x)); %create weighting matrix
P=inv(H'*W*H); %calculate covariance matrix
sigma_A1 = sqrt(P(1,1)); %extract out sigma_A
sigma_B1 = sqrt(P(2,2)); %extract out sigma_B
sigma_C1 = sqrt(P(3,3)); %extract out sigma_C
sigma_q2 = diag(sqrt(H*P*H')); %calculate sigma_q based on sigma_q^2 = HPH^T
plot(x,yfit+2*sigma_q2,'--','Color',cmap(3,:),'LineWidth',2) %plotting 95% confidence interval of linear fit data
plot(x,yfit-2*sigma_q2,'--','Color',cmap(3,:),'LineWidth',2, 'HandleVisibility','off') %plotting 95% confidence interval of linear fit data
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - mat.','Location','EastOutside');

%% plot the "exact" line
plot(x,-x.^2 + 2*x+4,'-k','LineWidth',2)
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - mat.','y=2x+4','Location','EastOutside');

%% use polyfit and polyval
[p,S]=polyfit(x,y,2); 
% [P,S] = polyfit(X,Y,N) returns the polynomial coefficients P and a
%     structure S for use with POLYVAL to obtain error estimates for
%     predictions.  S contains fields for the triangular factor (R) from a
%     QR decomposition of the Vandermonde matrix of X, the degrees of
%     freedom (df), and the norm of the residuals (normr).  If the data Y
%     are random, an estimate of the covariance matrix of P is
%     (Rinv*Rinv')*normr^2/df, where Rinv is the inverse of R.

% Y = polyval(P,X) returns the value of a polynomial P evaluated at X. P
%     is a vector of length N+1 whose elements are the coefficients of the
%     polynomial in descending powers:
%  
%         Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
%  
%     The polynomial P is evaluated at all points in X.

plot(x,polyval(p,x),'--','Color',cmap(6,:),'LineWidth',2)
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - mat.','y=-x^2+2x+4','polyval','Location','EastOutside');

polyfitCovar = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
disp('polyfit covariance: ')
fprintf('%d %d %d \n',[polyfitCovar]);
disp('matrix method covariance: ')
fprintf('%d %d %d\n',[P]);