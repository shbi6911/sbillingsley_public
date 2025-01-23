clear
clc
close all
cmap=colormap(lines(10));

%% Initialize random spread of data around line of y=4x+2
x=linspace(-50,50,101)';
spread = 25;
rng('default'); %this and next lines keep "random" data to be repeatable
rng(1); %comment out this line to repeat with new data each time
y=(2*x+4) + spread*randn(size(x)); %create new random data
figure(1);set(1,'Units','inches','Position',[ 1,1 ,10,6],'Color','w');hold on;
plot(x,y,'o','Color',cmap(1,:))
xlabel('x');ylabel('y');
legend('data','Location','EastOutside');

%% Calculate A and B using summation form
Delta = length(x)*sum(x.^2) - (sum(x)).^2;
Asum = (sum(x.^2)*sum(y) - sum(x)*sum(x.*y))/Delta;
Bsum = (length(x)*sum(x.*y) - sum(x)*sum(y))/Delta;

%% Calcualte A and B using matrix method
H = [ones(size(x)) x]; %create H for linear fit
xHat = inv(H'*H)*H'*y; %solve for xHat (weights are not necessary)
A=xHat(1); %intercept
B=xHat(2); %slope

%% Compare A and B from various methods!
checkA = round(A,6)==round(Asum,6);
checkB = round(B,6)==round(Bsum,6);
if(checkA&&checkB==1)
    disp('The methods are equivalent for A and B!')
else
    disp('The methods are not equivalent for A and B!')
end

%% Plot the line of best fit on top of the data
yfit=B*x+A; %line of best fit
plot(x,yfit,'Color',cmap(2,:),'LineWidth',2);
legend('data','best fit line','Location','EastOutside');

%% Calculate sigma_y (uncertainty in y data) and plot 95% confidence interval
sigma_y = sqrt(1/(length(x)-2) * sum((yfit-y).^2)); %uncertainty in y data
plot(x,yfit+2*sigma_y,'--','Color',cmap(5,:),'LineWidth',2) %plotting 95% confidence interval of y data
plot(x,yfit-2*sigma_y,'--','Color',cmap(5,:),'LineWidth',2, 'HandleVisibility','off') %plotting 95% confidence interval of y data
legend('data','best fit line','95% confidence \sigma_y','Location','EastOutside');
%% Test sigma_y function and compare
[sigma_ytesting, W_test] = sigma_y_test(x,y,xHat);
checky = round(sigma_y,6)==round(sigma_ytesting,6);

if(checky(:)==1)
    disp('The methods are equivalent for sigma_y!')
else
    disp('The methods are not equivalent for sigma_y!')
    disp(sigma_y)
    disp(sigma_ytesting)
end


%% Calculate sigma_A, sigma_B and sigma_q (uncertainty in A and B terms and of the line of bets fit) algebraically and plot 95% confidence interval of sigma_q
sigma_A = sigma_y*sqrt(sum(x.^2)/Delta); %calculate uncertainty in A term
sigma_B = sigma_y*sqrt(length(x)/Delta); %calculate uncertainty in B term
sigma_q = sqrt(sigma_A^2 + (x*sigma_B).^2 ); %calculate uncertainty in the line of best fit
plot(x,yfit+2*sigma_q,'Color',cmap(4,:),'LineWidth',2) %plotting 95% confidence interval of linear fit data
plot(x,yfit-2*sigma_q,'Color',cmap(4,:),'LineWidth',2, 'HandleVisibility','off') %plotting 95% confidence interval of linear fit data
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - alg.','Location','EastOutside');

%% Calculate sigma_A, sigma_B and sigma_q (uncertainty in A and B terms and of the line of bets fit) using matrices and plot 95% confidence interval of sigma_q
P=covariance(x,W_test); %calculate covariance matrix
sigma_A1 = sqrt(P(1,1)); %extract out sigma_A
sigma_B1 = sqrt(P(2,2)); %extract out sigma_B
sigma_q2 = sigma_q_test(x,P); %calculate sigma_q based on sigma_q^2 = HPH^T
plot(x,yfit+2*sigma_q2,'--','Color',cmap(3,:),'LineWidth',2) %plotting 95% confidence interval of linear fit data
plot(x,yfit-2*sigma_q2,'--','Color',cmap(3,:),'LineWidth',2, 'HandleVisibility','off') %plotting 95% confidence interval of linear fit data
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - alg.','95% confidence \sigma_q - mat.','Location','EastOutside');

%% compare sigma_q vectors
checkq = round(sigma_q,6)==round(sigma_q2,6);

if(checkq(:)==1)
    disp('The methods are equivalent for sigma_q!')
else
    disp('The methods are not equivalent for sigma_q!')
end

%% plot the "exact" line
plot(x,2*x+4,'-k','LineWidth',2)
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - alg.','95% confidence \sigma_q - mat.','y=2x+4','Location','EastOutside');

%% use polyfit and polyval
[p,S]=polyfit(x,y,1); 
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
legend('data','best fit line','95% confidence \sigma_y','95% confidence \sigma_q - alg.','95% confidence \sigma_q - mat.','y=2x+4','polyval','Location','EastOutside');

polyfitCovar = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
disp('polyfit covariance: ')
fprintf('%d %d \n',[polyfitCovar]);
disp('matrix method covariance: ')
fprintf('%d %d \n',[P]);

function P = covariance(x,W)
N = length(x);
H = [ones(N,1) x];
P = inv(H'*W*H);
end

function [sigma_y, W] = sigma_y_test(x,y,fit)
sigma_y = sqrt(1/(length(x)-2) * sum(((fit(2)*x+fit(1))-y).^2));
W = (1/sigma_y^2)*eye(length(x));
end

function sigma_q = sigma_q_test(x,P)
N = length(x);
H = [ones(N,1) x];
sigma_q = diag(sqrt(H * P * H'));
end