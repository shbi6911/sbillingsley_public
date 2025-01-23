%% ASEN 3801, Quadrotor Simulation and Control: Task 3
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/29/2024
% Last Modified: 4/7/2024


% K1 = angular rate gain (p, q)
% K2 = angular gain (phi, theta)
% K3 = inertial velocity gain (v, u)

% Function to calculated real eigenvalues and associated gains

function [eigVals, gains] = latGains(tau1, tau2, K3_max, K3_min, numPts)
    % Initializing necessary parameters
        [~, g, ~, ~, I, ~, ~]= getParams();
    % Moment of inertia about body x-axis
        Ix= I(1,1); % [kg*m^2]
    % Allocating memory
        gains= zeros(numPts,3); eigVals=gains;
    % Obtaining K2 & K1 from time constants and storing
        [K2, K1]= timeConst2gains(tau1, tau2, Ix);
        gains(:,1)= K1*ones(numPts,1); % roll rate gain
        gains(:,2)= K2*ones(numPts,1); % roll angle gain

    % Vector of K3 values to test
        K3_vec=linspace(K3_max, K3_min, numPts)';
        gains(:,3)= K3_vec; % 

% Looping through eigenvalues resulting from variation of K3
    for i=1:numPts
    % Full 4x4 system
        Acl= [0, 1, 0, 0;...
              0, 0, g, 0;...
              0, 0, 0, 1;...
            0, -K3_vec(i)/Ix, -K2/Ix, -K1/Ix] ;
    % Reducing to 3x3 system for obtaining eigenvalues
        Acl= Acl(2:4, 2:4);
    % Getting eigenvalues
        eigVals(i,:)= eig(Acl)'; 
    end

end








