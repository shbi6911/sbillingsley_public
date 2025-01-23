%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/25/2024
% Last Modified: 4/13/2024

%% Begin Function Definition
function plotMiniDrone(x, t, refMat, fig, plotOptions)
    arguments
        x double
        t double
        refMat double
        fig (1,3) double
        plotOptions.figLabel string = ''
        plotOptions.LineColor string = "#0072BD"
        plotOptions.sgTitleSz (1,1) double = 13
        plotOptions.LineSz (1,1) double = 1.5
        plotOptions.LabelSz (1,1) double = 12
        plotOptions.TexLabelSz (1,1) double = 13
        plotOptions.LegendSz (1,1) double = 12
    end

% Extracting inertial translation and rotation variables
    % Inertial positions in meters
        xE= x(:,1); yE= x(:,2); zE= x(:,3);
    % Roll, Pitch, & Yaw angles in radians
        psi= x(:,4); theta= x(:,5); phi= x(:,6);
% Extracting reference positions and angles to compare to in plots
    % Inertial reference positions in meters
        xE_ref= refMat(:,1); yE_ref= refMat(:,2); zE_ref= refMat(:,3);
    % Reference Roll, Pitch, & Yaw angles in radians
        psi_ref= refMat(:,4); theta_ref= refMat(:,5); phi_ref=refMat(:,6);

%% Task 1, Question 1.5: Plotting Minidrone Flight Path & Angles
% Checking data for NaNs prior to plotting
    nanCheck= sum(isnan([t, xE, yE, zE]));
if sum(nanCheck)==0 % Plotting if no NaNs in data

% Subplot of inertial positions vs. time
figure(fig(1))
sgtitle( strcat(plotOptions.figLabel, 'Minidrone Inertial Positions V_{E}^{E} vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1)
    hold on; grid on;
        xlabel('Time (s)', 'FontSize', plotOptions.LabelSz)
            xlim([t(1) t(end)])
        ylabel('x_{E}^{E} (m)', 'Rotation', 0, 'FontSize', plotOptions.LabelSz)
            plot(t, xE, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(xE_ref,'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Position', 'FontSize',plotOptions.LabelSz)
subplot(3,1,2)
    hold on; grid on;
        xlabel('Time (s)', 'FontSize', plotOptions.LabelSz)
            xlim([t(1) t(end)])
        ylabel('y_{E}^{E} (m)', 'Rotation', 0, 'FontSize', plotOptions.LabelSz)
            plot(t, yE, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(yE_ref, 'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Position', 'FontSize',plotOptions.LabelSz)
subplot(3,1,3)
    hold on; grid on;
        xlabel('Time (s)', 'FontSize', plotOptions.LabelSz)
            xlim([t(1) t(end)])
        ylabel('z_{E}^{E} (m)', 'Rotation', 0, 'FontSize', plotOptions.LabelSz)
            plot(t, -zE, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(-zE_ref, 'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Position', 'FontSize',plotOptions.LabelSz)

% Subplot of 3-2-1 Euler Angles vs. time
figure(fig(2))
sgtitle( strcat(plotOptions.figLabel, 'Minidrone Euler Angles vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1)
    hold on; grid on;
        xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
            xlim([t(1) t(end)])
        ylabel("\phi (rad)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
            plot(t, phi, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(phi_ref, 'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Angle', 'FontSize',plotOptions.LabelSz)
subplot(3,1,2)
    hold on; grid on;
        xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
            xlim([t(1) t(end)])
        ylabel("\theta (rad)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
            plot(t, theta, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(theta_ref,'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Angle', 'FontSize',plotOptions.LabelSz)
subplot(3,1,3)
    hold on; grid on;
        xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
            xlim([t(1) t(end)])
        ylabel("\psi (rad)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
            plot(t, psi, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            yline(psi_ref, 'g-.', 'LineWidth', plotOptions.LineSz)
        legend('Actual Position', 'Reference Angle', 'FontSize',plotOptions.LabelSz)

% 3-D trajectory of minidrone
    figure(fig(3))
    hold on; grid on; view(3)
        title( strcat(plotOptions.figLabel, 'Minidrone Flight Path in 3-D Space'), ...
            'FontSize', plotOptions.sgTitleSz)
        xlabel('(m)', 'FontSize', plotOptions.LabelSz)
        ylabel('(m)', 'Rotation', 0, 'FontSize', plotOptions.LabelSz)
        zlabel('(m)', 'Rotation', 0, 'FontSize', plotOptions.LabelSz)
            plot3(xE, yE, -zE, 'LineWidth', plotOptions.LineSz, 'Color', plotOptions.LineColor)
            plot3(xE(1), yE(1), -zE(1), 'g*', 'MarkerSize', 8)
            plot3(xE(end), yE(end), -zE(end), 'r*', 'MarkerSize', 8)
        legend('Flight Path', 'Initial Position', 'Final Position', ...
            'FontSize', plotOptions.LabelSz)

else
    disp('Data contains NaNs')
end % end statement for conditional

end % end statement for function