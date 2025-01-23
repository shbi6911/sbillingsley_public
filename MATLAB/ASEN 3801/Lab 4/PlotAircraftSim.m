%% ASEN 3801, Quadrotor Simulation and Control
% Group Number: 03
% Team Members: Shane Billingsley, Adrian Bryant, Kyle Goodall, Daniela
% Mohammadi
% Date Created: 3/13/2024
% Last Modified: 4/13/2024

% INPUTS
% time: row or column vector of length n (seconds)
% aircraft_state_array: 12 x n matrix of aircraft state at n points in
%      time, all in SI base units
% control_input_array: 4 x n matrix of the four control inputs 
%       [Zc; Lc; Mc; Nc] at n time points, where Zc is the thrust force 
%       (Newtons), and Lc, Mc, & Nc are the contol moments (Newton*meters) 
% fig: row or column vector of length 6 containing integer values of
%       each figure number
% LineColor: specifies line color using RGB triplet, hexadecimal code, or 
%       short name (ex. 'r', 'g', etc.)
% LineSz: specifies line thickness using a numeric value (default is 0.5)
% LabelSz: specifies label font size using a numeric value (default is 11)
%
% OUTPUTS: no variable outputs, produces 6 figures

%% Begin Function Definition
function PlotAircraftSim(time, aircraft_state_array, control_input_array, ...
    fig, LineColor, plotOptions)
% Axis scaling variables
arguments
    time double
    aircraft_state_array (12,:) double
    control_input_array (4,:) double
    fig (1,6) double
    LineColor string
    plotOptions.figLabel string = ''
    plotOptions.TitleSz (1,1) double = 12
    plotOptions.sgTitleSz (1,1) double = 13
    plotOptions.LineSz (1,1) double = 1.5
    plotOptions.LabelSz (1,1) double = 12
    plotOptions.TexLabelSz (1,1) double = 13
    plotOptions.LegendSz (1,1) double = 12
    plotOptions.LineStyle (1,1) string = '-'
    plotOptions.LegendEntries = 'off'
    plotOptions.LegendEntries_3D = 'off'
    plotOptions.xE_xlim = 'tickaligned'
    plotOptions.xE_ylim = 'tickaligned'
    plotOptions.yE_xlim = 'tickaligned'
    plotOptions.yE_ylim = 'tickaligned'
    plotOptions.zE_xlim = 'tickaligned'
    plotOptions.zE_ylim = 'tickaligned'

    plotOptions.phi_xlim = 'tickaligned'
    plotOptions.phi_ylim = 'tickaligned'
    plotOptions.theta_xlim = 'tickaligned'
    plotOptions.theta_ylim = 'tickaligned'
    plotOptions.psi_xlim = 'tickaligned'
    plotOptions.psi_ylim = 'tickaligned'

    plotOptions.u_xlim = 'tickaligned'
    plotOptions.u_ylim = 'tickaligned'
    plotOptions.v_xlim = 'tickaligned'
    plotOptions.v_ylim = 'tickaligned'
    plotOptions.w_xlim = 'tickaligned'
    plotOptions.w_ylim = 'tickaligned'

    plotOptions.p_xlim = 'tickaligned'
    plotOptions.p_ylim = 'tickaligned'
    plotOptions.q_xlim = 'tickaligned'
    plotOptions.q_ylim = 'tickaligned'
    plotOptions.r_xlim = 'tickaligned'
    plotOptions.r_ylim = 'tickaligned'

    plotOptions.xE_3D_xlim = 'tickaligned'
    plotOptions.yE_3D_ylim = 'tickaligned'
    plotOptions.zE_3D_zlim = 'tickaligned'
end

% Breaking out state array values for convenience in reference
xE = aircraft_state_array(1,:);  yE = aircraft_state_array(2,:);
zE = aircraft_state_array(3,:);  phi = aircraft_state_array(4,:);
theta = aircraft_state_array(5,:);  psi = aircraft_state_array(6,:);
u = aircraft_state_array(7,:);  v = aircraft_state_array(8,:);
w = aircraft_state_array(9,:);  p = aircraft_state_array(10,:);
q = aircraft_state_array(11,:); r = aircraft_state_array(12,:);
Zc = control_input_array(1,:);
Lc = control_input_array(2,:);
Mc = control_input_array(3,:);
Nc = control_input_array(4,:); 

% Plotting inertial positions vs. time
figure(fig(1)); 
sgtitle( strcat(plotOptions.figLabel, 'Inertial Positions V_{E}^{E} vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1); hold on; grid on;
        plot(time, xE, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("x_{E}^{E} (m)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.xE_xlim)
            ylim(plotOptions.xE_ylim)

subplot(3,1,2); hold on; grid on;
        plot(time, yE, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("y_{E}^{E} (m)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.yE_xlim)
            ylim(plotOptions.yE_ylim)

subplot(3,1,3); hold on; grid on;
        plot(time, -zE, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("z_{E}^{E} (m)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.zE_xlim)
            ylim(plotOptions.zE_ylim)

% Plotting Euler angles vs. time
figure(fig(2)); 
sgtitle( strcat(plotOptions.figLabel, '3-2-1 Euler Angles vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1); hold on; grid on;
        plot(time, phi, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("\phi (rad)", 'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.phi_xlim)
            ylim(plotOptions.phi_ylim)

subplot(3,1,2); hold on; grid on;
        plot(time, theta, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("\theta (rad)", 'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.theta_xlim)
            ylim(plotOptions.theta_ylim)

subplot(3,1,3); hold on; grid on;
        plot(time,psi, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("\psi (rad)", 'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.psi_xlim)
            ylim(plotOptions.psi_ylim)

% Plotting velocites in body coordinates vs. time
figure(fig(3));
sgtitle( strcat(plotOptions.figLabel,'Body Frame Velocities V_{B}^{E} vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1); hold on; grid on;
        plot(time, u, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time(s)", 'FontSize', plotOptions.LabelSz);
    ylabel("u^{E} (m/s)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.u_xlim)
            ylim(plotOptions.u_ylim)

subplot(3,1,2); hold on; grid on;
        plot(time, v, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("v^{E} (m/s)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.v_xlim)
            ylim(plotOptions.v_ylim)

subplot(3,1,3); hold on; grid on;
        plot(time, w, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time(s)", 'FontSize', plotOptions.LabelSz);
    ylabel("w^{E} (m/s)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.w_xlim)
            ylim(plotOptions.w_ylim)

% Plotting angular velocities in body coordinates vs. time
figure(fig(4));
sgtitle( strcat(plotOptions.figLabel,'Body Frame Rotation Rates \omega_{B}^{E} vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(3,1,1); hold on; grid on;
        plot(time, p, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("$p$ $(\frac{rad}{s})$", 'interpreter', 'latex', ...
        'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.p_xlim)
            ylim(plotOptions.p_ylim)

subplot(3,1,2); hold on; grid on;
        plot(time, q, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("$q$ $(\frac{rad}{s})$", 'interpreter', 'latex', ...
        'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.q_xlim)
            ylim(plotOptions.q_ylim)

subplot(3,1,3); hold on;    grid on;
        plot(time, r, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("$r$ $(\frac{rad}{s})$", 'interpreter', 'latex', ...
        'Rotation', 0, 'FontSize', plotOptions.TexLabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.r_xlim)
            ylim(plotOptions.r_ylim)

% Plotting control inputs vs tine
figure(fig(5));
sgtitle( strcat(plotOptions.figLabel,'Control Force and Moments vs. Time'), ...
    'FontSize', plotOptions.sgTitleSz, 'FontWeight', 'bold')
subplot(2,2,1); hold on; grid on;
        plot(time, Zc, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("Z_{c} (N)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)

subplot(2,2,2); hold on; grid on;
        plot(time, Lc, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("L_{c} (N\cdotm)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)

subplot(2,2,3); hold on; grid on;
        plot(time,Mc, 'Color', LineColor, 'LineWidth',  ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("M_{c} (N\cdotm)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)

subplot(2,2,4); hold on; grid on;
        plot(time, Nc, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
    xlabel("Time (s)", 'FontSize', plotOptions.LabelSz);
    ylabel("N_{c} (N\cdotm)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries, 'FontSize', plotOptions.LegendSz)

% Plotting path in 3D space
figure(fig(6));
hold on; grid on; view(3);
        plot3(xE, yE, -zE, 'Color', LineColor, 'LineWidth', ...
            plotOptions.LineSz, 'LineStyle', plotOptions.LineStyle);
        plot3(xE(1), yE(1), -zE(1), 'g*', 'MarkerSize', 8, 'LineWidth', 1.2)
        plot3(xE(end), yE(end), -zE(end), 'rx', 'MarkerSize', 8, ...
            'LineWidth', 1.4)
    title( strcat(plotOptions.figLabel,"Flight Path in 3-D Space"), ...
        'FontSize', plotOptions.TitleSz);
    xlabel("x_{E}^{E}  (m)", 'FontSize', plotOptions.LabelSz);
    ylabel("y_{E}^{E}  (m)", 'FontSize', plotOptions.LabelSz);
    zlabel("z_{E}^{E} (m)", 'Rotation', 0, 'FontSize', plotOptions.LabelSz);
    legend(plotOptions.LegendEntries_3D, 'FontSize', plotOptions.LegendSz)
        % Plot limits from name-value arguments
            xlim(plotOptions.xE_3D_xlim)
            ylim(plotOptions.yE_3D_ylim)
            zlim(plotOptions.zE_3D_zlim)
end


