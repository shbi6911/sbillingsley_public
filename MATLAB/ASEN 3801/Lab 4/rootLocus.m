% Function to plot root locus

% K3 is 1x2 and must have the larger (more positive) value as index 1
% K_inner contains the inner loop gains K1 and K2, where K1 is the angular
% derivative gain and K2 is the angular proportional gain

function rootLocus(eigVals, K_inner, K3, dynTypeStr, figNum, plotOptions)
    arguments
        eigVals
        K_inner (1,2) double
        K3 (1,2) double       
        dynTypeStr
        figNum (1,1) double = 1
        plotOptions.ptSz (1,1) double = 5
        plotOptions.lineSz (1,1) double = 0.5
        plotOptions.TitleSz (1,1) double = 12
        plotOptions.LabelSz (1,1) double = 12
        plotOptions.LegendSz (1,1) double = 12

    end
% Extracting inner-loop gains for title
    K1= K_inner(1); K2= K_inner(2);

% Splitting eigenvalues into real & imaginary components
    re= real(eigVals);
    im= imag(eigVals);

figure(figNum)
hold on; grid on; axis padded
    title([dynTypeStr ' Dynamics Root Locus with K_{1}=' num2str(K1) ...
        ' & K_{2}=' num2str(K2)], 'FontSize', plotOptions.TitleSz)
    xlabel('Re.', 'FontSize', plotOptions.LabelSz)
    ylabel('Im.', 'FontSize', plotOptions.LabelSz, 'Rotation', 0)

        plot(re, im, '*', 'MarkerSize', plotOptions.ptSz, 'LineWidth', ...
            plotOptions.lineSz)
        plot(re(1,:),im(1,:), 'g*', 'MarkerSize', plotOptions.ptSz, ...
            'LineWidth', plotOptions.lineSz)
        plot(re(end,:),im(end,:), 'r*', 'MarkerSize', plotOptions.ptSz, ...
            'LineWidth', plotOptions.lineSz)

    legend('','','', ['Initial K_{3}=' num2str(K3(1))], ...
        ['Final K_{3}=' num2str(K3(2))], 'FontSize', plotOptions.LegendSz)
hold off


end