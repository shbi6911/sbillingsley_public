function [F, alpha, beta] = Lambers_Equation_Solver_Function(dtheta_star_corrected, flightpath, a, s, c, mu, tof)

n = sqrt(mu / a^3);

alpha_0 = 2 * asin(sqrt(s / (2 * a)));
beta_0 = 2 * asin(sqrt((s - c) / (2 * a)));

    if dtheta_star_corrected < 180 && flightpath == 0
        alpha = alpha_0;
        beta = beta_0;
    elseif dtheta_star_corrected < 180 && flightpath == 1
        alpha = 2*pi - alpha_0;
        beta = beta_0;
    elseif dtheta_star_corrected > 180 && flightpath == 0
        alpha = alpha_0;
        beta = -beta_0;
    elseif dtheta_star_corrected > 180 && flightpath == 1
        alpha = 2*pi - alpha_0;
        beta = -beta_0;
    end
    
%     F = alpha - beta - (sin(alpha) - sin(beta)) - tof * n;
    F = 1/n * (alpha - beta - (sin(alpha) - sin(beta))) - tof;
end