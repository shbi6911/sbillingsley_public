%By:    Shane Billingsley
%Class: ASEN 6060 Advanced Astrodynamics
%Date:  Spring 2024

%% Define constants
function const = getConst(verbose)
%getConst defines a struct of useful constants.  This iteration provides
%constants useful to the Circular Restricted 3 Body Problem in the context
%of ASEN 6060 homework and projects
%INPUTS:    verbose     a Boolean flag (0 or 1) to control whether output
%                       statements such as plotting and printing are
%                       executed
%OUTPUTS:   const       structure of constant values, either given or
%                       derived
%
    arguments
        verbose logical =0
    end
    %set verbose flag to control output (plots & print statements)
    const.vrb = verbose;
    % gravitational parameter
    const.mu.Earth = 3.98600435507e5;
    const.mu.Moon = 4902.800118;
    % const.mu.Mars = 4.305e4;        %km^3/s^2
    const.mu.Sun = 1.32712440041279419e11;
    % const.mu.Saturn = 3.794e7;
    % const.mu.Jupiter = 1.268e8;
    
    %average planetary radii
    const.radius.Earth = 6378.1363;
    const.radius.Moon = 1738;
    % const.radius.Mars = 3397.2;     %km
    % const.radius.Saturn = 60268;
    % const.radius.Jupiter = 71492;
    
    %semi-major axes
    const.a.Earth = 149598023;   %km, Sun-relative
    const.a.Moon = 384400;      %km, Earth-relative
    % const.a.Mars = 1.52;
    % const.a.Saturn = 9.554909595;
    % const.a.Jupiter = 5.202603191;
    
    %eccentricities
    const.e.Earth = 0.016708617;
    const.e.Moon = 0.05490;
    
    % const.AU = 149597870.7;         %km
    const.G = 6.67408e-20;           %km^3/(kg*s^2)

    %CR3BP parameters
    % characteristic masses
    const.m_star.em = (const.mu.Earth/const.G) + (const.mu.Moon/const.G);
    const.m_star.se = (const.mu.Sun/const.G) + (const.mu.Earth/const.G);
    % mass ratios
    const.mu.em = (const.mu.Moon/const.G)/const.m_star.em;
    const.mu.se = (const.mu.Earth/const.G)/const.m_star.se;
    % find nondimensional length quantities
    const.l_star.em = const.a.Moon;
    const.l_star.se = const.a.Earth;
    % find nondimensional time quantities
    const.t_star.em = sqrt((const.l_star.em^3)/(const.G*const.m_star.em));
    const.t_star.se = sqrt((const.l_star.se^3)/(const.G*const.m_star.se));

    %equilibrium point locations
    const.eq_pts.em = eqPts(const.mu.em);
    const.eq_pts.se = eqPts(const.mu.se);

end