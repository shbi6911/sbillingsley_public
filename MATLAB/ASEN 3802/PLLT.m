function [e,c_l,c_di] = PLLT(const,N)
%PLLT applies the Fundamental Equation of Prandtl's Lifting Line Theory to
%a given airfoil, defined by input parameters
%   INPUTS: const:  A STRUCT containing the following fields:
%
%               const.b:      scalar value of wingspan (unit distance)
%               (degrees per unit distance)
%               const.a0_t:   scalar cross-sectional lift slope at wingtip
%               const.a0_r:   scalar cross-sectional lift slope at wing root
%               (unit distance)
%               const.c_t:    scalar value of chord at wingtip
%               const.c_r:    scalar value of chord at wing root
%               (degrees)
%               const.aero_t: scalar value of zero-lift angle of attack at
%               wingtip
%               const.aero_r: scalar value of zero-lift angle of attack at wing
%               root
%               (degrees)
%               const.geo_t:  scalar value of geometric angle of attack at the
%               wing tip
%               const.geo_r:  scalar value of geometric angle of attack at the
%               wing root
%           N:      scalar number of ODD terms to use in the Fourier
%               series expansion within PLLT
%
%   OUTPUTS:    e:      scalar value of span efficiency factor
%               c_l:    scalar value of coefficient of lift for the given
%               airfoil geometry
%               c_di:   scalar value of induced drag coefficient for the
%               given airfoil geometry

Nvec = linspace(1, (2*N)-1, N);     %linearly spaced vector of odd values of N
theta = (Nvec*pi)./(2*Nvec(end));       %corresponding angular locations along span
AR = (const.b^2)/const.S;      %calculate aspect ratio

%convert input values to radians
%a0_t = a0_t*(pi/180); a0_r = a0_r*(pi/180);
aero_t = const.aero_t*(pi/180); aero_r = const.aero_r*(pi/180);
geo_t = const.geo_t*(pi/180); geo_r = const.geo_r*(pi/180);

%get vectors of values at theta locations assuming linear variation along span
a0_vec = span(const.a0_t,const.a0_r,theta); 
c_vec = span(const.c_t,const.c_r,theta);
aero_vec = span(aero_t,aero_r,theta);
geo_vec = span(geo_t,geo_r,theta);

%generate matrix to solve for Fourier series coefficients (A_n)
%this matrix  is being generated as the sum of two matrices, one for each
%term in the Fundamental Equation of PLLT

%generate a matrix of size NxN of terms of the form n*theta
ntheta = ones(length(Nvec),length(theta));
ntheta = theta'.*ntheta;
ntheta = Nvec.*ntheta;

%generate a matrix of size NxN of terms of the form (4b/(a_0*c))sin(ntheta)
coeff_vec = (4*const.b)./(a0_vec.*c_vec); %leading term of Fund. Eq. of PLLT
mat1 = coeff_vec' .* sin(ntheta); %first coefficient matrix

%generate a matrix of size NxN of terms of the form
%n(sin(ntheta)/sin(theta))
mat2 = Nvec.*sin(ntheta);
mat2 = (1./sin(theta))' .* mat2;

%add these matrices and solve the system
finalmat = mat1 + mat2;
A_nvec = finalmat\(geo_vec - aero_vec)';

%calculate outputs using the Fourier series coefficients
e = 1/(1+sum(((A_nvec(2:end)./A_nvec(1)).^2).*Nvec(2:end)'));
c_l = A_nvec(1)*pi*AR;
c_di = (c_l^2)/(pi*e*AR);

function [y] = span(t, r, theta)
    % SPAN outputs an array of values at y locations, assuming linear
    % variation, decreasing from r at pi/2 to t at zero.
    %   INPUTS: t:      scalar value at theta=0
    %           r:      scalar value at theta=pi/2
    %           theta:  vector of angular locations between 0 and pi/2 to get values
    %                   assumed to be between zero and pi/2, increasing
    %   OUTPUTS:    span:   vector of values at y locations
    y = (t-r).*cos(theta) + r;  %calculate values
    
end

end