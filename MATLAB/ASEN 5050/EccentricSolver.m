%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function E = EccentricSolver(a,mu,e,delta_t,tol)
%EccentricSolver uses Newton's method to solve Kepler's equation to find
%eccentric anomaly given an elapsed time past periapsis
%INPUTS     a           semi-major axis km
%           mu          gravitational parameter km^3/s^2
%           e           eccentricity
%           delta_t     t-tp, time past periapsis at a point in the orbit
%                       in seconds
%           tol         tolerance to converge to
%
%OUTPUTS    E           eccentric anomaly at the time point given by
%                       delta_t in radians
n = sqrt(mu/a^3);   %find mean motion
period = (2*pi)/n;  %find orbit period
if delta_t > period         %correct delta_t in case it is more than one period
    numperiods = fix(delta_t/period);
    T = delta_t - (period*numperiods);
else
    T = delta_t;
end %period corrector

M = n*T;                %find mean anomaly
if (-pi <= M && M <=0) || (pi <= M && M < (2*pi))%form initial guess
    E = M - e;
else
    E = M + e;
end %initial guess if

Kepler = E - e*sin(E) - M;      %find g(E) based on initial guess
iterator = 0;
while abs(Kepler) >= tol       %iterate E until g(E) ~ 0
    E = E + ((M-E+e*sin(E))/(1 - e*cos(E)));
    iterator = iterator +1;
    if iterator >= 100
        break
    end                      %break if
    Kepler = E - e*sin(E) - M;
end %while loop

end %function