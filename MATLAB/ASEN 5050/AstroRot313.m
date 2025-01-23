%By:    Shane Billingsley
%Class: ASEN 5050 Space Flight Dynamics
%Date:  Fall 2024

function DCM = AstroRot313(o, i, t)
%AstroRot313 constructs a DCM to rotate an astrodynamical vector from
%inertial XYZ axes into r theta h axes
%
%INPUT      three angles in radians
%   o       right ascension of the ascending node
%   i       inclination
%   t       arg of periapsis + true anomaly   
co = cos(o);    ci = cos(i);    ct = cos(t);
so = sin(o);    si = sin(i);    st = sin(t);

DCM = [((co*ct)-(so*ci*st)),((-co*st)-(so*ci*ct)),(so*si);
       ((so*ct)+(co*ci*st)),((-so*st)+(co*ci*ct)),(-co*si);
       (si*st), (si*ct), (ci)];
end