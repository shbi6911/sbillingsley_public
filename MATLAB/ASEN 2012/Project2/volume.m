%By:        Shane Billingsley
%Class:     ASEN 2012 Experimental and Computational Methods
%Date:      Fall 2022

function vol = volume(t,v,c)
A_t = pi*c.bottle.d_e^2*(1/4);
sub_v = (c.params.P0*((c.params.vol0/v)^c.const.gamma))-c.const.P_a;
vol = c.bottle.c_dis*A_t*sqrt((2/c.const.rho_water)*sub_v);
if vol < c.bottle.V
    vol = 0;
end