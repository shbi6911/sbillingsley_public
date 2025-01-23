function const = getConst()
%getConst defines a struct of values describing the physical situation of a
%finite wing.  Geometry of the wing as well as some aspects of the flight
%condition are specified
%INPUTS:        None
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere 

    constPLLT = getConstPLLT();         %constants for PLLT
    %use PLLT to define constants for convergence, as well as some others
    constSYS = getConstSYS(constPLLT);
    %these lines concatenate the two preceding structs into one
    fn1 = fieldnames(constPLLT);
    fn2 = fieldnames(constSYS);
    fn = [fn1; fn2];
    c1 = struct2cell(constPLLT);
    c2 = struct2cell(constSYS);
    c = [c1;c2];
    const = cell2struct(c,fn,1);

function const = getConstPLLT()
%getConstPLLT is a subfunction that defines a struct of constants used by
%the PLLT function, in the context of Question 5 of the ASEN 3802 lab 3
%report.  Thus the wing geometry defined herein is that required by that
%section.
%INPUTS:        None
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere 

    % geometry of the wing
    const.b = 33 + (4/12);    %wingspan in feet
    const.c_r = 5 + (4/12);   %root chord in feet
    const.c_t = 3 +(8.5/12);  %tip chord in feet
    const.S = (const.b/2)*(const.c_r+const.c_t);    %planform area
    
    %root chord airfoil is a NACA 2412
    const.a0_r = 2*pi;    %lift slope in rad^-1
    const.aero_r = -1;    %zero-lift angle of attack in degrees
    const.wing_geo_r = 1;      %geometric angle of attack in degrees
    %tip chord airfoil is a NACA 0012
    const.a0_t = 2*pi;    %lift slope in rad^-1
    const.aero_t = 0;    %zero-lift angle of attack in degrees
    const.wing_geo_t = 0;      %geometric angle of attack in degrees

    const.aot = 4;    %angle of attack the wing is flying in degrees
     %add this angle to the geo angles of the wing
     const.geo_r = const.wing_geo_r + const.aot;
     const.geo_t = const.wing_geo_t + const.aot;
end

function const = getConstSYS(constPLLT)
%getConstSYS is a subfunction that defines a struct of constants used by
%the errorfinder function, in the context of Question 5 of the ASEN 3802 lab 3
%report.  Thus the wing geometry defined herein is that required by that
%section.  This subfunction is separate from getConstPLLT because
%getConstSYS must call the PLLT function in order to define convergence
%values.
%INPUTS:        constPLLT:  A struct of constants needed by the PLLT
%                           function
%OUTPUTS:       A struct containing various fields used by the PLLT
%               function, the errorfinder function, and elsewhere
 
    % values of the physical system
    const.vel = 100 *1.68781;   %velocity of the plane in ft/sec
    const.altitude = 10000;       %altitude in feet
    
    %find atmospheric values
    %convert altitude to meters to use atmoscoesa function
    const.altitude_m = const.altitude * 0.3048;     %altitude in meters
    %find pressure and density at given altitude
    [~,~,const.P,const.rho] = atmoscoesa(const.altitude_m); 
    const.P = const.P* 0.0208854342; %convert Pascals to lb/ft^2
    const.rho = const.rho * 0.00194032; %convert kg/m^3 to slug/ft^3
   
    %find lift and drag coefficients with a very large N
    %assume these are the convergence values
    [const.e_abs,const.c_l_abs,const.c_di_abs] = PLLT(constPLLT,1000);
    
    %calculate convergence values for lift and induced drag using coefficients
    const.q = 0.5*const.rho*const.vel^2;  %dynamic pressure in Pa
    const.lift_abs = const.c_l_abs*const.q*constPLLT.S;     %total lift
    const.drag_i_abs = const.c_di_abs*const.q*constPLLT.S;  %total induced drag
end

end