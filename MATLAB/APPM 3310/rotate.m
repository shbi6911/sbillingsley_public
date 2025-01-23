%By:        Shane Billingsley
%Class:     APPM 3310 Matrix Methods and Applications
%Date:      Spring 2023

function [outvec] = rotate(theta,axis,invec)
%ROTATE This function rotates a vector by an angle theta
%   Inputs: Angle (this is the SCALAR rotation angle)
%           Axis (this is a unit COLUMN VECTOR in the direction of the axis)
%           Vector (this is the COLUMN VECTOR to be rotated)
%   Outputs:    Outvec (this is the new vector after rotation)    
    axis = axis./norm(axis); %normalize axis vector

    nmatrix = axis*axis';           %n*n^T
    crossmatrix = [0 -axis(3) axis(2); axis(3) 0 -axis(1); %cross product
                    -axis(2) axis(1) 0];
    a = cos(theta);
    b = 1-a;
    c= sin(theta);
    %define axis-angle rotation matrix
    rotatematrix = (b*nmatrix)+(a*eye(3))+(c*crossmatrix);
    outvec = rotatematrix*invec; %mulitply for resultant vector
end