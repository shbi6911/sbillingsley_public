% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% RotationMatrix313
% Created: 1/31/24

function attitude321 = EulerAngles321(DCM)
% INPUTS:   DCM:    3x3 3-2-1 Rotation Matrix
%           
% OUTPUTS:  attitude312: 3x1 array of angles [deg]

alpha = atan2d(DCM(2,3),DCM(3,3));  %convert to Euler angles using formulas
beta = -asind(DCM(1,3));
gamma = atan2d(DCM(1,2),DCM(1,1));
attitude321 = [alpha; beta; gamma];
end


