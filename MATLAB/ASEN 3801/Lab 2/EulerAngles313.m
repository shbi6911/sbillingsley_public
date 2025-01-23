% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% EulerAngles313
% Created: 1/31/24

function attitude313 = EulerAngles313(DCM)
%INPUTS     DCM         a 3x3 Direction Cosine Matrix for a rotation
%
%OUTPUTS    attitude313 a 3x1 vector of Euler angles associated with the
%                       rotation in the order 3-1-3
%
%METHODOLOGY    EulerAngles313 converts a given DCM for a rotation into
%Euler angles in the 3-1-3 sequence.
alpha = atan2d(DCM(3,1),-DCM(3,2));
beta = -asind(DCM(3,3));
gamma = atan2d(DCM(1,3),DCM(2,3));
attitude313 = [alpha beta gamma];
end