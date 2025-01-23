%By:        Shane Billingsley
%Class:     ASEN 3700 Attitude Dynamics & Orbital Mechanics
%Date:      Spring 2024
%
function quaternion = Quat(axis,angle)
%Quat generates a quaternion from an input rotation axis and angle
%
%INPUTS     axis    vector representing axis of rotation
%           angle   scalar value of rotation angle in radians
%OUTPUTS    quaternion      associated quaternion for the rotation,
%                           represented as a 4-element column vector
x = sin(angle/2);
if isrow(axis)
    axis = axis';
end
quaternion = [x.*axis;cos(angle/2)];
end