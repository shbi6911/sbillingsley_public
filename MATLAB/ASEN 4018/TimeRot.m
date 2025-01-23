%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

function rotated_array = TimeRot(input_array, delta_t)
%TimeRot takes in a cell array of column vectors in Earth Fixed axes
%(x-axis on prime meridian, y-axis at 90 deg longitude, z axis north), and
%rotates this array's vectors by a time delta_t based on Earth's rotation
%
% INPUTS        input_array     cell array of column vectors
%               delta_t         time of rotation in seconds
% OUTPUTS       rotated_array   the same cell array of column vectors,
%                               rotated by an angle based on delta_t
%
rate = (2*pi)/86164;    %rotation rate rad/s based on sidereal day
angle = delta_t*rate;   %angle to rotate
%rotate about the z-axis by angle
rot_matrix = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
%preallocate
arr_size = size(input_array);
rotated_array = cell(size(input_array));
%loop through input array and rotate each vector
for i = 1:arr_size(1)
    for j = 1:arr_size(2)
        rotated_array{i,j} = rot_matrix*input_array{i,j};
    end
end

end
