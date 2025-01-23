% Shane Billingsley, Gabriel Law, Sean McCluskey
% ASEN 3801
% LoadASPENData
% Created: 1/31/24

function [t_vec, av_pos_inert, av_att, tar_pos_inert, tar_att] = LoadASPENData(filename)
%INPUTS     filename        a string variable of the filename to be loaded
%
%OUTPUTS    t_vec           1xn time vector in seconds
%           av_pos_inert    3xn matrix of position vector of vehicle in 
%                           meters in inertial frame (E)
%           av_att          3xn matrix of 3-2-1 Euler angles of vehicle in
%                           radians in inertial frame (E)
%           tar_pos_inert   3xn matrix of position vector of target in 
%                           meters in inertial frame (E)
%           tar_att         3xn matrix of 3-2-1 Euler angles of target in
%                           radians in inertial frame (E)
%METHODOLOGY    LoadASPENData takes in data from a previously cleaned csv
%file of data taken in the ASPEN facility.  See lab documentation for data
%sanitation instructions.  Note that ConvertASPENData must be used to
%convert given helical angles to Euler angles
data = readmatrix(filename);    %read in data
t_vec = data(:,1)'./100;             %time vector in seconds
%note, all position data is in mm, angular data is in radians
tar_att_temp = data(:,2:4)';    %target attitude, helical angles N frame
tar_pos_temp = data(:,5:7)';    %target position, N frame
av_att_temp = data(:,8:10)';    %vehicle attitude, helical angles N frame
av_pos_temp = data(:,11:13)';   %vehicle position, N frame

%remove data values where either target has NaN 
nan_cols = any(isnan(av_att_temp) | isnan(tar_att_temp), 1); %find NaNs
t_vec(:, nan_cols) = [];
av_att_temp(:, nan_cols) = [];
tar_att_temp(:, nan_cols) = []; %trim arrays indexed by nan_cols
av_pos_temp(:, nan_cols) = [];
tar_pos_temp(:, nan_cols) = [];
%convert data into proper format using ConvertASPENData
[av_pos_inert, av_att, tar_pos_inert, tar_att] = ConvertASPENData...
    (av_pos_temp,av_att_temp,tar_pos_temp,tar_att_temp);

%convert angle arrays to degrees
av_att = rad2deg(av_att);
tar_att = rad2deg(tar_att);


end