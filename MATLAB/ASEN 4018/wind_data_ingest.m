%By:    Shane Billingsley
%Class: ASEN 4018 Senior Projects
%Date:  Fall 2024

%% ingest wind data

rows = [80;60;40;20;0;-20;-40;-60;-80]; %latitude values as rows of data matrix
%longitude values as columns of data matrix
columns = [-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180];
data = readmatrix('wind.xlsx'); %read in data
%parse out angle data
angles = data(:, 2:2:end);
%parse out magnitude data
magnitudes = data(:,1:2:end);

%convert to cell array
angles_cell = num2cell(angles);
%loop through array to convert angles to unit vectors and multiply by
%magnitudes
arr_size = size(angles_cell);
for i = 1:arr_size(1)
    for j = 1:arr_size(2)
        angle = angles_cell{i,j};
        vector = [-sind(angle);cosd(angle);0];
        vector_2 = vector.*magnitudes(i,j);
        rot_matrix = WindRot(rows(i),columns(j));
        angles_cell{i,j} = rot_matrix'*vector_2;
    end
end

