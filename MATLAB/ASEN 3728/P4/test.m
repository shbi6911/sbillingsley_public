% Define the center point
center_x = 200;
center_y = 200;

% Define the side length
side_length = 5;

% Calculate the coordinates of the square
x1 = center_x - side_length/2;
y1 = center_y - side_length/2;
width = side_length;
height = side_length;

% Plot the square
rectangle('Position', [195, 195, 10, 10], 'EdgeColor', 'b', 'LineWidth', 2);

% Set axis limits and aspect ratio
axis([center_x - side_length, center_x + side_length, center_y - side_length, center_y + side_length]);
axis equal; % Equal aspect ratio
grid on; % Show grid