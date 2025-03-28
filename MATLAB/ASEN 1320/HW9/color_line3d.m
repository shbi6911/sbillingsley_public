%By:        Shane Billingsley
%Class:     ASEN 1320 Aerospace Computing and Engineering Applications
%Date:      Fall 2021

function h = color_line3d(c, x, y)
% color_line3 plots a 3-D "line" with c-data as color
%
%       color_line3d(c, x, y)
%
%  in:  x      x-data
%       y      y-data
%       c      coloring
%

z = zeros(size(x));

h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','interp',...
  'EdgeColor','interp',...
  'Marker','none', ...
  'LineWidth',2);

colorbar

end
