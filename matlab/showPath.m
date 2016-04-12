function showPath(shollfile, pathfile, varargin)
% Usage: showPath(shollfile, pathfile)
%
% Input:
%   shollfile is an N x 4 space-separated (x,y,z,color) txt file.
%   pathfile is an N x 3 space-separated (x,y,z) txt file or str to file.
%
% Given the complete sholl color file from sholl_color in
% neuron_getProperties.py and the path points from axon_path in same file,
% this shows the path from the soma to the axon tip superimposed on the
% sholl color plot.
%


% Inputs
if isstr(shollfile) ~= 1
  m = size(shollfile);
  if m(2) ~= 4
    disp('The sholl file must either be a string for a file (Nx4, txt)')
    disp(' or an Nx4 (x,y,z,col) matrix')
  else
    sholl = shollfile;
  end
else
  sholl = importdata(shollfile);
end
if isstr(pathfile) ~= 1
  m = size(pathfile);
  if m(2) ~= 3
    disp('The pathfile must either be a string for a file (Nx3, txt)')
    disp(' or an Nx3 (x,y,z) matrix')
  else
    path = pathfile;
  end
else
  path = importdata(pathfile);
end

if nargin > 2
  somapos = varargin{1};
  shollColor(sholl, somapos);
else
  % Plotting
  shollColor(sholl);
end


hold on;
scatter3(path(:,1), path(:,2), path(:,3), 20, 'ko', 'filled')
axis('equal');
set(gca, 'visible','off');


end