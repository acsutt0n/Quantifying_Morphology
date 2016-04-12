function shollAxons(shollfile, axonfile)
% Usage: shollAxons(shollfile, axonfile)
% Input:
%   shollfile is an N x 4 space-separated (x,y,z,color) txt file.
%   axonfile is an N x 7 list of potential axon branches from
%     axons_endpoints in neuron_getProperties.py. Can be str or matrix.
%     Format is: branchNum   x0    y0    z0    x1    y1    z1
%                   int    float float float float float float
%
% This function will plot the neuron and the possible branch numbers so the
% user can decide which is the real axon and send it back to
% neuron_getProperties for analysis.


% Inputs
if ischar(shollfile) ~= 1
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
if ischar(axonfile) ~= 1
  m = size(axonfile);
  if m(2) ~= 7
    disp('The axonfile must either be a string for a file (Nx4, txt)')
    disp(' or an Nx7 (branchNum, x0, y0, z0, x1, y1, z1) matrix')
  else
    axons = axonfile;
  end
else
  axons = importdata(axonfile);
end


% Plot 
[m, ~] = size(axons);
shollColor(sholl);
hold on;
for i = 1:m
  str = int2str(axons(i,1));
  text(axons(i,2),axons(i,3), axons(i,4), str);
  text(axons(i,5),axons(i,6), axons(i,7), str);
end
axis('equal');
set(gca, 'visible','off');


end
