function [gpts] = showEllipseFit(pts, reverse)
% usage: showEllipseFit(pts)
%
% Shows the scaled grid points to the fitted ellipse; use removePoints 
%   first (before calling this).
%
% Inputs:
%   pts - the points produced by neuron_resampleGrid.py
%   reverse - (optional) which axis to reverse (needed sometimes)
%          ex: [-1], or [-1,-1] or [-1,1,-1]
%
% Outputs:
%   gpts - the adjusted scaled grid points
%   plot of the gridpoints and fit points (pts)
%

% parse inputs
if nargin > 1
  reverse = reverse;
else
  reverse = 0;
end

[cent, ax, evecs, v] = ellipsoid_fit( pts );

gpts = getEllipseGrid(ax);

% now scale and translate gpts to be closer to data
if reverse ~= 0
  evecs(:,1) = evecs(:,1)*reverse(1);
end
if numel(reverse) > 1
  evecs(:,2) = evecs(:,2)*reverse(2);
end
if numel(reverse) > 2
  evecs(:,3) = evecs(:,3)*reverse(3);
end

gpts = gpts*evecs;

% translate gpts to centroid
for i = 1:3
  gpts(:,i) = gpts(:,i) + cent(i);
end

figure()
hold on;
scatter3(gpts(:,1), gpts(:,2), gpts(:,3), 'r.')
scatter3(pts(:,1), pts(:,2), pts(:,3), 'bo', 'FaceColor', 'b')
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
axis equal;
