function [newpts, epts, ax] = plotEllipseGrid(pts)
% usage: plotEllipseGrid(pts)
%
% newpts = scaled gridpoints that will fit ellipsoid
% epts = ellipsoid surface points
% this calls ellipsoid_fit.m and makes a grid to plot the ellipse
% around the origin

% get the evecs and evals (vectors and magnitudes of maj/min axes)
[center, ax, evecs, v] = ellipsoid_fit(pts);

% set radii for plotting
r = ax; % previously: r = ax/2 was too small

% subtract mean from points to center them at origin
oldpts = pts;
pts(:,1) = pts(:,1)-mean(pts(:,1));
pts(:,2) = pts(:,2)-mean(pts(:,2));
pts(:,3) = pts(:,3)-mean(pts(:,3));

% plotting
pplot = 1;
if pplot == 1
  dt = 0.05;
  a = r(1);
  b = r(2);
  c = r(3);
  x0 = 0;
  y0 = 0;
  z0 = 0;
  u = 0:dt:2*pi;
  v = 0:dt:pi;
  
  % pre-allocate plotting vectors 
  xs = zeros(length(u)*length(v),1);
  ys = zeros(length(u)*length(v),1);
  zs = zeros(length(u)*length(v),1); % but depends only on v
  
  % populate plotting vectors
 count = 0;
  for i = 1:length(u)
    for j = 1:length(v)
      count = count + 1;
      xs(count) = a * cos(u(i)) * sin(v(j));
      ys(count) = b * sin(u(i)) * sin(v(j));
      zs(count) = c * cos(v(j));
    end
  end
  
  % scatter3(xs, ys, zs, 'r.')
  % axis equal
  
  newpts = pts;
  epts = [xs, ys, zs];
  scatterSurface(epts, 0.2);
  hold on;
  % add ellipse gridpoints if desired %%%%%
  % scatter3(pts(:,1), pts(:,2), pts(:,3), 'k^', 'MarkerFaceColor','k');
  
end
