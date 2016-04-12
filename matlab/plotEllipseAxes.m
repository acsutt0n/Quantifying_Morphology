function plotEllipseAxes(pts)
% usage: plotEllipseAxes(pts)

% finds major and minor axes of ellipse and plots them with ellipse points

% get eigenvectors and eigenvalues for the pts
[coeff, ~, lat] = pca(pts);

% scale latent ("eigenvalues")
far = farthestPt(pts);
lat = lat/max(lat);
lat = lat*far;

% make eigenvectors (3 xyz triplets) and then move them
A1 = eigPoints(coeff, lat, 1);
A2 = eigPoints(coeff, lat, 2);
A3 = eigPoints(coeff, lat, 3);
ptmean = mean(pts);

for r = 1:3
  A1(:,r) = A1(:,r) + ptmean(r);
  A2(:,r) = A2(:,r) + ptmean(r);
  A3(:,r) = A3(:,r) + ptmean(r);
end
A = [A1; A2; A3];



% plot the points for the axes
hold on;
color = ['k','r','b'];

for c = 1:3 % m = AXIS
  plot3( A(c*3-2:c*3,1), A(c*3-2:c*3, 2), A(c*3-2:c*3,3), color(c));
  set( findobj(gca,'type','line'), 'LineWidth', 2);
end
  
  %plot3([-(coeff(c,1)*0.5*lat(c)), 0, (coeff(c,1)*0.5*lat(c))], ...
  %      [-(coeff(c,2)*0.5*lat(c)), 0, (coeff(c,2)*0.5*lat(c))], ...
  %      [-(coeff(c,3)*0.5*lat(c)), 0, (coeff(c,3)*0.5*lat(c))], color(c));


scatter3(pts(:,1), pts(:,2), pts(:,3), 'g.')






function far = farthestPt(pts)
far = 0;
for m = 1:length(pts)
  for n = 1:length(pts)
    dist = norm(pts(m)-pts(n));
    if dist > far
      far = dist;
    end
  end
end


function points = eigPoints(coeff, lat, v)
points = [ [-(coeff(v,1)*0.5*lat(v)), 0, (coeff(v,1)*0.5*lat(v))]; ...
         [-(coeff(v,2)*0.5*lat(v)), 0, (coeff(v,2)*0.5*lat(v))]; ...
         [-(coeff(v,3)*0.5*lat(v)), 0, (coeff(v,3)*0.5*lat(v))] ];
points = points';

% format: [ x1 y1 z1
%           x2 y2 z2
%           x3 y3 z3 ]








       