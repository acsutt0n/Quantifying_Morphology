function [xs, ys, xell, yell] = zslice(a, b, c, scaleOn, pplot)
% usage: zslice(a, b)
% 
% this plots a slice of an ellipsoid at a given z

if nargin > 3
  scaleOn = scaleOn;
else
  scaleOn = 0; % default - no scale
end

if nargin > 4
  pplot = pplot;
else
  pplot = 0; % default is no plot -- fxn is called in for loop
end

dt = 0.1;
u = [0:dt:2*pi];
v = [0:dt:pi];
scale = [0:0.1:1];
xs = zeros(length(u)*length(v)*length(scale),1);
ys = zeros(length(u)*length(v)*length(scale),1);
zs = zeros(length(u)*length(v)*length(scale),1);

if scaleOn == 1
  count = 0;
  for s = 1:length(scale)
    for i = 1:length(u)
      for j = 1:length(v)
        count = count+1;
        xs(count) = scale(s)*a * cos(u(i)) * sin(v(j));
        ys(count) = scale(s)*b * sin(u(i)) * sin(v(j));
        % zs(count) = c * cos(v(j));
      end
    end
  end
else
  count = 0;
  for i = 1:length(u)
    for j = 1:length(v)
      count = count+1;
      xs(count) = a * cos(u(i)) * sin(v(j));
      ys(count) = b * sin(u(i)) * sin(v(j));
      % zs(count) = c * cos(v(j));
    end
  end
end

% figure()
% plot(xs, ys, 'k-')
% axis equal

xell = zeros(length(u)*length(v)*length(scale),1);
yell = zeros(length(u)*length(v)*length(scale),1);
zell = zeros(length(u)*length(v)*length(scale),1);

if scaleOn == 1
  count = 0;
  for s = 1:length(scale)
    for i = 1:length(u)
      count = count+1;
      if u(i) > pi
        vell = u(i) - pi;
      else
        vell = u(i);
      end
        xell(count) = scale(s)*2*a * cos(vell) * sin(vell);
        yell(count) = scale(s)*2*b * sin(vell) * sin(vell);
    end
  end
else
  count = 0;
  for i = 1:length(u)
    count = count+1;
    if u(i) > pi
      vell = u(i) - pi;
    else
      vell = u(i);
    end
      xell(count) = 2*a * cos(vell) * sin(vell);
      yell(count) = 2*b * sin(vell) * sin(vell);
  end
end

% compensate for yell being shifted up by no sin > pi (no negative sin
% term)

R = unique([xell, yell],'rows');
yu = unique(yell);
yell = yell - mean(yu);

if pplot == 1
  figure()
  scatter(R(:,1), R(:,2), 'r.')
  alpha(0.5)
  axis equal
end

xell = R(:,1); yell = R(:,2);

%fprintf('Black: xmax %.3f, xmin %.3f, ymax %.3f, ymin %.3f \n', ...
%       max(xs), min(xs), max(ys), min(ys));
%fprintf('Red: xmax %.3f, xmin %.3f, ymax %.3f, ymin %.3f \n', ...
%       max(xell), min(xell), max(yell), min(yell));


end


