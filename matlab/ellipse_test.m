
D = 3;

% 2-D ellipse
if D == 2
  a=5; % horizontal radius
  b=10; % vertical radius
  c=5; % z axis
  x0=0; % x0,y0 ellipse centre coordinates
  y0=0;
  t=-pi:0.01:pi;
  x=x0+a*cos(t);
  y=y0+b*sin(t);
  plot(x,y)
end

% 3-D ellipsoid
if D == 3
  dt = 0.1;
  a = 5;
  b = 10;
  c = 5;
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
  
  scatter3(xs, ys, zs)
  axis equal
end
