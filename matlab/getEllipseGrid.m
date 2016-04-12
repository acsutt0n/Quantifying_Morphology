function gridpoints = getEllipseGrid(ax, pplot)
% usage: gridpoints = getEllipseGrid(epts)
%
% input: ax is a 3x1 of the magnitudes for each of the maj/min axes
% output: xyz-triplets of the entire grid of the ellipse volume
%

% [~, ~, xell, yell] = zslice(ax(1), ax(2), ax(3));

if nargin > 1
  pplot = pplot;
else
  pplot = 0;
end

v = [0:0.01:pi]; % make sure this matches (in some way) ellipsoid_fit.m dt
zrange = ax(3) .* cos(v);

ax = ax;% * 2; % expand for fun
dax=[ax(1), ax(2), ax(3)*2];
ax=dax;

zs = linspace(ax(3), 0, 10); % downsample here
zs = zs - mean(zs);

grid = zeros(3,1);

for z = 1:length(zs)
  s = zs(z)/max(zs); % scale to max z
  % get scaled ellipse points
  [~, ~, xell, yell] = zslice(s*ax(1), s*ax(2), s*ax(3)); 
  R = [xell, yell];
  R = unique(R, 'rows');
  %R = sortrows(R, 2);
  R_pos_x = find(R(:,2) > 0);
  R = R(R_pos_x,:);
  R = flip(R,1); % not sure why flipping over end here...
  
  for n = 1:length(R(:,1))/2
    % dx = 10*abs(mean(diff(R(:,1))));
    % samples = ceil(abs(R(n,1)-R(end-n+1,1))/dx);
    % xxx = linspace(R(n,1), R(end-n+1,1), samples);
    xcurr = R(n,1);
    dy = 5*abs(mean(diff(R(:,2)))); % downsample here
    samples = ceil(abs(R(n,2)-R(end-n+1,2))/dy);
    yyy = linspace(R(n,2), R(end-n+1,2), samples);
    
    % subtract mean to re-center
    yyy = yyy - mean(yyy);
    
    for j=1:length(yyy)
      %if length(xxx) == length(yyy)
        grid(end+1, 1) = xcurr; % xxx(j);
        grid(end, 2) = yyy(j);
        grid(end, 3) = zs(end-z+1);
        %disp(size(grid));
      %end
    end
    
  end
  
end

grid = grid(1:end/2,:);
grid(:,3) = grid(:,3)+abs(min(grid(:,3))); % translate z so all positive
grid = [grid; grid];
grid(1:end/2, 3) = -grid(1:end/2, 3);
grid = [grid;grid];
grid(end/2:end,1) = -grid(end/2:end,1);
gridpoints = grid;

% adjust grid manually
grid(:,3) = grid(:,3) * 2;

if pplot == 1
  scatter3(grid(:,1), grid(:,2), grid(:,3), 'k.');
  xlabel('X Axis')
  ylabel('Y Axis')
  zlabel('Z Axis')
  axis equal;
end

end


