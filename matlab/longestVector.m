function [v, dist] = longestVector(pts)
% usage: v = longestVector(pts)
%
% returns a unit vector in the direction of the line connecting the
% two farthest points from pts

ProgressBar('Checking for longest vector...', length(pts(:,1)));

dist = 0;
for m = 1:length(pts(:,1))
  for n = 1:length(pts(:,1))
    d = norm(pts(m,:) - pts(n,:));
    if d > dist
      pt0 = pts(m,:);
      pt1 = pts(n,:);
      dist = d;
    end
  end
  ProgressBar('Checking for longest vector...');
end

v = zeros(3,1);
for i = 1:3
  v(i) = pt0(i)-pt1(i);
end
% disp(dist)
end