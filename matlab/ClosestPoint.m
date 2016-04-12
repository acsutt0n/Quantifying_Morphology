% function closestPoint
% given one point and a matrix of N-by-3 triplets,

function current1_index = ClosestPoint(current0, all_points, version)

%note: drop current1 => edit calling functions
if nargin < 3
   version = 1;
end

if version == 1
  dist = sum(bsxfun(@minus, current0, all_points).^2, 2);
  [~, current1_index] = min(dist);
elseif version == 2
  dist = sum((all_points - ones(size(all_points,1),1) * current0).^2, 2);
  [~, current1_index] = min(dist);
elseif version == 3
  dist = sum((repmat(current0, size(all_points,1), 1) - all_points).^2, 2);
  [~, current1_index] = min(dist);
else
  minDist = Inf;
  for i = 1:size(all_points,1)
    dist = sum((all_points(i,:) - current0).^2);
    if dist < minDist
      minDist = dist;
      current1_index = i;
    end
  end
end
