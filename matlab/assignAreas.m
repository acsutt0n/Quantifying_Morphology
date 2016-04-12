% assign areas to segments contained in an array of cells with Nx3

function [segment_areas] = assignAreas (surface, segments)
% input     - surface: matrix surface (Nx4) from plotCents.m
%           xyz centroid (:,1:3); area (:,4)
%           - segments: direct output from orderSegments.m; matrix is a
%           Nx4 with N total skelpoints (no endpoints) in (:,2:4) and
%           segment number (from connections ascii Amira skeleton file)
%           in (:,1)

% this function will assign an area to each segment
% based on its centroid's proximity to the closest skelpoint in a given
% segment;
% output - Nx1 vector of N segments whose index corresponds to the segment
% # from connections data in ascii Amira skeleton file and the input to
% this function

ProgressBar('SurfacePoints',length(surface(:,1)));

% find total number of segments; create empty vector
numSegs = max(segments(:,1));
segment_areas = zeros(numSegs,1);

parfor n = 1:length(surface(:,1));
    
    % find closest point
    [current_point, current_seg] = getClosestPoint(surface(n,1:3), ...
        segments ); % go through all segment points
    
    % add area to segment (point is unnecessary)
    segment_areas(current_seg) = segment_areas(current_seg) + ...
        surface(n,4);
    
    ProgressBar('SurfacePoints');
    
end

end

% this sub-function is modified from closestPoint but modified to return
% the value from col 1 of segments (segNum) rather than the index
    function [current1, current1_seg] = getClosestPoint(current0, all_points)

% some huge distance
dist2 = 1000;

for i = 1:length(all_points(:,1))
    current_temp = all_points(i,2:4);
    dist1 = point_dist2(current0(1),current0(2),current0(3), ...
        current_temp(1),current_temp(2),current_temp(3));
    
    if dist1 < dist2 & dist1 ~= 0
        current1 = current_temp;
        current1_seg = all_points(i,1);
        dist2 = dist1;
    else
    end
    
end
    end
    