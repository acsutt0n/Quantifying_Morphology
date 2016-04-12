function [segAreas] = assignAreasParfor(orderedSegsFile, areasFile)
% [segVols] = assignVols(orderedSegsFile, imgDirectory)
% assigns volumes to segments
% input:
%       1. orderedSegsFile as output from amiraSegments.m ('*.txt')
%       2. areasFile is a .txt created by plotCents.m where cols 1:3 are
%       x-y-z centroid tuple and col 4 is the area of that midpoint
% output:
%       1. segAreas - numSegs x 1 array of areas for each segment
%       - currently this is only output
%       2. ptVols - numPoints x 5 array where col 1 = segNum, col 2:4 = xyz
%       tuple, col 5 = volume of *segment* that point belongs to

% load data
fprintf('Importing %s ...', orderedSegsFile)
segPoints = importdata(orderedSegsFile);
fprintf(' done.\n')
segFirstCol = segPoints(:,1);
segPoints = segPoints(:,2:4);
fprintf('Importing %s ...', areasFile)
centroids = importdata(areasFile);
fprintf(' done.\n')

% create matrix for areas
numSegs = max(segFirstCol);
segAreas = zeros(numSegs, 1);
numCentroids = length(centroids(:,1));
segInd = zeros(numCentroids,1);

MyBlock = ParallelBlock(); %#ok<NASGU>
ProgressBar('AreaCentroids',numCentroids);

parfor h = 1:numCentroids
    current_point = centroids(h,1:3);
    
    % find closest point & segment
    segInd(h) = ClosestPoint(current_point, segPoints);
    %segVols(close_ind) = segVols(close_ind) + 1;
    
     %current_seg = segPoints(h,1);
     %segVols(current_seg) = segVols(current_seg) + 1;
     
     ProgressBar('AreaCentroids');
end

ProgressBar('AssignAreas',numCentroids)

for i = 1:numCentroids; % for each volume coord
    
    segment = segFirstCol(segInd(i)); % find the segment it was positive for
    segAreas(segment) = segAreas(segment) +1;
    ProgressBar('AssignAreas');
end
