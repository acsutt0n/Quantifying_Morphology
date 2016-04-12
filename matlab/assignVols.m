function [segVols] = assignVols(orderedSegsFile, volCoordsFile)
% [segVols] = assignVols(orderedSegsFile, imgDirectory)
% assigns volumes to segments
% input:
%       1. orderedSegsFile as output from amiraSegments.m ('*.txt')
%       2. volCoords is a matlab matrix variable Nx3 produced after scaling the output from
%       matrix2coords.m (as determined by getBounds.m) with scaleCoords.m;
%       original output is from tiff_stack.m
% output:
%       1. segVols - numSegs x 1 array of volumes for each segment
%       2. ptVols - numPoints x 5 array where col 1 = segNum, col 2:4 = xyz
%       tuple, col 5 = volume of *segment* that point belongs to

% load data
fprintf('Importing %s ...', orderedSegsFile)
segPoints = importdata(orderedSegsFile);
fprintf(' done.\n')
segFirstCol = segPoints(:,1);
segPoints = segPoints(:,2:4);
fprintf('Importing %s ...', volCoordsFile)
volCoords = importdata(volCoordsFile);
fprintf(' done.\n')
% create matrix for areas
numSegs = max(segFirstCol);
segVols = zeros(numSegs, 1);
numVolCoords = length(volCoords(:,1));
segInd = zeros(numVolCoords,1);

MyBlock = ParallelBlock(); %#ok<NASGU>
ProgressBar('VolumeCoordinates',numVolCoords);

parfor h = 1:numVolCoords
    current_volpoint = volCoords(h,:);
    
    % find closest point & segment
    segInd(h) = ClosestPoint(current_volpoint, segPoints);
    %segVols(close_ind) = segVols(close_ind) + 1;
    
     %current_seg = segPoints(h,1);
     %segVols(current_seg) = segVols(current_seg) + 1;
     
     ProgressBar('VolumeCoordinates');
end

ProgressBar('AssignVols',numVolCoords)

for i = 1:numVolCoords; % for each volume coord
    
    segment = segFirstCol(segInd(i)); % find the segment it was positive for
    segVols(segment) = segVols(segment) +1;
    ProgressBar('AssignVols');
end
