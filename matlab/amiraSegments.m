function orderedSegs = amiraSegments(segmentLengthFile, skelpointFile)
% orderedSegs = amiraSegments(segmentLengthFile, skelpointFile)
% creates an Nx4 array where col 1 = seg number (same for all points in
% segment) and col 2:4 = xyz tuple of each point (includes endpoints)
% input:
%       1. segmentLengthFile ('*.txt') exactly as taken from amira spatial
%       graph file
%       2. skelpointFile ('*.txt') exactly as taken from amira spatial
%       graph file
% output:
%       1. ordered segments as Nx4 matrix


% importdata
segLengths = importdata(segmentLengthFile);
skelpoints = importdata(skelpointFile);

% allocate new matrix
[numPoints, col] = size(skelpoints);
orderedSegs = zeros(numPoints, col+1);
orderedSegs(:,2:4) = skelpoints;

ProgressBar('TotalPoints',numPoints);
% index col 1 for all points

points_so_far = 1;
for seg = 1:length(segLengths(:,1));
    
    orderedSegs(points_so_far:(points_so_far+segLengths(seg)-1),1) = seg;
    
    points_so_far = points_so_far + segLengths(seg);
    ProgressBar('TotalPoints');
end

% change save name in between
save new_orderedSegs_836 orderedSegs -ascii;
