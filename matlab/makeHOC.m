function hoc = makeHOC(amiraSegmentsFile, volumeFile, segLengthsFile)

% Creat a hoc file from Amira files
% Inputs (as strings '*.txt')
%       1. amiraSegmentsFile - taken from SpatialGraph.am
%       2. volumeFile - scaled to amiraSegments where each row is an x-y-z
%       triplet of coordinates of 1 volume unit (voxel); here taken as 1
%       3. segLengthsFile - taken from SpatialGraph.am
% Output
%       1. a hoc file (also saved as hoc.txt which should be renamed)

