function majoraxis = findMajorAxis(somaCoordsFile)
% usage: majoraxis = findMajorAxis(somaCoordsFile)
% This function finds the major (longest) axis through the soma to be used
% as the "skeleton" segment for the soma in a .hoc file
% Inputs: 'somaCoordsFile.txt' as a string - this is an N x 3 matrix of
% the coordinates of the soma surface centroids
% Outputs: majoraxis - an M x 3 matrix of (x,y,x) coordinates of the
% segment

% This function is best if run after plotCents.m (to use surface points to
% find major axis)

fprintf('Loading data file %s ...', somaCoordsFile)
coords = importdata(somaCoordsFile);
fprintf(' done.\n')

numCoords = length(coords(:,1));

dist0 = 

for coord_ind = 1:numCoords
    