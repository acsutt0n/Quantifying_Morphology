function [avgDist,distList] = avgDist(skelpointsFile,segLengthsFile)
% returns the quantal (smallest incremental) distance between points within
% the same segment (useful for clustering)
% Input: 1 - orderedSkelpointsFile - *.txt as exported by Amira
%        2 - segLengthsFile - *.txt as exported by Amira
%         -> both inputs as strings
% Output:  - avgDist with histogram for each segment's avg distance

% load data
fprintf('Loading %s ...', skelpointsFile)
skelpoints = importdata(skelpointsFile);
fprintf(' done. \n Loading %s ...',segLengthsFile)
lengths = importdata(segLengthsFile);
fprintf(' done.\n')

numSegs = length(lengths);
distList = zeros(numSegs,1);
ptsofar = 0;

ProgressBar('AverageDistances',numSegs);

% find distance between points in a given segment
for seg = 1:numSegs
    distSums = 0;
    
    for pt = 1:lengths(seg) - 1
        distSums = distSums + point_dist2(skelpoints(ptsofar+pt,:), ...
            skelpoints(ptsofar+pt+1,:));

            ptsofar = ptsofar + 1;

        segdist = distSums / (lengths(seg) - 1);
    end
    
    distList(seg) = segdist;
    ProgressBar('AverageDistances');
end

avgDist = mean(distList);
figure;
bins = [0.1:0.1:5];
hist(distList,bins);

end
