function [newSurfPoints, sectionData] = ...
          basicNeuronGeometry(skelFile, surfFile, volFile)
% usage: [secVolumes, surfPoints, secSurfs] = ...
%         basicNeuronGeometry(skelFile, surfFile, volFile)
% Output: volPoints - a list of volume points (millions) and which section
%                     they are assigned to (not a good idea to save or even
%                     return in most cases)
%                     VOLPOINTS output removed for time-saving
%         surfPoints - a list of surface points (lots) and which section
%                     they are assigned to; this is considerably smaller
%                     than volumes; can be saved if desired (N x 4)
%         sectionData - each section has a surface area, raw(int) and
%                     calculated (float) volume with it (N x 4)
% Input: skelFile - (str) a N x 5 skeleton file created from getSkeleton.py 
%                   on an Amira/Imaris Hoc file (sec, x, y, z, rad)
%        surfFile - (str) a N x 4 file created from plotCents.m where each
%                   row is a xyz centroid of a surface point; col 4 is area
%        volFile -  (str) a N x 3 file created from matrix2coords.m where
%                   each row is a xyz tuple of a voxel 


% data loader
fprintf(' Loading %s ...', skelFile)
skelmatrix = importdata(skelFile);
fprintf(' done. \n Loading %s ...', surfFile)
surfpoints = importdata(surfFile);
fprintf(' done. \n Loading %s ...', volFile)
volpoints = importdata(volFile);
fprintf(' done.\n')


% initializations

% volumeOn ???????
volumeOn = 0;

% skels
skelpoints = skelmatrix(:,2:4);
skelSections = skelmatrix(:,1);
skelRadius = skelmatrix(:,5);
numskelpoints = length(skelpoints(:,1));

  % check for section IDs = 0, if they exist increment all sections +1
  if any(skelSections == 0)
    skelSections = skelSections + 1;
  end

% surface & volume
numSurfPoints = length(surfpoints(:,1));
numvolpoints = length(volpoints(:,1));
newSurfPoints = [surfpoints, zeros(length(surfpoints(:,1)),1)];
% section
numSections = max(skelSections);
sectionVols = zeros(numSections, 1);
sectionData = zeros(numSections, 4);
sectionData(:,1) = [1:numSections]';

% VOXEL SIZE:
voxSize = 0.09 * 0.09 * 0.5; % um^3



%    -----------------------------------------------------------
%                             ASSIGN AREAS
%    -----------------------------------------------------------

fprintf('\nAssigning areas ...')

MyBlock = ParallelBlock(); %#ok<NASGU>
ProgressBar('SurfaceAreas',numSurfPoints);

% disp(size(skelpoints))
% disp(size(surfSegInd))

parfor h = 1:numSurfPoints             % for each surface point, write the
                                       % new surface area to col 2 of 
                                       % sectionData
                                       % and add the section # to
                                       % newSurfPoints
  current_point = surfpoints(h,1:3);
  
  % returns the index of the closest skelpoint (skelpoint(i,:))
  %surfSegInd(h) = AssignPoint(current_point, skelpoints(:,2:5));
  skelPointIndex = AssignPoint(current_point, skelpoints, skelRadius);
  % returns the section of the closes skelpoint
  segment = skelSections(skelPointIndex);
  newSurfPoints(h,5) = segment;
  
  ProgressBar('SurfaceAreas');
end

ProgressBar('AssignAreas',numSurfPoints)

for i = 1:numSurfPoints; % for each volume coord
  
  % update each sec for new area added
  %  newSurfPoints is:
  %  [    x      y     z     area     section    ]
  %  sectionData is:
  %  [   section   area    raw_volume    volume  ]
  
  sectionData(newSurfPoints(i,5),2) = sectionData(newSurfPoints(i,5),2) + ...
    newSurfPoints(i,4);
    
    ProgressBar('AssignAreas');
end

fprintf(' done.\nAssigning volumes ...')

% turn off volumes?
if volumeOn == 1
%    ------------------------------------------------------------
%                             ASSIGN VOLUMES
%    ------------------------------------------------------------

%MyBlock = ParallelBlock(); %#ok<NASGU>
ProgressBar('VolumePoints',numvolpoints);


parfor h = 1:numvolpoints
    current_volpoint = volpoints(h,:);
    
    % returns index of closest skelpoints
    volPointIndex = AssignPoint(current_volpoint, skelpoints, skelRadius);
    segment = skelSections(volPointIndex);
    sectionVols(h) = segment;
     
     ProgressBar('VolumePoints');
end

ProgressBar('AssignVolumes',numSections)

for i = 1:Sections;         % for each voxel, write raw to col 3 and
                            % calculated volume to col 4 of sectionData
    
    sectionData(sectionVols(i),3) = numel(find(sectionVols == i));
    sectionData(sectionVols(i),4) = sectionData(sectionVols(i),3) * voxSize;

    ProgressBar('AssignVolumes');
end

fprintf(' done.\n')


elseif volumeOn == 0
  disp('volume operation skipped')
end

end



% modified closest point function with hoc radius limits
% assumes operating under version 1
function current1_index = AssignPoint(current0, all_points, ...
  points_radius, version)
  % all_points = [x,y,z,rad]
  %note: drop current1 => edit calling functions
if nargin < 4
   version = 1;
end

if version == 1
  all_points_copy = all_points;
  % all_points_dist = zeros(length(all_points(:,1)),4);
  dist = sum(bsxfun(@minus, current0, all_points_copy).^2, 2);
  % all_points_dist = [all_points, dist];
    % all_points_dist = [x,y,z,rad,dist] sorted by dist
  [~, current1_index] = min(dist); %[close_dist, current1_index] = min(dist)
  % disp(close_dist)
  radius = points_radius(current1_index);
  
  % disp(size(dist))
  % disp(size(close_dist))
  
%   while close_dist > 1.5 * radius
%     % current1_index
%     dist(current1_index) = [];                 % delete invalid distance
%     all_points_copy(current1_index,:) = [];    % delete invalid point
% 
%     [close_dist, current1_index] = min(dist);  % find new close point
%     radius = points_radius(current1_index);    % find & check new radius
%     
%   end
  % disp(size(all_points_copy(current1_index,:)))
  % disp(size(all_points))
    [~, current1_index] = ismember(all_points_copy(current1_index,:), ...
      all_points, 'rows');
   % current1_index

  
% Other Versions:
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
end
