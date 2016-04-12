function [newSurfPoints, sectionData, troublemakers] = ...
          basicNeuronGeometryRadius(skelFile, surfFile, volFile)
% usage: [secVolumes, surfPoints, secSurfs] = ...
%         basicNeuronGeometry(skelFile, surfFile, volFile)
% Output: volPoints - a list of volume points (millions) and which section
%                     they are assigned to (not a good idea to save or even
%                     return in most cases)
%                     -> VOLPOINTS output removed for time-saving
%         surfPoints - a list of surface points (lots) and which section
%                     they are assigned to; this is considerably smaller
%                     than volumes; can be saved if desired (N x 5)
%                     [ x  y  z  area  section ]
%         sectionData - each section has a surface area, raw(int) and
%                     calculated (float) volume with it (N x 4)
%                     [ section  x  y  z ]
%         flags       - flags indexed surfpoints with distances > radius
%                     allowance
% Input: skelFile - (str) a N x 5 skeleton file created from getSkeleton.py 
%                   on an Amira/Imaris Hoc file (sec, x, y, z, rad)
%        surfFile - (str) a N x 4 file created from plotCents.m where each
%                   row is a xyz centroid of a surface point; col 4 is area
%        volFile -  (str) a N x 3 file created from matrix2coords.m where
%                   each row is a xyz tuple of a voxel 

volumeOn = 0;
trouble = 0;

if ischar(skelFile)
  % data loader
  fprintf(' Loading %s ...', skelFile)
  skelmatrix = importdata(skelFile);
  fprintf(' done. \n Loading %s ...', surfFile)
  surfpoints = importdata(surfFile);
  fprintf(' done. \n ')
  
  if volumeOn == 1
    fprintf('Loading %s ...', volFile)
    volpoints = importdata(volFile);
    fprintf(' done.\n')
    numvolpoints = length(volpoints(:,1));
  end
  
else
  skelmatrix = skelFile;
  surfpoints = surfFile;
  if volumeOn == 1
    volpoints = volFile;
    numvolpoints = length(volpoints(:,1));
  end
  fprintf('Data loaded.\n')
end

% initializations

% volumeOn ???????   troublemakers on ???????


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

newSurfPoints = zeros(numSurfPoints,5);
newSurfPoints = [surfpoints, zeros(length(surfpoints(:,1)),1)];
troublemakers = zeros(numSurfPoints,1);
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
  if length(skelpoints) == length(skelRadius)
    % do nothing
  else
    fprintf('points and radius not same size - trial %i',h)
    print('points is size: ')
    size(skelpoints)
    print('radius is size: ')
    size(skelRadius)
  end
  
%  [skelPointIndex, flag, ~] = AssignPoint ...
%    (current_point, skelpoints, skelRadius);
  [skelPointIndex, flag, ~] = AssignPoint ...
    (current_point, skelpoints, 0, 2);
  
  if flag == 1              % flag sections that are outside radius tolerance
   troublemakers(h) = h;
 %  troublemakers(h,2) = Dist;
   segment = skelSections(skelPointIndex);
   newSurfPoints(h,5) = segment; % assign bogus section to flag troublemaker
   ProgressBar('SurfaceAreas');
   
  else
    % returns the section of the closes skelpoint
    segment = skelSections(skelPointIndex);
    newSurfPoints(h,5) = segment;
  
    ProgressBar('SurfaceAreas');
  end
  
end

save newSurfPoints.txt newSurfPoints -ascii;
print('Saved surf point progress.')


% %%%%%%%%%%%%%%%%%%%%%%%%% deal with troublemakers %%%%%%%%%%%%%%%%%%%

% deal with troublemakers?
if trouble == 1
  fprintf('%d \% of points were troublemakers.\n',...
    numel(troublemakers)/numSurfPoints)

  % for troublemakers, find the nearest surfpoint that is assigned and assign
  % it to the same section

  % only keep previously-assigned surfpoints
  assignedInds = setdiff(1:numSurfPoints, troublemakers);

  % assignedInds = find(newSurfPoints(:,5)>0);
  assignedPoints = newSurfPoints(assignedInds, 1:3);

  % troubleInds = find(troublemakers);
  numTroublemakers = numel(troublemakers);
  troublePoints = surfpoints(troublemakers(troublemakers>0),1:3);
  assignedRads = zeros(length(assignedInds),1);
  % troublemakers(troubleInds, :);

  for i = 1:length(assignedInds(:,1))
    try
      assignedRads(i) = skelRadius(assignedInds(i));
    catch
      fprintf('Exception in skelRadius(assignedInds( %d ))\n', i)
      break
    end
  end


  ProgressBar('Troublemakers',numTroublemakers);

  fprintf('Num of assignedPoints is %i', length(assignedPoints))
  fprintf('Num of troublemakers is %i', length(troublePoints))

  for t = 1:numTroublemakers

    % only send the troublePoints and the radius values associated with
    % previously assigned points

    currentTrouble = troublePoints(t,:);

    % find closest assigned point  end
    [assignedPointIndex, ~, DDist] = AssignPoint ...
      (currentTrouble, assignedPoints);

    if DDist < troublemakers(t,2)  
      % find its associated section
      [~, assignedSurfInd] = ismember(assignedPoints(assignedPointIndex,1:3), ...
        newSurfPoints(:,1:3),'rows');
      [~, troubleSurfInd] = ismember(troublePoints(t,1:3), ...
        newSurfPoints(:,1:3),'rows');
      newSurfPoints(troubleSurfInd,5) = ...
        newSurfPoints(assignedSurfInd,5);
    end

    assignedPoints(end+1,:) = currentTrouble;
    % add new point to assigned points?
    ProgressBar('Troublemakers');

  end

  % make sure all points are assigned
  if any(find(newSurfPoints(:,5)<0))
    disp('Error - could not assign all troublemakers')
  else
    disp('All troublemakers assigned.')
  end
end % end of troublemaker stuff

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
fixedsurfpoints

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



%%%%%%%%%%%%%%%%%%%%%%%%%%% AssignPoint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modified closest point function with hoc radius limits
% assumes operating under version 1
function [current1_index, flag, close_dist] = AssignPoint(current0, all_points, ...
  points_radius, version)

  % all_points = [x,y,z,rad]
  %note: drop current1 => edit calling functions
if nargin < 3
  version = 2;
end
if nargin < 4
   version = 1;
end

if version == 1
 % start = clock;                                    % CLOCK
 % TimeOut = 0.1; % secs

  dist = sum(bsxfun(@minus, current0, all_points).^2, 2);
  distrad = dist ./ points_radius;

  [close_dist, current1_index] = min(distrad); 
  % if the point is > 2 times the ratio of the distance to the radius, flag
  if close_dist > 2
    flag = 1;
  else
    flag = 0;
  end

elseif version == 2

  dist = sum(bsxfun(@minus, current0, all_points).^2, 2);
  [close_dist, current1_index] = min(dist);
  flag = 0;
  
 %   if (etime(clock, start) > TimeOut)                % CLOCK
 %     flag = 1;
 %     return
 %   end


 end
%   
%   if close_dist < close_dist2
%     current1_index = current2_index;
%   else
%     [~, current1_index] = ismember(all_points_copy(current1_index,:), ...
%       all_points, 'rows');
%   end

  
end

