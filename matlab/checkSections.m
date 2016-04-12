function badSections = checkSections(surfPointsFile, ...
                                     skelPointsFile, varargin)
% usage: badSections = checkSections(surfPointsFile, skelPointsFile, 
%                                    specSections, pplot)
%
% Input:
%    surfPointsFile - (str) as output newSurfPoints from
%      basicNeuronGeometry.m
%      [   x     y     z     area     section ]
%    skelPointsFile - (str) input from getSkeleton.py
%      [ section  x    y    z    radius ]
%    specSections - (optional, array) specified sections for the function
%      to check; if blank or 0 the function will check them all
%    pplot - (optional, 1/0) if 1 the function will plot the bad sections
%      using simpleShowSections; if 0 or blank it will not plot (but the
%      output badSections can be passed directly to simpleShowSections
%      along with surfPointsFile to easily plot the sections)
% Output:
%    badSections - simple array of the sections that were found to be bad
%      or questionable (flagged)
% 


% data loader
fprintf('Loading %s ...', surfPointsFile)
surfMatrix = importdata(surfPointsFile);
fprintf(' done. Loading %s ...', skelPointsFile)
skelMatrix = importdata(skelPointsFile);
fprintf(' done.')

% parse inputs
if nargin < 3
  specSections = 0;
  pplot = 0;
elseif nargin < 4
  specSections = varargin{1};
  pplot = 0;
elseif nargin == 4
  specSections = varargin{1};
  pplot = varargin{2};
else
  disp('bad input arguments')
end

if pplot ~= 1 || pplot ~= 0
  disp('bad pplot (3rd) input; must be 1 or 0');
end

% initializations
% numSections = max(surfMatrix(:,5));
surfSections = surfMatrix(:,5);
surfPointsOnly = surfMatrix(:,1:3);

% check input and section names
skelSections = skelMatrix(:,1);  
skelPoints = skelMatrix(:,2:4);
if any(skelSections == 0)                % ensure skel & surf sections
  skelSections = skelSections + 1;       % match up
end
if max(skelSections) == max(surfSections)
  numSections = max(skelSections);
else
  disp('Warning: unequal number of sections in surf and skel points')
  numSections = min([max(skelSections), max(surfSections)]);
end






% -------------------- CHECK ALL SECTIONS --------------------

if specSections == 0
  
  parfor sec = 1:numSections
    
    tempInds = surfSections == sec;
    tempSurfPoints = surfPointsOnly(tempInds,:);
    tempInds = skelSections == sec;
    tempSkelPoints = skelPoints(tempKinds,:);
    for i = 1:length(tempInds)
      
    
    
    
    
  end
  
end



























