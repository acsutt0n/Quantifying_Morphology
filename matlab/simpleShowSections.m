function simpleShowSections(SurfPointsFile, SkelPointsFile, varargin)
% usage: simpleShowSections(SurfPointsFile, numPlots)
% This function simply plots a few of the sections from SurfPointsFile,
% which takes the format shown below (entered as a string).
%  Inputs:
%    SurfPointsFile - (str) from basicNeuronGeometry (N x 5)
%    SkelPointsFile - (str or matrix) N x 5 created by getSkeleton.py
%       (optional -- use 0 is no file desired)
%    numPlots - if int, then the number of 3-D sections plotted is the
%       number given (if left blank, 10 is used, so the first 10 sections
%       are plotted)
%       - if array, then those particular sections are plotted; but the
%       sections will be +1 incremented from the rawHoc / getSkeleton.py
%       files, but will match the sections from basicNeuronGeometry.m
%  Outputs:
%    graphs only
%
% newSurfPoints (from basicNeuronGeometry) is of format:
%  [  x    y     z    area    section ]


% data loader
if ischar(SurfPointsFile)
  fprintf('Loading data from %s ...', SurfPointsFile)
  surfPoints = importdata(SurfPointsFile);
  fprintf(' done.\n')
else
  surfPoints = SurfPointsFile;
end

if ischar(SkelPointsFile)
  fprintf('Loading data from %s ...', SkelPointsFile)
  skelpoints = importdata(SkelPointsFile);
  fprintf(' done.\n')
elseif SkelPointsFile == 0
  skelpoints = 0;
else
  skelpoints = SkelPointsFile;
end
disp(size(skelpoints))

% correct python->matlab section numbering scheme so no sections are 0
if skelpoints ~= 0
  Zs = find(skelpoints(:,1) == 0);
  if any(Zs)
    skelpoints(:,1) = skelpoints(:,1) + 1;
  end
end

% verify point ranges are identical (for sections)
[minsurf, maxsurf, minskel, maxskel] = getRanges(surfPoints, skelpoints);
diffmax = maxsurf - maxskel;
diffmin = minsurf - minskel;
if diffmax > 0
  skelpoints(:,1) = skelpoints(:,1) + 1;
elseif diffmax < 0
  surfPoints(:,5) = surfPoints(:,5) + 1;
end
[minsurf, maxsurf, minskel, maxskel] = getRanges(surfPoints, skelpoints);
if maxskel ~= minskel
  disp('range is off: max')
end
if minskel ~= minsurf
  disp('range is off: min')
end
  

% parse input arguments
if nargin < 3
  numPlots = 10;
elseif nargin == 3
  arg2 = varargin{1};
  if length(arg2) == 1
    numPlots = varargin{1};
  elseif length(arg2) > 1
    numPlots = 0;
    arrayPlots = varargin{1};
  end
else
  disp('cannot parse input argument; maybe third input argument is')
  disp(' not enclosed in []')
  numplots = -1;
end


%  ---------------------------------------------------
%                 PLOTTING FIRST FEW
%  ---------------------------------------------------

if numPlots > 0
  
  for p = 1:numPlots
    
    newInds = find(surfPoints(:,5) == p);
    newPoints = surfPoints(newInds, 1:3);
    figure(p);
    hold on;
    scatter3(newPoints(:,1), newPoints(:,2), newPoints(:,3),'k.');
    if length(skelpoints) > 1
      skelInds = find(skelpoints(:,1) == p);
      plotskel = skelpoints(skelInds, 2:4);
      scatter3(plotskel(:,1), plotskel(:,2), plotskel(:,3),'r*');
    end
    hold off;
  end
  

%  ---------------------------------------------------
%               PLOT SPECIFIED SECTIONS
%  ---------------------------------------------------
  
elseif numPlots == 0
  
  for p = 1:length(arrayPlots)
    
    newInds = find(surfPoints(:,5) == arrayPlots(p));
    newPoints = surfPoints(newInds, 1:3);
    figure();
    scatter3(newPoints(:,1), newPoints(:,2), newPoints(:,3),'k.');
    
  end
  
elseif numPlots == -1
  % do nothing -- warning was already returned
end
  
end










function [minsurf, maxsurf, minskel, maxskel] = ...
  getRanges(surfMatrix, skelmatrix)

% surf:
minsurf = min(surfMatrix(:,5));
maxsurf = max(surfMatrix(:,5));
% skel:
minskel = min(skelmatrix(:,1));
maxskel = max(skelmatrix(:,1));

end


















