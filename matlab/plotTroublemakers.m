function sectionPoints = plotTroublemakers(surfPoints, troubleMakers, varargin)

% usage: sectionPoints = plotTroublemakers(troublemakers, surfPoints, plot(op))
%{
Output - sectionPoints - a matrix of (x, y, z, secNum, troublemaker(1/0))
         that keeps track of a few sections and their troublemakers
       - a few representative sections are also plotted by default
Input  - surfPoints - a matrix or string for a filename with the surfpoints
         from basicNeuronGeometryRadius; this is an N x 5 matrix or file
         (x,y,z,area, segNum)
       - troublemakers - a vector of indices of troublemakers in surfPoints
         these will be plotted a different color to show what they look
         like
       - skelPoints (optional) - a string or matrix to load the skelpoints
         if want to specify plot arg but not give a skelPoint file, can
         simply put 0 here
       - plot (optional) - defaults to 10 where it plots 10 different
         segments highlighting the troublemakers; if a vector, only the
         enumerated segments will be plotted; if a scalar, plots that many
         segments (from the beginning)

%}


% Data loader / arg parser
if ishcar(surfPoints) == 1
  fprintf('Loading surfpoints from %s ...', surfPoints)
  surfpoints = importdata(surfPoints);
  fprintf(' ... done.\n')
else
  surfpoints = surfPoints;
  fprintf('Surfpoints loaded from workspace.')
end

if ishcar(troubleMakers) == 1
  fprintf('Loading troubleMakers from %s ...', troubleMakers)
  troublemakers = importdata(troubleMakers);
  fprintf(' ... done.\n')
else
  troublemakers = troubleMakers;
  fprintf('Troublemakers loaded from workspace.')
end

if nargin > 2
  skel = varargin{1};
  if ischar(skel) == 1
    fprintf('Loading skelPoints data from %s ...', skelPoints)
    skelpoints = importdata(skelPoints);
    skel = 1;
    fprint(' ... done.\n')
  elseif length(skel) > 1 % if not a string or a scalar
    skelpoints = skel;
    skel = 1;
  end
end

if nargin > 3
  pplot = varargin{2};
  print('Plotting user-defined segments.\n')
else
  pplot = [1:10];
  print('Plotting default of first 10 segments.\n')
end
if length(pplot) > 1
  vplot = 1;
  print('Multiple plots is on.\n')
else
  vplot = 0;
  print('Plotting only 1 graph.\n')
end


%%%%%%% split data into assigned points and troublepoints %%%%%%%%

allInds = 1:length(surfpoints(:,1));
assignInds = setdiff(allInds, troublemakers);
assignPoints = surfpoints(assignInds, 1:5);
troublePoints = surfpoints(troublemakers, 1:5);

if vplot == 0
  
  
end








end


%%% plotting function
function multiScatter(groupA, groupB)
% parse inputs
[mA,nA] = size(groupA);
[mB,nB] = size(groupB);
if nA ~= nB
  print('cols of groupA must be equal to cols of groupB \n')
end
if nA ~= 3
  print('Expected an object with 3 rows \n')
end
% condition inputs
Ax = groupA(:,1); Ay = groupA(:,2); Az = groupA(:,3);
Bx = groupB(:,1); By = groupB(:,2); Bz = groupB(:,3);
% plot
figure();
hold on;
scatter3d(Ax, Ay, Az, 'r.')
scatter3d(Bx, By, Bz, 'b.')




