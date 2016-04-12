function volpoints = scaleMatrix(matrix, varargin)
% usage: volpoints = scaleMatrix(matrix, [voxel size], skelpointsFile)
%
% Inputs: matrix - output from tiff_stack.m, a M x N x Z matrix 
%                  computed from each M x N of Z images
%         voxel size (optional) - a 1 x 3 vector [x y z] of each 
%                  dimension's voxel size. If none is provided, [1x1x1]
%                  is used and the actual volume can be scaled after
%         skelpointsFile (optional) - given by getSkeleton.py, this .txt file is
%                  N x 5 of (segNum, x, y, z, hocRadius). Here, the
%                  x,y,z is used to scale the matrix coordinates. If
%                  none is given the matrix is not scaled and it just
%                  returned
% Output: volpoints - N x 3 of coordinates of the voxels in the neuron
%

% data loader
if nargin == 2
  pixWidth = varargin{1};
else
  pixWidth = [1 1 1];
end

if nargin == 3
  skelFile = varargin{2};
  if ischar(skelFile) == 0
    fprintf('skeleton points file should be a string.\n')
  else
    fprintf('Loading $s ... ', skelFile)
    skelpoints = importdata(skelFile);
    fprintf(' done.\n')
  end
else
  skelFile = 0;
end


% matrix to coordinate conversion
fprintf('Extracting indices from matrix ...')
[xco, yco, zco] = ind2sub(size(matrix), find(matrix));
fprintf(' done.\n')

xco = xco*pixWidth(1);
yco = yco*pixWidth(2);
zco = zco*pixWidth(3);

fprintf('Running correctDim ... ')
tempVols = [xco, yco, zco];
[volpoints, skelpoints2, ~] = correctDim(tempVols, skelpoints);
fprintf(' done.\n')

if prod(skelpoints2==skelpoints) < 1
  fprintf('Warning: new and old skelpoints appear different')
end

end

