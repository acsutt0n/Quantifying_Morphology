% function scaled_stl = scaleSTL(unscaled_stl, scaled_volume, varargin)
% scale STL coordinates
% this function operates on the assumption that the max/min of the volume
% coordinates will be almost equal to the STL max/min and fits them
% accordingly
% input: 
%   1. the format of the STL is an N x 3 x 3 with N triangles from READ_stl
%   2. scaled volume (taken from matrix2coords)
%   there should be no need to switch the inputs in this case

if nargin == 3
    Switch = varargin{1};
else
    Swtich = 1; % default 
end

% get bounds for all, compare

volBounds = getBounds(scaled_volume);
stlBounds = getBounds(unscaled_stl(:,:,1));

fprintf('Bounds for %d are: \n', scaled_volume);
fprintf('Bounds for %d are: \n', unscaled_stl);


%% INCOMPLETE -
% abandoned because only the volume was off, the stl is to scale (and
% proper x-y-z)

    
