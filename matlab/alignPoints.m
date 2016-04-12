function newpts = alignPoints(pts, do, scale)
%
% Flip or rotate points to align with those exported from the neuron
% function
% usage: newpts = alignPoints(pts, do)
%   input: N x 3 pts, do, scale
%   where do = 'rotate': changes the axis alignment
%              'translate': turns hour-glass into sphere
%              'move': shift in x-y-z to new center location
%              'scale': multiply all values by scalar
%   where scale = 'rotate': [x-y axis, x-z axis]
%                 'translate': 0 (no input required)
%                 'move': [new_x, new_y, new_z]
%                 'scale': scalar
%  output: newpts (N x 3)
%

