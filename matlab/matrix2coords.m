function newpoints = matrix2coords(matrix, varargin)
% newpoints = matrix2coords(matrix, varargin[ (x,y,z), Switch 1/0 ]
% inputs: matrix (N x M x P), pixel width array, switch - swtiches x and y
% coordinates (as Amira seems to do), default is SWITCH ON (=1)
% convert index back to coordinate system
% 2nd argument is option pixel width array as [x,y,z] pixel widths; if none
% are given then 1 is assumed;

if nargin == 2
    pixWidth = varargin{1};
else
    pixWidth = [1 1 1];
end

if nargin == 3
    Switch = varargin{2};
else
    Switch = 1;
end

fprintf('All data loaded correctly.\n');

% pixWidth = 0.0002974; % mm
if Switch == 1
    fprintf('Reversing X- and Y-coordinates...');
    [yco, xco, zco] = ind2sub(size(matrix), find(matrix));

    xco = xco*pixWidth(1);
    yco = yco*pixWidth(2);
    zco = zco*pixWidth(3);
    fprintf(' done.\n');

else
    fprintf('Not reversring X- and Y-coordinates...');
    [xco, yco, zco] = ind2sub(size(matrix), find(matrix));

    xco = xco*pixWidth(1);
    yco = yco*pixWidth(2);
    zco = zco*pixWidth(3);
    fprintf(' done.\n');
end

%if pplot == 1

    % plot
    % scatter3(xco, yco, zco, 'k.')
%else
%end
fprintf('Creating new coordinate matrix...');
    newpoints = zeros(length(xco),3);
    newpoints(:,1) = xco;
    newpoints(:,2) = yco;
    newpoints(:,3) = zco;
    fprintf(' done. \n');
    
