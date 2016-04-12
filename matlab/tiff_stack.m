function [data, matrix] = tiff_stack (varargin)
% load up split files for volume matrix
% only works with tiffs (but can be modified as needed)
% fileNames must be a string; easiest to run in current directory

if nargin == 1;
  tifFiles = dir(varargin{1});
else
  tifFiles = dir('*.tif');
end

numfiles = length(tifFiles);
data = cell(1,numfiles);
%names = cell(1,numfiles);

ProgressBar('ImportImages',numfiles+1);

for i=1:numfiles
  data{i} = imread(char(tifFiles(i).name));
  %[m,n] = size(data{i})
  ProgressBar('ImportImages');      
end

% print('Data is size:');
% [rows, cols] = size(data);
save data ;
fprintf('Data saved as matlab cell entries.\n')
    
matrix = cat(3,data{1:numfiles});
matrix = matrix ~= 0; % turn all positive integers into 1

end