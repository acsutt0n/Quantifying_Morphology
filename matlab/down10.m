function downsampledData = down10(dataFile)
% downsampledData = down10(dataFile)
% downsample by 10
fprintf('Loading %s ...', dataFile)
normData = importdata(dataFile);
fprintf(' done. \nDownsampling ...')
downsampledData = downsample(normData,10);
fprintf(' done.')