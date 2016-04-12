
function [SWC, map] = orderSWC(swcFile)



if ischar(swcFile)
  swc = importdata(swcFile);
else
  swc = swcFile;
end


numSegs = length(swc(:,1));
remainingSegs = 1:numSegs;

% map = (N x 2) matrix of (current node, old node) as a map
map = zeros(numSegs, 2);
map(:,1) = (1:numSegs)';
fixSWC = zeros(size(swc));
root = find(swc(:,7)==-1);
fixSWC(1,:) = swc(root,:);
fixSWC(1,1) = 1; % set node ID of root seg to 1
map(1,2) = swc(root,1);
swc(root,:) = [];


% MyBlock = ParallelBlock(); %#ok<NASGU>
ProgressBar('OrderingSegments',numSegs);

while any(swc(:,1))
  
  
  for q = 1:length(swc(:,1))
    % if an unknown seg's parent have already been mapped, then map that
    % seg
    
    [qYN, qIndex] = ismember(swc(q,7), map(:,2));
    % qYN - is q in map? qIndex - where?
    
    if sum(qYN) > 0
      [row, ID] = findRow(fixSWC);  % find the next row & ID to add to fixSWC
      fixSWC(row,:) = swc(q,:);
      fixSWC(row,1) = ID;                         % set current ID
      fixSWC(row,7) = map(qIndex, 1);             % set parent ID
      map(ID,2) = swc(q,1);
      swc(q,:) = [];                              % remove that segment
      ProgressBar('OrderingSegments')
    end
    
  end
    
end


SWC = fixSWC;

end



function [row, ID] = findRow(swcmatrix)
current_rows = find(swcmatrix(:,1)>0);
[m,n] = size(swcmatrix);
if m < current_rows
  disp('pre-allocated matrix is too small')
else
  row = max(current_rows) + 1;
  ID = max(swcmatrix(:,1)) + 1;
end
end
  






