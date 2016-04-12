function [fixA, fixB, fix] = correctDim(groupA, groupB)
% usage: [fixA, fixB, fix] = correctDim(groupA, groupB)
%
% This function receives two M x N matrices and, if they have the same
%   N (columns) dimensions, makes sure their ranges are similar (within
%   a given tolerance) in each dimension. It infers that col1 = x, 
%   col2 = y, col3 = z, and col4 = t.
%
% Input: groupA , groupB - M x N matrices with the same N dimensions
%
% Output: fixA, fixB - returns the new dimensions, may be swtiched 
%                      but the code for 'fix' tells you which
%                      dimensions were switched
%
%                      As a default, only A is changed; B is 'reference'
%
%         fix - which dimensions were switched?
%               Codes:
%                      'XY' - X and Y were switched
%                      'XZ' - ...
%                      'YZ' - ...
%

fix=0;
D=0;
while D == 0
  try
    [mA, nA] = size(groupA);
    [mB, nB] = size(groupB);
    D=2;
  catch d2err
    fprintf('2-D error. Trying 3-D...\n')
  end

  try
    [mA, nA, zA] = size(groupA);
    [mB, nB, zB] = size(groupB);
    D=3;
  catch d3err
    fprintf('3-D error. Trying 4-D...\n')
  end

  try
    [mA, nA, zA, tA] = size(groupA);
    [mB, nB, zB, tB] = size(groupB);
    D=4;
  catch d3err
    fprintf('4-D error. Out of options. \n')
  end
  D=1;
end

if D == 2
  boundsA = getBounds(groupA);
  boundsB = getBounds(groupB);
  
  diffx = (boundsA.xmin-boundsB.xmin) + (boundsA.xmax-boundsB.xmax);
  diffy = (boundsA.ymin-boundsB.ymin) + (boundsA.ymax-boundsB.ymax);
  diffxy = (boundsA.xmin-boundsB.ymin) + (boundsA.xmax-boundsB.ymax);
  
  if diffx > diffxy | diffy > diffxy
    xA = groupA(:,1); yA = groupA(:,2); zA = groupA(:,3);
    newA = [xA,yA,zA];
    newB = groupB;
    fix = 'XY';
  end
end

if fix ~= 0
  fixA = newA; fixB = newB;
else
  fixA = groupA; fixB = groupB;
end


end
