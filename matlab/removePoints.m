function newpts = removePoints(oldpts, xlim, ylim, zlim)
% usage: newpts = removePoints(oldpts, xlim (opt), ylim (opt), zlim (opt))

% removes points that lay beyond (greater than) a lim, 
% unless the lim is (-) then remove points less than that lim
% can give as many lims as want, in order x, y, z; use lim=1 to skip
% a limit and define a later one

if nargin > 1
  xlim = xlim;
else
  xlim = 1;
end

if nargin > 2
  ylim = ylim;
else
  ylim = 1;
end

if nargin > 3
  zlim = zlim;
else
  zlim = 1;
end

newpts = zeros(1,3);

%
for p = 1:length(oldpts)
  add = 1;
  % xlim first
  if xlim ~= 1
    % positive
    if xlim > 0
      if oldpts(p,1) < xlim
        xcurr = oldpts(p,1);
        add = add*1;
      else
        add = add*0;
      end
    % negative
    else
      if oldpts(p,1) > abs(xlim)
        xcurr = oldpts(p,1);
        add = add*1;
      else
        add = add*0;
      end
    end
  % no xlim bounds, so add point
  else
    add = add*1;
    xcurr = oldpts(p,1);
  end
  
  % ylim
  if ylim ~= 1
    % positive
    if ylim > 0
      if oldpts(p,2) < ylim
        ycurr = oldpts(p,2);
        add = add*1;
      else
        add = add*0;
      end
    % negative
    else
      if oldpts(p,2) > abs(ylim)
        ycurr = oldpts(p,2);
        add = add*1;
      else
        add = add*0;
      end
    end
  % no ylim bounds, so add point
  else
    add = add*1;
    ycurr = oldpts(p,2);
  end
  
  % zlim
  if zlim ~= 1
    % positive
    if zlim > 0
      if oldpts(p,3) < zlim
        zcurr = oldpts(p,3);
        add = add*1;
      else
        add = add*0;
      end
    % negative
    else
      if oldpts(p,3) > abs(zlim)
        zcurr = oldpts(p,2);
        add = add*1;
      else
        add = add*0;
      end
    end
  % no zlim bounds, so add point
  else
    add = add*1;
    zcurr = oldpts(p,3);
  end
  
  % if the point satisfies the requirements, add it
  if add == 1
    newpts(end+1, :) = [xcurr, ycurr, zcurr];
  end
end



