# connect nearest nodes as segments instead



def add_pts(pts, splitline):
  



def load_hoc_pts(filein):
  pts = []
  with open(filein, 'r') as fIn:
    for line in fIn:
      if line:
        splitline = line.split(None)
        if len(splitline) >= 1:
          if splitline[0].split('(').[0] == 'pt3dadd':
            pts = add_pts(pts, splitline)
  return pts



