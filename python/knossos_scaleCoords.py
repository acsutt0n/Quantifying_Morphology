# knossos_scaleCoords.py - scales the hoc pt3dadd coordinates
# usage: python knosso_scaleCoords.py hocfile majorvox minorvox zvox
#        Where _vox is the dimension of the voxel;
#        'x' and 'y' are usually the same, but this finds which is major
#         instead of differentiating 'x' and 'y'.



import numpy as np



################# Hoc loading functions ################

def add_point(line, points):
  h = line.split('(')[1]
  c = h.split(',')
  temp = [float(a) for a in c[:3]]
  temp.append(float(c[-1].split(')')[0]))
  points.append(temp)
  return points



def load_hoc(hocfile):
  """
  Load the hocfile and return an Nx4 (x,y,z,rad) list. Not currently used.
  """
  connections, points = [], []
  with open(hocfile, 'r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        if splitLine:
          if len(splitLine) == 1:
            l = splitLine[0]
            if l.split('(')[0] == 'pt3dadd':
              points = add_point(line, points)
          else:
            if splitLine[0] == 'connect':
              connections.append(line)
  return connections, points



################### Analysis ###################

def get_dims(points):
  """
  Find out which is x, y, etc. Not sure if this is necessary.
  """
  dims = {'xmin': np.inf, 'xmax': 0, 'ymin': np.inf, 'ymax': 0}
  return



def load_and_fix(hocfile, scales, outfile=None):
  """
  This version assumes a lot but it will probably work.
  """
  if len(scales) != 3:
    print('scales must have 3 values: x, y, z voxel dims')
    return
  scales = [float(i) for i in scales]
  if outfile is None:
    outfile = hocfile.split('.')[0] + '_scaled.hoc'
  #
  def fix_it(line, scales, fOut):
    # Fix the coordinates
    ll = line.split(',')
    pre, pt1 = ll[0].split('(')[0]+'(', float(ll[0].split('(')[1])
    pt2 = float(ll[1])
    pt3 = float(ll[2]) # Everything else is the same (rad1, rad2)
    rest = ll[3:]
    #print(pt1, scales[0], pt2, scales[1], pt3, scales[2])
    #print(type(pt1), type(scales[0]), type(pt2), type(scales[1]), 
    #      type(pt3), type(scales[2]))
    pt1, pt2, pt3 = pt1*scales[0], pt2*scales[1], pt3*scales[2]
    thing = [pre + str(pt1), str(pt2), str(pt3)]
    for i in rest:
      thing.append(i)
    fOut.write(', '.join(thing))
  # Go through each line and correct the coordinates
  with open(hocfile, 'r') as fIn:
    with open(outfile, 'w') as fOut:
      for line in fIn:
        if line:
          splitLine = line.split(None)
          # print(splitLine)
          if len(splitLine) >= 1:
            if splitLine[0].split('(')[0] == 'pt3dadd':
              fix_it(line, scales, fOut)
            else:
              fOut.write(line)
          else:
            fOut.write('\n')
  return
  



###########################################################################

if __name__ == "__main__":
  import sys
  args = sys.argv
  hocfile, scales = args[1], args[2:5]
  load_and_fix(hocfile, scales)
















