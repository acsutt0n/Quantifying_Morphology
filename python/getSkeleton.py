# getSkeleton.py -- get skeleton from hoc
"""
usage: $ python getSkeleton.py hocFile skeletonFile
"""

def getSkeleton(hocFile, skeletonFile=None, wwrite=0):
  
  points = []
  
  if wwrite == 1:
    outfile = open(skeletonFile, 'w')
  
  
  ### parse hoc file
  with open(hocFile, 'r') as infile:
    
    openSection = 0 # initialize openSection
    
    for line in infile:
      cols = line.strip().split('_')
      
      # start of new section
      if cols[0] == 'section':
        sec_cols = cols[1].split('{')
        current_section = int(sec_cols[0])
        openSection = 1
      
      elif cols[0] == '}' and openSection:
        openSection = 0
        current_section = -1
      
      else:
        pt3d_cols = cols[0].split('(')
        
        # find new point
        if pt3d_cols[0] == 'pt3dadd' or pt3d_cols[0] == '  pt3dadd' \
          and openSection and current_section >= 0:
          point_cols = pt3d_cols[1].split(',')
          X = point_cols[0]
          Y = point_cols[1]
          Z = point_cols[2]
          rad_cols = point_cols[3].split(')')
          rad = rad_cols[0]
          points.append([float(i) for i in [X,Y,Z,rad]])
          
          if wwrite == 1:
            # log new point
            print(current_section)
            skel_line = [current_section, X, Y, Z, rad]
            skel_string = ['0','0','0','0','0']
            for c in range(len(skel_line)):
              skel_string[c] = str(skel_line[c])
            printskel_string = ' '.join(skel_string)
            outfile.write(printskel_string)
            outfile.write('\n')
  
  return points



def simpleSkeleton(hocfile, skelfile):
  # No section logging, just the pt3dadd points.
  with open(skelfile, 'w') as fOut:
    with open(hocfile, 'r') as fIn:
      for line in fIn:
        if line:
          splitline = line.split(None)
          if splitline:
            if line.split('(')[0] == 'pt3dadd' or line.split('(')[0] == '  pt3dadd':
              stuffs = line.split(',')
              x = stuffs[0].split('(')[1]
              y = stuffs[1]
              z = stuffs[2]
              rad = stuffs[3].split(')')[0]
              fOut.write(' '.join([x, y, z, rad]))
              fOut.write('\n')
  return

  
#########################################################
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]
  if len(sys.argv) > 2:
    skeletonFile = sys.argv[2]
  else:
    skeletonFile = 'skeleton_points.txt'
  
  getSkeleton(hocFile, skeletonFile, wwrite=1)
