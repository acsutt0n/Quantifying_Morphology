# scatterplot 

import numpy as np
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D



class hocScatter(fName, outFileName=False, options=False):
  
  def __init__(self, fName):
    # initialize
    self.inFile = fName
    self.outFile = outFileName
    
    # run stuff
    # populate points
    self.readfile()
    
  
  def readfile(self):
    if self.outFile is not False:
      outfile = open(outfileName, 'w')
    else:
    Xs = Ys = Zs = [0.0]
    
    # parse hoc input and plot
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    
    with open(fName, 'r') as infile:
      for line in infile:
        splitLine = line.split(None)
        cols = splitLine.split('(')
        
        if cols[0] == 'pt3dadd':
          coords = cols[1].split(',')
          
          Xs.append(float(coords[0]))
          Ys.append(float(coords[1]))
          Zs.append(float(coords[2]))
          
          if self.outFile is not False:
            outfile.write(','.join([coords[0],coords[1],coords[2]]))
            outfile.write('\n')
          # ax.scatter(Xs[-1], Ys[-1], Zs[-1])
          # print('points added')
  
  # plot scatter
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for i in xrange(len(Xs)):
    ax.scatter(Xs[i], Ys[i], Zs[i])

  """ I tried to exclude the first point (0,0,0) by doing Xs[1:-1]
      but I think that actually sorted the array instead....
  """
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  
  plt.show()


if __name__ == '__main__':
  import sys
  fName = sys.argv[1]
  if len(sys.argv) < 3:
    wwrite = 0
  else:
    wwrite = 1
    outfileName = sys.argv[2]
  hocScatter(fName)
