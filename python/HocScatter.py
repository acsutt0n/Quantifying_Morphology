# scatterplot 

def hocScatter(fName):
  import numpy as np
  import matplotlib.pyplot as pyplot
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  
  if wwrite == 1:
    outfile = open(outfileName, 'w')
  
  Xs = Ys = Zs = [0.0]
  
  # parse hoc input and plot
  # fig = plt.figure()
  # ax = fig.add_subplot(111, projection='3d')
  
  with open(fName, 'r') as infile:
    for line in infile:
      cols = line.split('(')
      
      if cols[0] == '  pt3dadd':
        coords = cols[1].split(',')
        
   #     if not Xs:
        Xs.append(float(coords[0]))
   #     else:
   #       Xs = float(coords[0])
   #     if Ys:
        Ys.append(float(coords[1]))
   #     else:
   #       Ys = float(coords[1])
   #     if Zs:
        Zs.append(float(coords[2]))
   #     else:
   #       Zs = float(coords[2])
        
        
     #   new_coords = [x, y, z]
     #   scatter = np.vstack((scatter, new_coords))
        if wwrite == 1:
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
