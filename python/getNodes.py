# getNodes.py -- get the nodes from a hoc file, then do something 
#                with them
# usage:         python getNodes.py hocFile


import numpy as np
from loadHoc import *








def nodesControl(hocFile):
  """
  Stuff.
  """
  geometry = HocGeometry(hocFile)
  print('Num. of sections: %i' %len(geometry.secRads.keys()))
  
  return




#########################################
if __name__ == '__main__':
  arguments = sys.argv
  hocFile = arguments[1]
  nodesControl(hocFile)
