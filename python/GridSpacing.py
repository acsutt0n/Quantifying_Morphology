# GridSpacing.py -- this program analyzes the gridpoints created from
#                   a fitted ellipsoid to a neuropil and the hoc file
#                   used to create them. It examines how close each
#                   grid point is to the closest hoc skeleton point.
#
# usage: python GridSpacing.py hocFile.hoc gridpoints.txt
#

import numpy as np
import matplotlib.pyplot as plt
from NeuronGeometry import *
from getSkeleton import *
import sys
# import numbapro as pro


# python 2.7 holdover compatibility
xrange=np.arange


def loaddata(hocFile, gridFile):
  # import all the hocFile nodes
  # this returns 3 nodes per segment; need them all
  nodes = getSkeleton(hocFile,0) # all other arguments are optional
  print('Loaded %i nodes from %s' %(len(nodes), hocFile))
  newnodes = []
  for i in range(len(nodes)):
    newnodes.append(nodes[i][:3])
  nodes = newnodes
  
  grids = []
  lineNum = -1
  # import all gridpoints
  with open(gridFile, 'r') as fIn:
    for line in fIn:
      lineNum = lineNum + 1
      splitLine = line.split(None)
      grids.append([float(i) for i in splitLine])
  print('Loaded %i gridpoints from %s' %(len(grids), gridFile))
  
  return nodes, grids

"""
# GPU
@pro.vectorize(["float32(float32[:], float32[:])"])
def returnDist(point, pts):
  # vectorized with GPU
  dists = []
  for n in range(len(pts)):
    dists.append(np.linalg.norm(point - pts[n]))
  return min(dists)
"""



def closestPoint(rectpoint, nodes):
  # find the closest neuron node to a rectangle point
  dists = []
  for n in range(len(nodes)):
    dists.append(np.linalg.norm(rectpoint - nodes[n]))
  # print('Dist %.3f found' %dist)
  return min(dists)



def distanceHisto(nodes, grids):
  # make a histogram of how far each gridpoint is from each neuron pt
  dists = []
  nodes, grids = np.array(nodes), np.array(grids)
  # dists = [closestPoint(i, nodes)[2] for i in grids]
  
  for i in range(len(grids)):
    curr_gpt = grids[i]
    curr_dist = closestPoint(curr_gpt, nodes)
    # curr_dist = returnDist(curr_gpt, nodes)
    dists.append(curr_dist)
    if i%100 == 0:
      print('%i/%i done' %(i, len(grids)))
  
  print('Found %i distances' %len(dists))
  return dists
  


def writeHistogram(dists, newFile='histogram_distances.txt'):
  # write histo distances
  with open(newFile, 'w') as fOut:
    for i in range(len(dists)):
      fOut.write(str(dists[i]))
      fOut.write('\n')
  fOut.close()
  print('%s written' %newFile)
  return



def gridspacingcontrol(hocFile, gridFile, options=False):
  nodes, grids = loaddata(hocFile, gridFile)
  print('Nodes is length: %i' %len(nodes))
  print('Grids is length: %i' %len(grids))
  dists = distanceHisto(nodes, grids)
  writeHistogram(dists, '786_histogram_dists.txt')
  
  return




########## control ##########
if __name__ == "__main__":
  arguments = sys.argv
  hocFile, gridFile = arguments[1], arguments[2]
  gridspacingcontrol(hocFile, gridFile)
