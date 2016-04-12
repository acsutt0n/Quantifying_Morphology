# neuron_skeletonPlot.py - a suite for plotting skeletons -- 
# but not dendrograms; see radial_dendrogram.py and axon_dendrogram.py
# for dendrogram stuff

import matplotlib.pyplot as plt
from mpltoolkits.mplot3d import Axes3D
import numpy as np



#######################################################################
# NetworkX skeletons (very coarse except for knossos skeletons)
import networkx as nx

def get_nodes(geo):
  # Imports "nodes" -- really just midpoints of sections and their connections
  nodes, sources, targets = {}, [], []
  for s in geo.segments:
    nodes[len(nodes.keys())] = (s.coordAt(0))[:2] # sometimes coordAt(0.5) doesn't exist
    for n in s.neighbors:
      sources.append(s)
      targets.append(n)
  # get index for neighbors
  for s in range(len(geo.segments)):
    for t in range(len(targets)):
      if targets[t] == geo.segments[s]:
        targets[t] = s
    for j in range(len(sources)):
      if sources[j] == geo.segments[s]:
        sources[j] = s
  # should now have all ints for segments
  betch = [1 if type(i) == int else 0 for i in sources]
  betch2 = [1 if type(i) == int else 0 for i in targets]
  if np.prod(betch) == 0 or np.prod(betch2) == 0:
    print('Error! Some elements of sources or targets still not integers!')
    return None
  if len(targets) != len(sources):
    print('Error! Length of targets should equal length of sources!')
    return None
  return nodes, sources, targets


def simple_nxplot(geo,title=None):
  nodes, sources, targets = get_nodes(geo)
  g = nx.Graph()
  for e in range(len(sources)):
    g.add_edge(sources[e], targets[e])
  edgelist = list(zip(sources, targets))
  nxobj = nx.draw_networkx_edges(g, pos=nodes, edgelist=edgelist)
  plt.title(title)
  plt.show()
  return
  

def many_skeletons(geofiles, labels=None):
  sizes = {2: [211,212], 3: [221,222,223], 4: [221,222,223,224],
           5: [321,322,323,324,325], 6: [321,322,323,324,325,326],
           7: [331,332,333,334,335,336,337],
           8: [331,332,333,334,335,336,337,338],
           9: [331,332,333,334,335,336,337,338,339]}
  if len(geofiles) > 9:
    print('Only using first 9 (of %i) geo files' %len(geofiles))
    geofiles = geofiles[:9]
  fig = plt.figure()
  for range(len(geofiles)):
    ax = fig.add_subplot(sizes[len(geofiles)][g])
    nodes, sources, targets = get_nodes(geofiles[g])
    return
  


#######################################################################
# 3-D skeleton plotting, plotting all the nodes

def simple_skeleton(geo):
  return # Not finished






























