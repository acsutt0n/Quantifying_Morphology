import numpy as np
import matplotlib.pyplot as plt


from math import cos, sin
from numpy import r_
from numpy import atleast_2d as a2d
def from_polar(p):
  """(theta, radius) to (x, y)."""
  return _a(cos(p[0])* p[1], sin(p[0])* p[1])

def _a(a0, a1):
  return r_[a2d(a0), a2d(a1)]

def to_polar(c):
  """(x, y) to (theta, radius)."""
  return _a(arctan2(c[1], c[0]), (c** 2).sum(0)** .5)



# working version
def get_layers(geo):
  seglist = [geo.soma]
  layers, count = [], -1
  while len(seglist) > 0:
    count = count + 1
    sofar, newseglist = 0, []
    layers.append([])
    for s in range(len(seglist)):
      layers[count].append([]) # add the s-th list
      nc = 0
      for n in seglist[s].neighbors:
        if n not in newseglist:
          newseglist.append(n)
          nc = nc + 1
      if nc > 0:
        layers[count][s].append(nc)
        
    print(len(layers[-1]))
    if len(newseglist) == len(seglist):
      break
    seglist = [i for i in newseglist]
  return layers



def angle_span(angle, num):
  return np.linspace(angle[0],angle[1],num)



def stem(ang_bounds, num, rad):
  # return new base points
  num = num + num + 1
  ang = np.linspace(ang_bounds[0], ang_bounds[1],num)
  angs = [ang[i] for i in range(len(ang)) if i%2 != 0]
  new_bounds = []
  for i in range(len(ang)-2):
    if i % 2 == 0:
      new_bounds.append([ang[i],ang[i+2]])
  new_pts = [[i,rad] for i in angs]
  return new_pts, new_bounds
  


def concat_points(new_pts, new_bounds, prevpt, pts, connections, bounds):
  for p in range(len(new_pts)):
    pts.append(new_pts[p])
    connections.append([pts.index(prevpt),pts.index(new_pts[p])])
    bounds.append(new_bounds[p])
  return pts, connections, bounds



def many_completed(layers, l):
  # l is current layer, so layers[:l] are completed
  comp = 1 # start at 1 for origin
  for j in layers[:l]:
    for i in j:
      if type(i) is int:
        comp = comp + i
      elif type(i) is list:
        for k in i:
          if k:
            if type(k) is int:
              comp = comp + k
            elif type(k) is list:
              for h in k:
                if h:
                  if type(h) is int:
                    comp = comp + h
  return comp



##########################################################################
#                       Start here for constant length                   #
##########################################################################

from neuron_readExportedGeometry import *
import numpy as np
import matplotlib.pyplot as plt
geo = demoRead('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/803_151_63x_IM_69.hoc')
layers = get_layers(geo)
from numpy import pi



## for 2-point starters:
pts = [ [0,0],[3*pi/2,1],[pi/2,1] ]
connections = [[0,1],[0,2]]
prevpt = [pi/2,1]
bounds = [ [0,2*pi], [pi,2*pi], [0,pi] ]

## for 4-point starters
pts = [ [0,0],[pi/4,1],[3*pi/4,1],[7*pi/4,1], [5*pi/4,1] ]
connections = [[0,1],[0,2], [0,3], [0,4]]
prevpt = [pi/4,1]
# new_pts, new_bounds = stem([0,pi],layersB[1][0][0], 2)
bounds = [[0,2*pi],[0,pi/2],[pi/2,pi],[3*pi/2,2*pi],[pi,3*pi/2]]



# sample
pts = [ [0,0],[pi/4,1],[3*pi/4,1],[7*pi/4,1], [5*pi/4,1] ]
connections = [[0,1],[0,2], [0,3], [0,4]]
prevpt = [pi/4,1]
new_pts, new_bounds = stem([0,pi],layers[1][0][0], 2)
bounds = [[0,2*pi],[0,pi/2],[pi/2,pi],[3*pi/2,2*pi],[pi,3*pi/2]]

for p in range(len(new_pts)):
  pts.append(new_pts[p])
  connections.append([pts.index(prevpt),pts.index(new_pts[p])])
  bounds.append(new_bounds[p])
  
prevpt = pts[2]
new_pts, new_bounds = stem(bounds[2], layers[1][1][0],2)
pts, connections, bounds = concat_points(new_pts, new_bounds, prevpt, pts,
                                         connections, bounds)

completed = many_completed(layers, 1) # working on pts for layer 2
for l in range(len(layers[2])):
  if layers[2][l]:
    prevpt = pts[completed+l]
    new_pts, new_bounds = stem(bounds[completed+l], layers[2][l][0], 3)
    pts, connections, bounds = concat_points(new_pts, new_bounds, prevpt,
                                             pts, connections, bounds)



####################### money function ######################
def next_layer(layers, n, pts, bounds, connections):
  completed = many_completed(layers, n-1)
  for l in range(len(layers[n])):
    if layers[n][l]:
      prevpt = pts[completed+l] # changed from pts[completed+l]
      new_pts, new_bounds = stem(bounds[completed+l], layers[n][l][0], n+1)
      pts, connections, bounds = concat_points(new_pts, new_bounds, prevpt,
                                               pts, connections, bounds)
  return pts, bounds, connections





##################### plotting #######################
# convert to rectalinear points from polar
rpts = [from_polar(p) for p in pts]
# colors:
cols = ['b','g','r','c','m','y','k']*40 # 7*40 must be less than max radius
gols = [cols[p[-1]] for p in pts] # polar points

def plot_radial_dend(connections, pts, colors=None):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(0.5,0, c='k', s=50)
  ax.plot([0.5,0],[0,0], c='k')
  for c in connections:
    if colors:
      col = colors[connections.index(c)]
    else:
      col='k'
    ax.plot([pts[c[0]][0], pts[c[1]][0]], 
            [pts[c[0]][1], pts[c[1]][1]], c=col) # omit c for multicolored (kinda fun)
  plt.show()



###########################################################################
# Simple skeleton plot (for eve and morpho paper)
##########################################################################

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
  
  
  
  



# bad version ################################
def get_layers(geo):
  seglist = [geo.soma]
  layers, count = [], -1
  while len(seglist) > 0:
    count = count + 1
    sofar, newseglist = 0, []
    layers.append([])
    for s in range(len(seglist)):

      nc = 0
      for n in seglist[s].neighbors:
        if n not in newseglist and n not in seglist:
          newseglist.append(n)
          nc = nc + 1
          layers[count].append([]) # add the s-th list
      if nc > 0:
        if len(layers[count]) < s+1:
          layers[count].append([])
        else:
          layers[count][s].append(nc)
        
    print(len(layers[-1]))
    if len(newseglist) == len(seglist):
      break
    seglist = [i for i in newseglist]
  return layers







