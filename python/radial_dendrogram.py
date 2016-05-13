# radial_dendrogram.py
# This creates a radial dendrogram of a hoc/geometry object
# Lengths are constant for the moment, might try to change this
#   in the future



import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin
from numpy import r_
from numpy import atleast_2d as a2d
from neuron_readExportedGeometry import *
from numpy import pi



########################################################################
# Smaller functions 

def angle_span(angle, num):
  """ """
  return np.linspace(angle[0],angle[1],num)

def from_polar(p):
  """(theta, radius) to (x, y)."""
  return _a(cos(p[0])* p[1], sin(p[0])* p[1])

def _a(a0, a1):
  return r_[a2d(a0), a2d(a1)]

def to_polar(c):
  """(x, y) to (theta, radius)."""
  return _a(arctan2(c[1], c[0]), (c** 2).sum(0)** .5)


def get_branch(geo, seg):
  return [branch for branch in geo.branches if seg.name in branch.tags][0]



def branch_layers(geo, nlayers=None, id_axons=None):
  """ 
  Same as above but uses branches, much more better.
  branchlist    : branches to be examined for neighbors
  used          : branches that are spent (already examined)
  newbranchlist : neighbors that will be checked upon next iteration
                  when they become branchlist
  """
  axon_dict = []
  if nlayers is None:
    nlayers = np.inf # Get all layers
  branchlist = [get_branch(geo, geo.soma)]
  layers, used, cnt = [], [get_branch(geo, geo.soma)], -1
  while len(branchlist) > 0 and cnt <= nlayers:
    cnt += 1 # Layer count
    newbranchlist = []
    layers.append([]) # Add a new layer
    
    for s in range(len(branchlist)): # Go through each of the "new" branches
      used.append(branchlist[s]) # Don't want this one included in future neighbors
      layers[cnt].append([]) # Add the s-th item of this layer
      nebcnt = 0 # Which is initialized with 0 branches (neighbors)
      
      for n in branchlist[s].neighbors: # For each neighboring branch
        if id_axons is None: # Use native axon tags
          if 'Axon' in n.tags:
            axon_dict.append([cnt, s, geo.branches.index(n)])
          else: # Use the geo.branches index
            for ax in id_axons:
              if geo.branches.index(n) == ax:
                axon_dict.append([cnt, s, geo.branches.index(n)])
        if n not in used and n not in newbranchlist:
          newbranchlist.append(n)
          nebcnt += 1 # Increment neighbor count
      if nebcnt > 0:
        layers[cnt][s].append(nebcnt) # Add any new layers
    
    # print(len(layers[-1])) # Still in while loop
    if len(newbranchlist) == len(branchlist): # No new layers found
      break
    # Else, replace branchlist and continue to next layer
    branchlist = [i for i in newbranchlist]
  return layers, axon_dict





def stem(ang_bounds, num, rad):
  """"
  """
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
  """
  """
  for p in range(len(new_pts)):
    pts.append(new_pts[p])
    connections.append([pts.index(prevpt),pts.index(new_pts[p])])
    bounds.append(new_bounds[p])
  return pts, connections, bounds



def many_completed(layers, l):
  """
  """
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

# RADIAL

def get_started_radial(geo, rrange=[0,2*pi]):
  """
  rrange always starts at 0.
  """
  # Figure out how many branches come off the the primary neurite
  used, nebs, currseg = [], [geo.soma], geo.soma
  while len(nebs) < 2:
    currseg = nebs[0] # len(nebs) must be 1
    used.append(currseg)
    nebs = [] # Only un-used segments count as neighbors here
    for n in currseg.neighbors:
      if n not in used:
        nebs.append(n)
    # Should check len(nebs), is == 1 continue, else exit loop
  
  if len(nebs) == 2: # for 2-point starters:
    pts = [ [0,0],[3*(rrange[1])/4,1],[1*rrange[1]/4,1] ]
    connections = [[0,1],[0,2]]
    prevpt = [pi/2,1]
    bounds = [ rrange, [rrange[1]/2,rrange[1]], [0, rrange[1]/2] ]
    
  elif len(nebs) == 3: # for 3-point starters
    pts = [ [0,0], [(rrange[1])/3,1], [2*(rrange[1])/3,1], [0,1] ]
    connections = [[0,1],[0,2], [0,3]]
    prevpt = [2*pi/3,1]
    bounds = [ rrange, [0.5*(rrange[1])/3, 1.5*(rrange[1])/3], 
                       [(rrange[1])/2, 2.5*(rrange[1])/3], 
                       [2.5*(rrange[1])/3, 0.5*(rrange[1])/3] ]
  
  elif len(nebs) == 4: # for 4-point starters
    pts = [ [0,0],[1*(rrange[1])/8,1],[3*(rrange[1])/8,1],
            [5*(rrange[1])/8,1], [7*(rrange[1])/8,1] ]
    connections = [[0,1],[0,2], [0,3], [0,4]]
    prevpt = [pi/4,1]
    # new_pts, new_bounds = stem([0,pi],layersB[1][0][0], 2)
    bounds = [rrange, [0, 2*(rrange[1])/8], [2*(rrange[1])/8, 4*(rrange[1])/8],
                      [4*(rrange[1])/8, 6*(rrange[1])/8],
                      [6*(rrange[1])/8, 0]]
  
  else:
    print("Don't know how to handle %i starting points" %len(nebs))
    return None
  return pts, connections, prevpt, bounds



def next_layer(layers, n, pts, bounds, connections):
  """ 
  Map up to _n_ layers.
  """
  completed = many_completed(layers, n-1)
  for l in range(len(layers[n])):
    if layers[n][l]:
      prevpt = pts[completed+l] # changed from pts[completed+l]
      new_pts, new_bounds = stem(bounds[completed+l], layers[n][l][0], n+1)
      pts, connections, bounds = concat_points(new_pts, new_bounds, prevpt,
                                               pts, connections, bounds)
  return pts, bounds, connections



def radial_dendrogram(geo, nlayers=10, show=True, colors='trippy',
                      rrange=[0,2*pi], id_axons=None):
  """
  Can be run for many if type(geofiles) is list.
  """
  if type(geo) is str:
    geo = [demoReadsilent(geo)]
  # Get started
  pts, connections, prevpt, bounds = get_started_radial(geo, rrange)
  # Get the layers
  layers, axon_dict = branch_layers(geo, nlayers, id_axons)
  # Create the radial structure for _nlayers_
  for layer in range(1, len(layers)):
    pts, bounds, connections = next_layer(layers, layer, pts, bounds,
                                          connections)
  
  # Then plot
  # convert to rectalinear points from polar
  rpts = [from_polar(p) for p in pts]
  if colors is None:
    colors = ['k' for i in range(3000)] # Should be large enough
  if colors is 'trippy':
    colors = ['r', 'orange', 'y', 'g', 'b', 'purple']*600
  plot_radial_dend(connections, rpts, colors, axon_dict)
  
  return

# RECTANGULAR


def get_started_rect(geo, rrange=[0,2*pi]):
  """
  rrange always starts at 0.
  """
  # Figure out how many branches come off the the primary neurite
  used, nebs, currseg = [], [geo.soma], geo.soma
  while len(nebs) < 2:
    currseg = nebs[0] # len(nebs) must be 1
    used.append(currseg)
    nebs = [] # Only un-used segments count as neighbors here
    for n in currseg.neighbors:
      if n not in used:
        nebs.append(n)
    # Should check len(nebs), is == 1 continue, else exit loop
  
  if len(nebs) == 2: # for 2-point starters:
    pts = [ [0,0],[3*(rrange[1])/4,1],[1*rrange[1]/4,1] ]
    connections = [[0,1],[0,2]]
    prevpt = [pi/2,1]
    bounds = [ rrange, [rrange[1]/2,rrange[1]], [0, rrange[1]/2] ]
    
  elif len(nebs) == 3: # for 3-point starters
    pts = [ [0,0], [(rrange[1])/3,1], [2*(rrange[1])/3,1], [0,1] ]
    connections = [[0,1],[0,2], [0,3]]
    prevpt = [2*pi/3,1]
    bounds = [ rrange, [0.5*(rrange[1])/3, 1.5*(rrange[1])/3], 
                       [(rrange[1])/2, 2.5*(rrange[1])/3], 
                       [2.5*(rrange[1])/3, 0.5*(rrange[1])/3] ]
  
  elif len(nebs) == 4: # for 4-point starters
    pts = [ [0,0],[1*(rrange[1])/8,1],[3*(rrange[1])/8,1],
            [5*(rrange[1])/8,1], [7*(rrange[1])/8,1] ]
    connections = [[0,1],[0,2], [0,3], [0,4]]
    prevpt = [pi/4,1]
    # new_pts, new_bounds = stem([0,pi],layersB[1][0][0], 2)
    bounds = [rrange, [0, 2*(rrange[1])/8], [2*(rrange[1])/8, 4*(rrange[1])/8],
                      [4*(rrange[1])/8, 6*(rrange[1])/8],
                      [6*(rrange[1])/8, 0]]
  
  else:
    print("Don't know how to handle %i starting points" %len(nebs))
    return None
  return pts, connections, prevpt, bounds



def rect_dendrogram(geo, nlayers=10, show=True, colors='trippy',
                      rrange=[0,2*pi], id_axons=None):
  """
  Can be run for many if type(geofiles) is list.
  """
  if type(geo) is str:
    geo = [demoReadsilent(geo)]
  # Get started
  pts, connections, prevpt, bounds = get_started_rect(geo, rrange)
  # Get the layers
  layers, axon_dict = branch_layers(geo, nlayers, id_axons)
  # Create the radial structure for _nlayers_
  for layer in range(1, len(layers)):
    pts, bounds, connections = next_layer(layers, layer, pts, bounds,
                                          connections)
  
  # Then plot
  # convert to rectalinear points from polar
  rpts = [from_polar(p) for p in pts]
  if colors is None:
    colors = ['k' for i in range(3000)] # Should be large enough
  if colors is 'trippy':
    colors = ['r', 'orange', 'y', 'g', 'b', 'purple']*600
  plot_rect_dend(connections, rpts, colors, axon_dict)
  
  return







#################################################################### 
# plotting



def plot_radial_dend(connections, pts, colors=None, axon_dict=None):
  """
  Returns a plot -- but does not show it (allows for subplots).
  """
  #fig = plt.figure()
  # ax = fig.add_subplot(111)
  plt.scatter(0.5,0, c='k', s=50)
  plt.plot([0.5,0],[0,0], c='k')
  for c in connections:
    if colors:
      col = colors[connections.index(c)]
    else:
      col='k'
    plt.plot([pts[c[0]][0], pts[c[1]][0]], 
             [pts[c[0]][1], pts[c[1]][1]], c=col) # omit c for multicolored (kinda fun)



def many_radial_dend(geofiles, nlayers=None):
  """
  Plot many radial dendrograms in a subplot.
  """
  # 4 in a row
  dims = [int(len(geofiles)/4.+1),4]
  for g in range(len(geofiles)):
    plt.subplot(dims[0], dims[1], g+1)
    radial_dendrogram(geofiles[g], nlayers=nlayers, colors=None)
    plt.title(geofiles[g].name)
    # Get rid of side labels
    plt.axis('off')
  plt.show()






########################################################################

"""
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

"""



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
  for g in range(len(geofiles)):
    ax = fig.add_subplot(sizes[len(geofiles)][g])
    nodes, sources, targets = get_nodes(geofiles[g])
  return
  
  
  
  
#########################################################################

if __name__ == "__main__":
  if len(sys.argv) < 2:  # python radial(0) geofile(1)
    geo = demoReadsilent('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/803_151_63x_IM_69.hoc')
  else:
    try:
      geo = demoReadsilent(sys.argv[1])
    except:
      print("Couldn't use %s, using sample file instead" %sys.argv[1])
      geo = demoReadsilent('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/803_151_63x_IM_69.hoc')
  radial_dendrogram(geo, show=True_)
  
  
  


#########################################################################
#########################################################################
# old code -- NOT IN USE


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




def get_layers(geo):
  """
  """
  seglist = [geo.soma]
  layers, count = [], -1
  while len(seglist) > 0:
    count = count + 1 # Keep track of current layer
    sofar, newseglist = 0, []
    layers.append([]) # Add a new layer list
    for s in range(len(seglist)):
      layers[count].append([]) # add the s-th list
      nc = 0 # Neighbor count
      for n in seglist[s].neighbors: # For each segment ...
        if n not in newseglist:  #        ... if we haven't gotten its neighbors....
          newseglist.append(n)   #        ... add those neighbors
          nc = nc + 1  # increment the count
      if nc > 0: # If there were any new neighbors (new layer) ...
        layers[count][s].append(nc) #     ... add them
        
    print(len(layers[-1])) # Number of NEWEST neighbors found
    if len(newseglist) == len(seglist): # If we haven't found any new layers
      break
    # Else, re-wite seglist to prepare for next layer
    seglist = [i for i in newseglist]
  return layers





