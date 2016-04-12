
## imports
import numpy as np
import matplotlib.pyplot as plt
from math import *
from neuron_readExportedGeometry import *


##################################### dev ################################
geo = demoRead('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/803_151_63x_IM_69.hoc')
pts = [0,0]
connections, bounds = [], []
axonInds, axonTipPos = geo.getAxonIndices()
pDF = PathDistanceFinder(geo, geo.soma, 1)
path = pDF.pathTo(geo.segments[axonInds[1]])


## ALTERNATIVE
# get it started
pts = [[0,0]]
bounds = [[-0.5, 0.5]]
connections = []



#################################### helpers ##############################
def dist_to(geo, seg1, seg2, seg1end, seg2end):
  # segs must be class Segment instances
  pDF = PathDistanceFinder(geo, seg1, seg1end)
  return pDF.distanceTo(seg2, seg2end)


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


def concat_points(new_pts, new_bounds, prevpt, pts, connections, bounds):
  for p in range(len(new_pts)):
    pts.append(new_pts[p])
    connections.append([pts.index(prevpt),pts.index(new_pts[p])])
    bounds.append(new_bounds[p])
  return pts, connections, bounds



###########################################################################
def get_layers_path(geo, path, N):
  # given geo and a specific path segment, get the layers
  # path = list of segment objects, n = int (which seg is current)
  prev = path
  potential_seglist = path[N].neighbors
  seglist = []
  layers, count = [], -1
  lengths = [] # each layer has a length from the axon
  retired_segs = []
  
  for s in potential_seglist:
    if s not in path:
      seglist.append(s)
  
  while len(seglist) > 0:
    count = count + 1
    newseglist = []
    layers.append([])
    lengths.append([]) # lengths = [ 
    for s in range(len(seglist)):
      layers[count].append([]) # add the count-th layer to layers
      lengths[count].append([])
      nc = 0
      for n in seglist[s].neighbors:
        if n not in newseglist and n not in seglist and n not in retired_segs:
          newseglist.append(n)
          nc = nc + 1
        if nc > 0:
          layers[count][s].append(nc)
          # getting lengths is the super-slowest part (>99% of time)
          #lengths[count][s].append(dist_to(geo, path[N], n, 0, 0))
    
    print(len(layers[-1]))
    if len(newseglist) == len(seglist):
      break
    for i in seglist:
      retired_segs.append(i)
    seglist = [i for i in newseglist]
  return layers, lengths
      



def stem(old_bounds, num, dists):
  # given the old bounds and the number of requested bound, 
  # give new pts and bounds
  if type(dists) is not list:
    print(dists)
    print('This version takes multiple radii, one for each pt')
    return
  num = 2*num + 1
  npt = np.linspace(old_bounds[0], old_bounds[1], num)
  npts = [npt[i] for i in range(len(npt)) if i%2 != 0] # these are new pts x-vals
  new_bounds = []
  for i in range(len(npt)-2):
    if i % 2 == 0:
      new_bounds.append([npt[i],npt[i+2]])
  new_pts = [[npts[i],dists[i]] for i in range(len(npts))]
  return new_pts, new_bounds


##########################################################################
#                          money functions                               #
##########################################################################

def next_layer(layers, lengths, n, pts, bounds, connections):
  completed = many_completed(layers, n)
  print('Completed is: %i' %completed)
  for l in range(len(layers[n])):
    if layers[n][l]:
      prevpt = pts[completed-1]
      print(len(pts),pts) ##
      print(layers[n][l][0], lengths[n][l][0])
      new_pts, new_bounds = stem(bounds[completed-1], 
                                 layers[n][l][0], 
                                 lengths[n][l])
      pts, connections, bounds = concat_points(new_pts, new_bounds, prevpt,
                                               pts, connections, bounds)
  return pts, bounds, connections



def base_layer(layers, lengths, base, layer):
  



def create_pts(layers, lengths):
  # list of layers, list of lengths
  pts = [[0,0]]
  pcount = 0
  connections = []
  bounds = []
  for l in range(len(layers)):
    # temp_bounds = []
    pcount = pcount + 1
    # make new base point
    pts.append([pcount, 0])
    # add the new connection
    connections.append([0, len(pts)-1])
    # start the base bounds
    bounds.append([pcount-0.5, pcount+0.5])
    # initiate previous bounds
    prevbound = bounds[-1]
    # now set up all the layers for this given base point
    for m in range(len(layers[l])):
    # there are m layers for the l-th base point
      if len(layers[l][m]) > 0: # if this layer branches
        # create new bounds
        #newbounds = np.linspace(prevbounds[0], prevbounds[1],
        #                        len(layers[l][m]))
        print('Calling next_layer with args: %i, %i, %i, %i, %i, %i' \
              %(m, m, m, len(pts), len(bounds), len(connections)))
        pts, bounds, connections = next_layer(layers[l], lengths[l], m,
                                              pts, bounds, connections)
  return pts, bounds, connections
        


######################### plotting ###############################3
def plot_radial_dend(connections, pts, colors=None):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(-1,0, c='k', s=50)
  ax.plot([-1,0],[0,0], c='k')
  for c in connections:
    if colors:
      col = colors[connections.index(c)]
    else:
      col='k'
    ax.plot([pts[c[0]][0], pts[c[1]][0]], 
            [pts[c[0]][1], pts[c[1]][1]], c=col) # omit c for multicolored (kinda fun)
  plt.show()



################### development as of 4.10.2015 #######################
from neuron_readExportedGeometry import *
geo = demoRead('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/803_151_63x_IM_69.hoc')
axonInds, axonTips = geo.getAxonIndices()
pDF = PathDistanceFinder(geo, axonInds[1], axonTips[1])
path = pDF.pathTo(geo.soma)

pts = [[0,0],[0,1],[0,2],[-.25,3],[.25,3]]   
bounds=[[-.5,.5],[-.5,.5],[-.5,.5],[-.5,0],[0,.5]]
connections  = [[0,1],[1,2],[2,3],[2,4]]
pts, bounds, connections = next_layer(layers0, 2,pts, bounds, connections) 
## here there are some problems, the graph looks funky





"""
####################################### dev layers
def get_layers_2(path, N):
  seglist = [path[N]]
  layers, count = [], -1
  breakit = 0
  while len(seglist) > 0:
    count = count + 1
    sofar, newseglist = 0, []
    layers.append([])
    for s in range(len(seglist)):
      layers[count].append([]) # add the s-th list
      nc = 0
      for n in seglist[s].neighbors:
        if n not in newseglist and n not in path:
          newseglist.append(n)
          nc = nc + 1
      if nc > 0:
        layers[count][s].append(nc)
        
    print(len(layers[-1]))
    if len(newseglist) == len(seglist) and breakit < 2:
      breakit = breakit+1
    else:
      breakit = 0
    if breakit == 2:
      break
    
    seglist = [i for i in newseglist]
  return layers

## toy function
pcount = 0
for p in path:
  pcount = pcount + 1
  nebs = p.neighbors
  for n in nebs:
    if n in path:
      nebs.pop(nebs.index(n))
  # only go one direction, lose neighbors also in path
  pts.append([pcount,0])
  connections.append([0,len(pts)-1])
  bounds.append([pcount-0.5, pcount+0.5])
  # cycle through neighbors and add pts and connections
  while len(nebs) > 0:
    for n in nebs:
      
  
"""
