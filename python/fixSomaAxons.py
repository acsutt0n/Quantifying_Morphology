# show_skeleton.py -- can show a downsampled skeleton with or 
# without the soma highlighted
# usage: python show_skeleton.py hocFile


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from neuron_readExportedGeometry import *
import sys



def default_options(geo):
  # this does nothing right now
  numsegs = len(geo.segments)
  numnodes = len(geo.nodes)
  if numnodes / 10 >= numsegs * 3:
    how = 'segs'
  else:
    how = 'nodes'
  ops = {'how': how}
  return



def skeleton(geo, version='segs'):
  # return some reduced set of points based on the version for plotting
  if version == 'segs':
    pts = []
    for s in geo.segments:
      pts.append(s.coordAt(0))
      pts.append(s.coordAt(0.5))
      pts.append(s.coordAt(1))
  elif version == 'nodes':
    pts = []
    count = 0
    for s in geo.segments:
      for n in s.nodes:
        count = count + 1
        if count % 10 == 0:
          pts.append([n.x, n.y, n.z])
  soma = geo.soma.coordAt(0.5)
  return pts, soma



def axon_skeleton(geo, version='segs'):
  # return a reduced set of points for plotting plus axons and soma labeled
  axons = []
  for a in geo._axons:
    axons.append(a)
  # add reduced segment points but omit segments that are possible axons
  if version == 'segs':
    pts = []
    for s in geo.segments:
      if s not in axons:
        pts.append(s.coordAt(0))
        pts.append(s.coordAt(0.5))
        pts.append(s.coordAt(1))
  elif version == 'nodes':
    pts = []
    count = 0
    for s in geo.segments:
      if s not in axons:
        for n in s.nodes:
          count = count + 1
          if count % 10 == 0:
            pts.append([n.x, n.y, n.z])
  soma = geo.soma.coordAt(0.5)
  return axons, soma, pts
  


def plot_soma(hocgeo, version='segs'):
  # just plot the graph with a big red soma for geo.soma
  if type(hocgeo) is str: # if it's a filename (string)
    geo = demoRead(hocgeo)
  else: # assume it's a hocgeometry class
    geo = hocgeo
  pts, soma = skeleton(geo, version)
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for p in pts:
    ax.scatter(p[0],p[1],p[2], 'b.', edgecolor='b', s=1, alpha=0.1)
  ax.scatter(soma[0], soma[1], soma[2], 'ro', edgecolor = 'r', s=50)
  plt.show()
  return



def plot_possibles(hocgeo, which_soma=None, version='segs', axis='equal'):
  # it's likely the program thinks the true soma is an axon, so label
  # the axons and soma so the user can choose
  if type(hocgeo) is str: # if it's a filename (string)
    geo = demoRead(hocgeo)
  else: # assume it's a hocgeometry class
    geo = hocgeo
  axons, soma, pts = axon_skeleton(geo, version)
  if which_soma:
    for s in geo.segments:
      if s.name == which_soma:
        soma = s.coordAt(0.5)
        soma_name = s.name
    if not soma_name:
      print('Could not find chosen soma; reverting to default')
      soma_name = 'default'
  else:
    soma_name = geo.soma.name
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for p in pts:
    ax.scatter(p[0],p[1],p[2], 'b.', edgecolor='b', s=1, alpha=0.1)
  ax.scatter(soma[0], soma[1], soma[2], 'ro', edgecolor = 'r', s=50)
  ax.text(soma[0], soma[1], soma[2], soma_name, 'x')
  for a in axons:
    for n in a.nodes:
      ax.scatter(n.x, n.y, n.z, 'k.', s=2, alpha=0.5)
    ax.text(a.coordAt(0.5)[0], a.coordAt(0.5)[1], a.coordAt(0.5)[2],
            a.name, 'x')
  ax.axis(axis)
  plt.show()
  return




















