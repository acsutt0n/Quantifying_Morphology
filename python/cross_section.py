# cross_section.py -- optimize cross section using imageMatrix and spiral


from imageMatrix import *
from spiral import *
from neuron_readExportedGeometry import *



######################### helper functions #############################

def min_max(stack, isstack=True):
  def minime(slic):
    _min, _max = np.inf, 0
    for i in range(len(slic)):
      for j in range(len(slic[i])):
        obj = slic[i][j]
        if obj > _max:
          _max = obj
        if obj < _min:
          _min = obj
    return _min, _max
  if isstack:
    agg = [minime(slic) for slic in stack]
  else:
    agg = minime(stack)
  return agg
    

def beta_load(hocfile, folder):
  """
  Load stuff to keep it in memory so it doesn't go out of scope after
  function ends.
  """
  stack = import_images(folder, par=False)
  geo = demoReadsilent(hocfile)
  return stack, geo
  



def get_ten(geo, stack, N=10):
  """
  Returns ten raw cross-sections.
  """
  # Get the first N nodes
  def get_nodes(geo, N):
    nodes, ncount = [], -1
    for b in geo.branches:
      t_nodes = []
      for n in range(len(b.nodes)):
        if ncount >= N:
          return nodes
        else:
          t_nodes.append([b.nodes[n].x, b.nodes[n].y, b.nodes[n].z])
          ncount = ncount + 1
      nodes.append(t_nodes)
    print('something went wrong, num nodes: %i' %len(nodes))
    return None
  # 
  if type(geo) is list:
    nodes = geo[:N]
    llist = True
  else:
    nodes = get_nodes(geo, N)
    llist = False
  cross_secs, startpts = [], []
  if llist is False:
    for nlist in nodes:
      for n in range(len(nlist)-1):
        cs, startpt = get_cross_sec(nlist[n], nlist[n+1], stack)
        cross_secs.append(cs)
        startpts.append(startpt)
  else:
    for n in range(len(nodes)-1):
      cs, startpt = get_cross_sec(nodes[n], nodes[n+1], stack)
      cross_secs.append(cs)
      startpts.append(startpt)
  # Should have all (N) cross-sections
  return cross_secs, startpts
  
  


