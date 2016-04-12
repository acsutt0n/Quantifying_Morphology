# usage: python neuron_Asymmetry.py hocFile
# calculates partition asymmetry by magic

import sys
from neuron_readExportedGeometry import *



def get_segment(geo, segname):
  for s in geo.segments:
    if s.name == segname:
      return s



def tips_asymmetry(geo):
  # Get the tip asymmetry of the neuron. Follow the soma's neighbors
  # until there are more than 1, then start there.
  prevsegs = [geo.soma]
  newsegs = [i for i in geo.soma.neighbors if i not in prevsegs]
  go = True
  while go:
    if len(newsegs) > 1:
      nebs = newsegs
      go = False
    else:
      for k in newsegs:
        prevsegs.append(k)
        for j in k.neighbors:
          newsegs.append(j)
        # not sure if this is allowed, but should be since not referencing by index
        newsegs.pop(newsegs.index(k))
        
  pDF = PathDistanceFinder(geo, geo.soma, 0.5)
  # nebs = geo.soma.neighbors
  tips, tipPositions = geo.getTips()
  seg_names = {}
  seg_tips = {}
  for n in nebs:
    seg_names[n.name] = []
    seg_tips[n.name] = []
  seg_lengths = {}
  for t, pos in zip(tips, tipPositions):
    curr_path = pDF.pathTo(t, pos)
    for n in seg_names.keys():
      # if the bifurcation in question is contained in the path soma->tip
      if get_segment(geo,n) in curr_path:
        # add this tip to n
        seg_tips[n].append(t)
        for c in curr_path:
          if c not in seg_names[n]:
            seg_names[n].append(c)
  # now should have all of the segments that lead to the tips in each key
  return seg_names, seg_tips
  


# for dev
geo = demoRead('/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/786_062_63xgly_LG_ACS_10.29.2014.hoc')
seg_names, seg_tips = tips_asymmetry(geo)




def lengths_asymmetry(geo, seg_names):
  # For lengths, fill in the rest of of the segments
  # First, make a list of segments that aren't allowed to be added to
  # each subsequent tree
  from itertools import cycle
  seg_keys = list(seg_names.keys())
  done = [] # not used in this version of the function
  _seg = cycle(seg_keys)
  # need to cycle through until no_good isn't changing anymore?
  iterations = 0
  while iterations < 5:
    print('Iteration # %i' %iterations)
    # reset everything, move to next segment
    curr_seg = next(_seg)
    no_good = []
    # create no_good, a list of segments that can't be added to the 
    # current tree's neighbors
    for k in seg_keys:
      if k != curr_seg:
        for s in seg_names[k]:
          no_good.append(s)
    # add the current tree's segments' neighbors to the current seg tree
    #  if they're not connected to other trees
    # this may 'steal' some neighbors of neighbors of another tree,
    # but shit, you have to start somewhere; also soma needs to be excluded
    for s in seg_names[curr_seg]:
      for n in s.neighbors:
        if n not in no_good:
          seg_names[curr_seg].append(n)
    iterations = iterations + 1
    # add the current segment to done
    # done.append(curr_seg)
  # len(done) == len(seg_keys) -- all the segments are done
  # now, get the lengths
  print('Adding lengths...')
  seg_lengths = {}
  for k in seg_keys:
    seg_lengths[k] = 0
    for l in seg_names[k]:
      seg_lengths[k] = seg_lengths[k] + l.length
  # now should have all the segments in the entire geometry, more or less
  return seg_names, seg_lengths



# for dev
seg_names2, seg_lengths = lengths_asymmetry(geo, seg_names)




















###################################################
if __name__ == '__main__':
  arguments = sys.argv
  if len(arguments) < 2:
    print('Need a hoc file, dumbass!')
  else:
    hocFile = arguments[1]
