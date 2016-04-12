# neuron_Subtrees.py -- Subtree analysis was getting quite long
# and involved, especially with plotting, so it has been relocated here.


import collections # For namedtuple



def dist3(pt0, pt1):
  if len(pt0) == len(pt1) and len(pt0) == 3:
    return math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
  else:
    print('dimension mismatch')
    print(pt0, pt1)





######################################################################    
# Sub-tree analysis

# Note that sub-tree analysis works in tandem with Matlab for easier
#   3-D plotting and is somewhat interactive (it requires actions from
#   the user). The necessary Matlab scripts are ShollColor and ShollAxons.


def sholl_color(geo, outfile, interdist=1.):
  """
  This color-codes hoc points based on sholl distances.
  Pass interdist = np.inf to skip interpolation (large, detailed files)
  """
  pDF = PathDistanceFinder(geo, geo.soma)
  tips, tipPositions = geo.getTips()
  paths, _ = path_lengths2(geo)
  maxp = max(paths)
  nodes, dists = [], []
  for s in geo.segments:
    if s.length > interdist:
      x0, y0, z0 = s.coordAt(0)[0],s.coordAt(0)[1], s.coordAt(0)[2]
      x1, y1, z1 = s.coordAt(1)[0],s.coordAt(1)[1], s.coordAt(1)[2]
      temp_nodes = [np.linspace(x0, x1, int(s.length/interdist)),
                    np.linspace(y0, y1, int(s.length/interdist)),
                    np.linspace(z0, z1, int(s.length/interdist))]
      #print(np.shape(temp_nodes))
      for t in range(np.shape(temp_nodes)[1]):
        nodes.append([temp_nodes[0][t],temp_nodes[1][t],temp_nodes[2][t]])
        dists.append(pDF.distanceTo(s, 0.5)/maxp)
    else: # Add at least one node per segment
      nodes.append(s.coordAt(0.5))
      dists.append(pDF.distanceTo(s, 0.5)/maxp)
  if len(dists) != len(nodes):
    print('Warning! No. nodes: %i, No. distances: %i' %(len(nodes), len(nodes)))
  with open(outfile, 'w') as fOut:
    for i in range(len(nodes)):
      fOut.write('%.5f %.5f %.5f %.5f' %(nodes[i][0], nodes[i][1], 
                                         nodes[i][2], dists[i]))
      fOut.write('\n')
  print('%s file written' %outfile)
  return
  

def axons_endpoints(geo, outfile=None, Format='matlab'):
  # This prints possible axons and their start and end points. Feed
  # the output into shollAxons.m to select only the "true" axon(s).
  axons = []
  for s in geo.branches: # Operates on BRANCHES!
    if "Axon" in s.tags or "axon" in s.tags:
      axons.append([geo.branches.index(s), s.coordAt(0), s.coordAt(1)])
  if outfile is not None:
    print('Found %i potential branch axons' %len(axons))
    with open(outfile, 'w') as fOut:
      for a in axons:
        fOut.write('%i %.5f %.5f %.5f\n%i %.5f %.5f %.5f' 
                   %(a[0], a[1][0], a[1][1], a[1][2],
                     a[0], a[2][0], a[2][1], a[2][2]))
    return
  else:
    for a in axons:
      if Format == 'matlab':
        print('%i %.5f %.5f %.5f\n%i %.5f %.5f %.5f' 
              %(a[0], a[1][0], a[1][1], a[1][2],
                a[0], a[2][0], a[2][1], a[2][2]))
      else:
        print(a[0], a[1], a[2])
    print('Found %i potential branch axons' %len(axons))
    return a



def axon_help(geo, x=None, y=None, z=None):
  for s in geo.segments:
    for g in [0,1]:
      if x is not None:
        if int(s.coordAt(g)[0]) == x:
          if y is not None:
            if int(s.coordAt(g)[1]) == y:
              print(geo.segments.index(s), s.coordAt(g))
          else:
            print(geo.segments.index(s), s.coordAt(g))
  return



def simple_axon(geo, axon, thing='branch'):
  # Turn the axon into a NeuronGeometry.Segment (but NOT a branch!)
  def mid_pt(geo):
    midPt = [np.mean([s.coordAt(0)[0] for s in geo.segments]),
             np.mean([s.coordAt(0)[1] for s in geo.segments]),
             np.mean([s.coordAt(0)[2] for s in geo.segments])]
    return midPt
  #
  def distal_seg(geo, axon, midPt): # If axon==branch, return the most distal segment
    dists = [dist3([n.x,n.y,n.z], midPt) for n in axon.nodes]
    node = [x for (y,x) in sorted(zip(dists, axon.nodes))][-1]
    for s in geo.segments:
      if node in s.nodes:
        return s
  #
  midPt = mid_pt(geo)
  if thing == 'branch':
    return distal_seg(geo, geo.branches[axon], midPt)
  elif thing == 'seg' or thing == 'segment':
    return geo.segments[axon]



def axon_path(geo, axons=None, things=None, outfile=None, interdist=1.):
  """
  Given a geofile and an axon (default axon if none provided), this
  returns the subtrees between the soma and that axon(s).
  """
  print(axons, things)
  if axons is not None and things is not None:
    if type(axons) is not list:
      axons = [axons]
    if things is not None:
      if type(things) is not list:
        things = [things]
    if len(axons) != len(things):
      if len(things) == 1:
        things = [things[0] for t in axons]
      else:
        print('_Things_ should be len(axons) or length 1 (if all are same type)')
        return
    beh = [simple_axon(geo, a, t) for a, t in zip(axons, things)]
  elif axons is not None and things is None:
    beh = [a for a in axons]
  else:
    print(axons, things)
    for ax in beh:
      print(ax.filamentIndex)
    return
  # Show the filamentIndex of the axons (for reference)
  print('found %i axons' %len(beh))
  for ax in beh:
    print(ax.filamentIndex)
  # Now get the main paths
  pDF = PathDistanceFinder(geo, geo.soma)
  paths = [pDF.pathTo(b) for b in beh]
  # If no outfile, just return the path
  if outfile is None:
    return paths
  else: # Else, write it to a file (with interpolation as before)
    nodes, dists = [], []
    for path in paths:
      for p in path:
        if p.length > interdist:
          x0, y0, z0 = p.coordAt(0)[0],p.coordAt(0)[1], p.coordAt(0)[2]
          x1, y1, z1 = p.coordAt(1)[0],p.coordAt(1)[1], p.coordAt(1)[2]
          temp_nodes = [np.linspace(x0, x1, int(p.length/interdist)),
                        np.linspace(y0, y1, int(p.length/interdist)),
                        np.linspace(z0, z1, int(p.length/interdist))]
          #print(np.shape(temp_nodes))
          for t in range(np.shape(temp_nodes)[1]):
            nodes.append([temp_nodes[0][t],temp_nodes[1][t],temp_nodes[2][t]])
            #dists.append(pDF.distanceTo(s, 0.5)/maxp)
        else: # Add at least one node per segment
          nodes.append(p.coordAt(0.5))
          #dists.append(pDF.distanceTo(s, 0.5)/maxp)
    with open(outfile, 'w') as fOut:
      for i in range(len(nodes)):
        fOut.write('%.6f %.6f %.6f' %(nodes[i][0], nodes[i][1], 
                                           nodes[i][2]))#, dists[i]))
        fOut.write('\n')
    print('%s file written' %outfile)
  return



def mainpath_color(geo, outfile, axon=None, things=None, interdist=1.):
  # This will color segments based on their distance from the mainpath
  print(axon)
  paths = np.array(axon_path(geo, axon, things))
  bpaths = []
  for path in paths:
    for p in path:
      bpaths.append(p)
  downpaths = bpaths[::20]
  segdists = []
  count = -1
  for s in geo.segments:
    count = count + 1
    if s in paths:
      segdists.append(0)
    else:
      pDF = PathDistanceFinder(geo, s)
      segdists.append(min([pDF.distanceTo(p) for p in downpaths]))
    if count % 10 == 0:
      print('%i (of %i) segments processed' %(count, len(geo.segments)))
  if len(segdists) != len(geo.segments):
    print('Num segdists (%i) not equal to num segments (%i)'
          %(len(segdists), len(geo.segments)))
    return
  pDF = PathDistanceFinder(geo, geo.soma)
  maxp = max(segdists)
  nodes, dists = [], []
  for s in geo.segments:
    if s.length > interdist:
      x0, y0, z0 = s.coordAt(0)[0],s.coordAt(0)[1], s.coordAt(0)[2]
      x1, y1, z1 = s.coordAt(1)[0],s.coordAt(1)[1], s.coordAt(1)[2]
      temp_nodes = [np.linspace(x0, x1, int(s.length/interdist)),
                    np.linspace(y0, y1, int(s.length/interdist)),
                    np.linspace(z0, z1, int(s.length/interdist))]
      #print(np.shape(temp_nodes))
      for t in range(np.shape(temp_nodes)[1]):
        nodes.append([temp_nodes[0][t],temp_nodes[1][t],temp_nodes[2][t]])
        dists.append(segdists[geo.segments.index(s)])
    else: # Add at least one node per segment
      nodes.append(s.coordAt(0.5))
      dists.append(segdists[geo.segments.index(s)])
  with open(outfile, 'w') as fOut:
    for n in range(len(nodes)):
      fOut.write('%.6f %.6f %.6f %.6f\n' %(nodes[n][0], nodes[n][1], 
                                         nodes[n][2], dists[n]))
  print('File %s written with %i nodes' %(outfile, len(nodes)))
  return




# Helper function for get_subtrees
def add_unidirect(geo, seg, prevseg, path):
  # Add all neighbors in a unidirectional manner.
  subtree = [] # This is a collection of segments in the subtree
               # There is no inherent structure to the subtree
  for n in seg.neighbors:
    if n != seg and n != prevseg and n not in path:
      subtree.append(n)
  same = 0
  while same < 5:
    prevlength = len(subtree)
    for s in subtree:
      for n in s.neighbors:
        if n != seg and n != prevseg and n not in subtree and n not in path:
          subtree.append(n)
    if len(subtree) == prevlength:
      same = same + 1
    else:
      same = 0
      print('subtree is length: %i' %len(subtree))
  return subtree



def get_subtrees(geo, axons, things=None):
  # This returns the subtrees coming off the 'main' path (soma->axon)
  # Should work with multiple axons.
  def is_subtree_root(geo, seg, path):
    subtree = []
    for n in seg.neighbors:
      if n not in np.array(path).flatten():
        subtree.append(add_unidirect(geo, n, seg, path)) ## 
    return subtree
  
  # Condition axon input
  paths = axon_path(geo, axons, things)
  flatpath = []
  for path in paths:
    for p in path:
      if p not in flatpath:
        flatpath.append(p)
  # Assume path[0] == geo.soma
  
  subtrees = []
  for p in flatpath:
    # Loop over the segments in the path to create the subtrees
    if len(p.neighbors) > 2:
      subs = is_subtree_root(geo, p, flatpath)
      for sub in subs:
        if sub not in subtrees:
          if len(sub) not in [len(i) for i in subtrees]:
            subtrees.append(sub)
  # Should have all unique subtrees now
  return subtrees



def show_subtrees(geo, outfile, axon=None, things=None, subtrees=None, interdist=1.):
  # Make a txt file to plot individual subtrees in matlab
  if subtrees is None:
    subtrees = get_subtrees(geo, axon, things)
    path = axon_path(geo, axon, things)
    subtrees.append(path)
  treenum = np.linspace(0,1,len(subtrees))
  # Interpolate for plotting
  nodes, dists = [], []
  for tree in range(len(subtrees)):
    for seg in subtrees[tree]:
      if type(seg) is list:
        for s in seg:
          if s.length > interdist:
            x0, y0, z0 = s.coordAt(0)[0],s.coordAt(0)[1], s.coordAt(0)[2]
            x1, y1, z1 = s.coordAt(1)[0],s.coordAt(1)[1], s.coordAt(1)[2]
            temp_nodes = [np.linspace(x0, x1, int(s.length/interdist)),
                          np.linspace(y0, y1, int(s.length/interdist)),
                          np.linspace(z0, z1, int(s.length/interdist))]
            #print(np.shape(temp_nodes))
            for t in range(np.shape(temp_nodes)[1]):
              nodes.append([temp_nodes[0][t],temp_nodes[1][t],temp_nodes[2][t]])
              dists.append(treenum[tree])
          else: # Add at least one node per segment
            nodes.append(s.coordAt(0.5))
            dists.append(treenum[tree])
      else:
        s = seg
        if s.length > interdist:
          x0, y0, z0 = s.coordAt(0)[0],s.coordAt(0)[1], s.coordAt(0)[2]
          x1, y1, z1 = s.coordAt(1)[0],s.coordAt(1)[1], s.coordAt(1)[2]
          temp_nodes = [np.linspace(x0, x1, int(s.length/interdist)),
                        np.linspace(y0, y1, int(s.length/interdist)),
                        np.linspace(z0, z1, int(s.length/interdist))]
          #print(np.shape(temp_nodes))
          for t in range(np.shape(temp_nodes)[1]):
            nodes.append([temp_nodes[0][t],temp_nodes[1][t],temp_nodes[2][t]])
            dists.append(treenum[tree])
        else: # Add at least one node per segment
          nodes.append(s.coordAt(0.5))
          dists.append(treenum[tree])
  with open(outfile, 'w') as fOut:
    for n in range(len(nodes)):
      fOut.write('%.6f %.6f %.6f %.6f\n' %(nodes[n][0], nodes[n][1], 
                                         nodes[n][2], dists[n]))
  print('File %s written with %i nodes' %(outfile, len(nodes)))
  return



### Wiring
def subtree_wiring(geo, subtrees=None, path=None, axons=None, things=None):
  # Returns the percent of total wiring of each subtree (and the path)
  if subtrees is None:
    ssubtrees = get_subtrees(geo, axons, things)
  else:
    ssubtrees = subtrees
  # Get rid of zero subtree
  pico = [len(t) for t in ssubtrees]
  if min(pico) == 0:
    ssubtrees.pop(pico.index(0))
  if path is None:
    paths = axon_path(geo, axons, things)
  else:
    paths = path
  pathwire = sum([sum([p.length for p in pa]) for pa in paths])
  allsegs = []
  for t in ssubtrees:
    for s in t: # Add each segment
      allsegs.append(s)
  total = sum([s.length for s in np.array(allsegs)])
  total = total + pathwire
  for p in paths:
    ssubtrees.append(p)
  subwiring = []
  for t in ssubtrees:
    subwiring.append(sum([s.length for s in t])/total)
  print('Subtree length: %i' %len(ssubtrees))
  print('Max: %.5f, Min: %.5f' %(max(subwiring), min(subwiring)))
  print('Mean: %.5f, Med: %.5f' %(np.mean(subwiring), np.median(subwiring)))
  print('Path wiring: %.5f, Path wiring percent: %.5f' %(pathwire, 
                                                         pathwire/total))
  print('Total tree wiring: %.5f' %total)
  print('Axons: '); print(axons)


def summary_subtree(prop_dict, types=['GM']):
  """
  Report summary statistics on subtrees for desired cell types.
  """
  for t in types:
    for u in range(len(prop_dict['cellTypes'])):
      if prop_dict['cellTypes'][u] == t:
        
        print(prop_dict['files'][u])
        print('Max: %.5f, Min: %.5f' %(max(prop_dict['subtree_wires'][u]),
                                       min(prop_dict['subtree_wires'][u])))
        print('Mean: %.5f, Med: %.5f' %(np.mean((prop_dict['subtree_wires'][u])),
                                        np.median((prop_dict['subtree_wires'][u]))))
        try:
          print('Path wiring: %.5f, Path wiring percent: %.5f'
                %(sum(prop_dict['subtree_pathlengths'][u]),
                  sum(prop_dict['subtree_pathlengths'][u])/
                  sum(prop_dict['subtree_wires'][u])))
        except:
          print('Path wiring: %.5f, Path wiring percent: %.5f'
                %(prop_dict['subtree_pathlengths'][u],
                  prop_dict['subtree_pathlengths'][u]/
                  sum(prop_dict['subtree_wires'][u])))
        print('%i subtrees' %len(prop_dict['subtree_wires'][u]))
  return

         
  



### Get filamentInds for subtrees
def subtree_filaments(geo, axons, things, retsegs=False):
  """
  This allows for rapid reconstruction of subtrees with subtree_fromfilaments.
  """
  subtrees = get_subtrees(geo, axons, things)
  # Assume first segment in subtree is subtree root
  # Get the filamentIndex for all segments and save by subtree
  subtrees = [s for s in subtrees if len(s) > 0]
  sub_dict = {s[0].filamentIndex: [seg.filamentIndex for seg in s]
              for s in subtrees}
  
  # return the NeuronGeometry segments?
  if retsegs:
    return sub_dict, {s[0].filamentIndex: s for s in subtrees}
  return sub_dict



def subtree_fromfilaments(geo, fil_dict, asdict=True, version=1):
  """
  Given a geo object and the dictionary of subtree roots and all the
  segments contained in that subtree, this populates the segment dict.
  """
  if version == 0:
    seg_dict = {seg:[] for seg in geo.segments if seg.filamentIndex in
                fil_dict.keys()}
    for k in seg_dict.keys():
      seg_dict[k] = [seg for seg in geo.segments if seg.filamentIndex in
                     fil_dict[k.filamentIndex]]
    
    if asdict:
      return seg_dict
    return seg_dict.values()

  elif version == 1:  
    # Alternatively...
    seg_dict = {}
    for k in fil_dict.keys():
      this_key = [seg for seg in geo.segments if seg.filamentIndex == int(k)][0]
      #print(this_key)
      seg_dict[this_key] = [seg for seg in geo.segments if seg.filamentIndex
                            in [int(u) for u in fil_dict[str(this_key.filamentIndex)]]]
      print('Found %i segments for %i' %(len(seg_dict[this_key]), 
                                         this_key.filamentIndex))
    if asdict:
      return seg_dict
    return list(seg_dict.values())
  



def subtree_segments(geo, subtree_inds):
  """
  Return the subtrees as NeuronGeometry segments instead of filaments.
  These subtrees are location dicts. (NO. USE SUBTREE_TIPS)
  """
  subtrees = {k: [seg for seg in geo.segments if seg.filamentIndex in
                  subtree_inds[k]] for k in subtree_inds.keys()}
  return subtrees
  


def subtree_tips(geo, subtrees, asindex=False):
  """
  Given the subtrees (either asindex=True (numbers only) or False 
  (NeuronGeometry segment objects), get the tip clusters.
  """
  if asindex:
    subtrees = subtree_segments(geo, subtrees)
  
  # Get the subtree start location and the tip locations
  pDF = PathDistanceFinder(geo, geo.soma)
  tipInds, tipLocs = geo.getTipIndices()
  
  # A dict of dicts
  locations = {k: {'root': [seg.coordAt(1) for seg in geo.segments if
                            seg.filamentIndex == int(k)],
                   'tips': [seg.coordAt(tipLocs[tipInds.index(seg.filamentIndex)])
                            for seg in subtrees[k] if seg.filamentIndex
                            in tipInds]} for k in subtrees.keys()}
  
  # Quality testing
  for k in list(locations.keys()):
    if len(locations[k]['root']) == 0: # If no root was found, check all segments
      locations[k]['root'] = [seg.coordAt(1) for seg in geo.segments if
                              seg.filamentIndex == k][0]
    elif len(locations[k]['root']) > 0: # Otherwise, use the first (only) root
      locations[k]['root'] = locations[k]['root'][0]
    else:
      print(locations[k]['root']) # Something strange happened
  
  return locations



def randomize_subtrees(locations, N=2000, keep=0):
  """
  Assign tip clusters to different subtrees; repeat N times and send to
  subtree_statistic for each tree.
  """
  from scipy.special import comb
  # Condition the tree -- make each tip clustered mean-centered
  import copy
  tiplocs, rootlocs = [], []
  for k in locations.keys():
    for tip in locations[k]['tips']:
      tiplocs.append([tip[0]-locations[k]['root'][0],
                      tip[1]-locations[k]['root'][1],
                      tip[2]-locations[k]['root'][2]])
    rootlocs.append([locations[k]['root'], len(locations[k]['tips'])]) # Location of tip AND how many it gets
  
  # An original root with 7 tips can get ANY 7 new random tips, but only 7
  
  # If there are >5 times fewer possible options than requested trials
  numcombs = comb(len(rootlocs), len(rootlocs), repetition=False)
  if N/5. > numcombs:
    print('Requested %i trials but only %.2f combinations are possible'
          %(N, float(numcombs)))
    N = int(numcombs*5)
  
  # Run N times
  stat, ret, retain = [], [int(i) for i in np.random.rand(keep)], []
  for n in range(N):
    if n%200==0:
      print('%i / %i subtrees constructed so far' %(n, N))
    
    # Generate the dictionary
    trialtips = [int(i) for i in np.random.random(len(tiplocs))*len(tiplocs)] # Num tips
    trialroots = [int(i) for i in np.random.random(len(rootlocs))*len(rootlocs)] # Num roots
    # Make the randomized test dictionary
    this_dict = {}
    for t in range(len(trialroots)): # For each root (NOT each tip)
      this_dict[t] = {'root': rootlocs[t][0],
                      'tips': [[i[0]+rootlocs[t][0][0], # Already root-mean-centered,
                               i[1]+rootlocs[t][0][1], # So add it back
                               i[2]+rootlocs[t][0][2]] for i in 
                               # Only use the next tips, however many this cluster should have
                               [tiplocs[o] for o in trialtips[0:rootlocs[t][1]]]]}
                               #tiplocs[trialtips[0]:trialtips[rootlocs[t][1]]]]}
      
      # Pop the used nodes from trial tips
      for u in range(rootlocs[t][1]):
        _ = trialtips.pop(0)
    if n in ret:
      retain.append(this_dict)
    
    # Get the test statistic for this new dict
    try:
      stat.append(subtree_statistics(this_dict))
    except:
      pass
  if keep > 0:
    return stat, retain
  return stat



def subtree_statistics(testdict, cluster_rad=False, intercluster_dist=False,
                       cluster_overlap=True):
  """
  Given a subtree of format {ind: {'root': (xyz), 'tips': [(xyz,xyz...)]}}
  This calculates the avg distance from the center of mass (of each cluster)
  to its tips (cluster_rad) and the avg distance from this cluster's center 
  of mass to the next-closest cluster's center of mass. (ind is meaningless)
  """
  clustrads, clustcenter, stds = [], [], []
  
  # Find cluster center (of mass) and avg 'radius' of each cluster
  for k in testdict.keys():
    if len(testdict[k]['tips']) > 0:
      clustcenter.append([ np.mean([x[0] for x in testdict[k]['tips']]),
                           np.mean([x[1] for x in testdict[k]['tips']]),
                           np.mean([x[2] for x in testdict[k]['tips']]) ])
      
      clustrads.append( np.mean([dist3(clustcenter[-1], tip) for tip in
                                 testdict[k]['tips']]) )
  if cluster_rad:
    return np.mean(clustrads)
    
  # Find distance to closest cluster
  elif intercluster_dist:
    interclust = [min([dist3(i, u) for u in clustcenter if
                       u != i]) for i in clustcenter]
    return np.mean(interclust)
  #print(clustcenter)
  
  # Overlap of standard deviation
   # Each entry is [xmin, xmax, ymin, ymax]
  elif cluster_overlap:
    Bounds = collections.namedtuple('Bounds', 'xmin xmax ymin ymax') # Named tuple
    for k in testdict.keys():
      if len(testdict[k]['tips']) > 0:
        mn = [np.mean([x[0] for x in testdict[k]['tips']]),
              np.mean([x[1] for x in testdict[k]['tips']])]
        st = [np.std([x[0] for x in testdict[k]['tips']]),
              np.std([x[1] for x in testdict[k]['tips']])]
        stds.append(Bounds(xmin=mn[0]-st[0], xmax=mn[0]+st[0], 
                           ymin=mn[1]-st[1], ymax=mn[1]+st[1]))
    
    # Detect overlaps
    overlaps = [[] for i in stds]
    for ref in stds:
      for st in stds:
        if ref != st: # Don't check the same cluster
          if st.xmin < ref.xmax and st.xmax > ref.xmin: # X lower overlap
            overlaps[stds.index(ref)].append(st)
          if st.xmax > ref.xmin and st.xmin < ref.xmax: # X upper overlap
            overlaps[stds.index(ref)].append(st)
          if st.ymin < ref.ymax and st.ymax > ref.ymin: # Y lower overlap
            overlaps[stds.index(ref)].append(st)
          if st.ymax > ref.ymin and st.ymin < ref.ymax: # Y upper overlap
            overlaps[stds.index(ref)].append(st)
    
    # Prune copies
    # print(overlaps[0])
    overlaps = [list(set(u)) for u in overlaps]
    overlaps = [len(u) for u in overlaps]
    return np.mean(overlaps)
  
  else:
    print('must set ONLY ONE to true at a time: cluster_rad, intercluster_dist,cluster_overlap')
  return None





def interclust_dist_test():
  """
  See if interclust distance is a valuable metric.
  """
  def gen_cloud(cent=(0,0), std=2):
    x, y = np.random.random(20)*std+cent[0], np.random.random(20)*std+cent[1]
    return list(zip(x, y))
  
  # Generate the random clouds
  def randomize_clouds(trueclouds, numclouds=2):
    centers = [[np.mean([pt[0] for pt in cloud]), np.mean([pt[1] for pt in cloud])]
               for cloud in trueclouds]
    all_pts = [[i[0]-centers[0][0], i[1]-centers[0][1]] for i in c1]
    for pt in c2:
      all_pts.append([pt[0]-centers[1][0], pt[1]-centers[1][1]]) # Get the list of all points
    clouds = []
    for u in range(numclouds): # Generate the clouds from random data
      rlist = [int(i) for i in np.random.random(int(sum([len(y) 
               for y in trueclouds]))/numclouds)*sum([len(y) for y in trueclouds])]
      rlist = rlist # Downsample to number of pts for this cloud [::numclouds]
      tpts = [[all_pts[r][0]+centers[u][0], all_pts[r][1]+centers[u][1]] for r in rlist] # Add the center back
      clouds.append(tpts)
    return clouds
  
  # Show the clouds
  def show_clouds(realclouds, randclouds): # Pts are zipped within cloudlists
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1) # First plot original clouds
    rcol = []
    for u in range(len(realclouds)): # Plot each cloud
      rcol.append(np.random.random(3))
      for pt in realclouds[u]:
        ax1.plot(pt[0], pt[1],'o', color=rcol[u], alpha=0.5)
    ax1.set_title('True clouds')
    ax2 = fig.add_subplot(1,2,2) # Then plot random clouds
    for u in range(len(randclouds)):
      for pt in randclouds[u]:
        ax2.plot(pt[0], pt[1], 'o', color=rcol[u], alpha=0.5)
    ax2.set_title('Random clouds')
    plt.show(); return
  
  # Get interclust distance
  def ic_dist(clouds):
    #print(clouds)
    cloud_cents = [ [np.mean([pt[0] for pt in cloud]),
                     np.mean([pt[1] for pt in cloud])] for cloud in clouds ]
    dists = [min([sum([np.sqrt(a[i]**2+b[i]**2) for i in range(2)])
                  for b in cloud_cents if b != a]) for a in cloud_cents]
    return dists
  
  # Run it through
  c1, c2 = gen_cloud((0,0), 3), gen_cloud((5,0), 2) # "real" clouds
  r_clouds = randomize_clouds([c1,c2], 2) # Generate random clouds
  #print(r_clouds)
  print('True IC distances:')
  print(ic_dist([c1,c2]))
  print('Randomized IC distances:')
  print(ic_dist(r_clouds))
  show_clouds([c1,c2], r_clouds)
  return



  


def show_tip_clusters(tipdict, showtips=False, showroots=True):
  """
  Tips and roots can be shown as scatter3d pts with the sd projected
  onto z=0 plane (not optional). Colored by subtree identity.
  """
  fig = plt.figure(figsize=[5,5], dpi=150)
  ax = fig.add_subplot(111, projection='3d')
  from matplotlib.patches import Ellipse
  # For each subtree, plot its root/tips and its projection
  for k in tipdict.keys():
    thic = np.random.rand(3) # thic = random color
    if showtips: # Show the tip locations
      for t in tipdict[k]['tips']:
        ax.scatter(t[0], t[1], t[2], color=thic, edgecolor=thic, s=10,
                   alpha=0.2, )
    if showroots: # Show the root location
      ax.scatter(tipdict[k]['root'][0], tipdict[k]['root'][1], 
                 tipdict[k]['root'][2], color=thic, #edgecolor='black', 
                 s=30, alpha=0.8)
    # Now do ellipses -- only in X-Y (z values ignored)
    dat1 = [t[0] for t in tipdict[k]['tips']]
    dat2 = [t[1] for t in tipdict[k]['tips']]
    cov = np.cov(dat1, dat2)
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    ell = Ellipse(xy=(np.mean(dat1), np.mean(dat2)), # Initial offset (center)
                  width=lambda_[0]*2, height=lambda_[1]*2,
                  angle=np.rad2deg(np.arccos(v[0,0])), alpha=0.1)
    ell.set_facecolor(thic)
    ell.set_edgecolor(thic)
    ell.set_linewidth(1)
    # ell.set_linestyle('--')
    ax.add_patch(ell)
    pathpatch_2d_to_3d(ell, z=0, normal='z')
    pathpatch_translate(ell, (0, 0, 0)) # Here can add _additional_ offset in x,y,z
  # Plot cosmetics
  ax.grid(False) # Remove the grids
  ax.w_xaxis.set_pane_color((0.,0.,0.,0.1))
  ax.w_yaxis.set_pane_color((0.,0.,0.,0.15))
  ax.w_zaxis.set_pane_color((1.,1.,1.,0.))
  ax.w_xaxis.line.set_lw(.2)
  ax.w_yaxis.line.set_lw(.2)
  ax.w_zaxis.line.set_lw(0.)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_zticks([])
  plt.show()
  return



def all_tip_clusters(locDicts, labelsin, switch=False, plotdims=[4,4],
                     showtips=False, showroots=True):
  """
  Plot the roots + tip cluster SD (not individual tips) for each cell.
  Pass as location dictionaries {'root': __, 'tips': {..}}
  """
  if switch:
    for i in range(len(locDicts)-1): # locDicts is actually a list of dicts
      locDicts.append(locDicts.pop(0))
      labelsin.append(labelsin.pop(0))
  colors = ['darkkhaki', 'royalblue', 'forestgreen','tomato']
  # colors = ['forestgreen','tomato', 'darkkhaki', 'royalblue', ] ### !!!!!
  L = list(np.unique(labelsin))
  C = [L.index(i) for i in labelsin]
  # Plot time
  fig = plt.figure(figsize=[5,5], dpi=150)
  from matplotlib.patches import Ellipse
  # For each subtree, plot its root/tips and its projection
  for p in range(len(locDicts)):
    ax = fig.add_subplot(plotdims[0], plotdims[1], p, projection='3d')
    tipdict, cnt = locDicts[p], 0
    thics = colormaker(C[p], len(tipdict.keys())) # Get colors
    for k in tipdict.keys():
      if showtips: # Show the tip locations
        for t in tipdict[k]['tips']:
          ax.scatter(t[0], t[1], t[2], color=thics[cnt], edgecolor=thics[cnt], 
                     s=10, alpha=0.2, fill='none')
      if showroots: # Show the root location
        ax.scatter(tipdict[k]['root'][0], tipdict[k]['root'][1], 
                   tipdict[k]['root'][2], color=thics[cnt], #edgecolor='black', 
                   s=10, alpha=0.8)
      # Now do ellipses -- only in X-Y (z values ignored)
      dat1 = [t[0] for t in tipdict[k]['tips']]
      dat2 = [t[1] for t in tipdict[k]['tips']]
      cov = np.cov(dat1, dat2)
      try:
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)
        ell = Ellipse(xy=(np.mean(dat1), np.mean(dat2)), # Initial offset
                      width=lambda_[0]*2, height=lambda_[1]*2,
                      angle=np.rad2deg(np.arccos(v[0,0])), alpha=0.1)
        ell.set_facecolor(thics[cnt])
        ell.set_edgecolor(thics[cnt])
        ell.set_linewidth(1)
        # ell.set_linestyle('--')
        ax.add_patch(ell)
        pathpatch_2d_to_3d(ell, z=0, normal='z')
        pathpatch_translate(ell, (0, 0, 0)) # Here can add _additional_ offset in x,y,z
      except:
        print('Subtree with only %i tips, no cov!' %len(tipdict[k]['tips']))
        print(cov)
      cnt = cnt + 1
      # Plot cosmetics
      ax.grid(False) # Remove the grids
      ax.w_xaxis.set_pane_color((0.,0.,0.,0.1))
      ax.w_yaxis.set_pane_color((0.,0.,0.,0.15))
      ax.w_zaxis.set_pane_color((1.,1.,1.,0.))
      ax.w_xaxis.line.set_lw(.2)
      ax.w_yaxis.line.set_lw(.2)
      ax.w_zaxis.line.set_lw(0.)
      ax.set_zlim([0, 100])
      ax.set_xticks([])
      ax.set_yticks([])
      ax.set_zticks([])
  plt.show()

def colormaker(whichColor, N=40, amp=True):
  """
  Creates a colormap of N colors.
  """
  wcolor = {0: 'darkkhaki', 1: 'royalblue', 2: 'forestgreen', 3: 'tomato'}
  colorDict = {'darkkhaki': [[250, 161], [217, 140], [85, 59]],
               'royalblue': [[179, 17], [186, 43], [245, 242]],
               'forestgreen': [[37, 68], [120, 245], [32, 59]],
               'tomato': [[217, 217], [90, 151], [0, 104]]}
  these = colorDict[wcolor[whichColor]]
  for t in range(len(these)):
    these[t] = np.linspace(these[t][0]/255., these[t][1]/255., N)
  if amp:
    these.append(np.random.rand(N))
  newt = [[these[0][i], these[1][i], these[2][i], these[3][i]] for i in 
          range(len(these[0]))]
  return newt




########################### Ellipse helpers! ################################

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import art3d
from mpl_toolkits.mplot3d import proj3d
import numpy as np

def rotation_matrix(v1,v2):
  """
  Calculates the rotation matrix that changes v1 into v2.
  """
  v1/=np.linalg.norm(v1)
  v2/=np.linalg.norm(v2)
  cos_angle=np.dot(v1,v2)
  d=np.cross(v1,v2)
  sin_angle=np.linalg.norm(d)
  if sin_angle == 0:
    M = np.identity(3) if cos_angle>0. else -np.identity(3)
  else:
    d/=sin_angle
    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)
    M = ddt + cos_angle * (eye - ddt) + sin_angle * skew
  return M


def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
  """
  Transforms a 2D Patch to a 3D patch using the given normal vector.
  The patch is projected into they XY plane, rotated about the origin
  and finally translated by z.
  """
  if type(normal) is str: #Translate strings to normal vectors
    index = "xyz".index(normal)
    normal = np.roll((1,0,0), index)
  path = pathpatch.get_path() #Get the path and the associated transform
  trans = pathpatch.get_patch_transform()
  path = trans.transform_path(path) #Apply the transform
  pathpatch.__class__ = art3d.PathPatch3D #Change the class
  pathpatch._code3d = path.codes #Copy the codes
  pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color    
  verts = path.vertices #Get the vertices in 2D
  M = rotation_matrix(normal,(0, 0, 1)) #Get the rotation matrix
  pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])


def pathpatch_translate(pathpatch, delta):
  """
  Translates the 3D pathpatch by the amount delta.
  """
  pathpatch._segment3d += delta


def eigsorted(cov):
  """
  Necessary for the 2d ellipses created in pretty_2d scatter plot.
  """
  vals, vecs = np.linalg.eigh(cov)
  order = vals.argsort()[::-1]
  return vals[order], vecs[:,order]


