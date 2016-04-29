# retrieve simple properties from a geo instance

from neuron_readExportedGeometry import *
import numpy as np
import math, os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import networkx as nx
import seaborn as sns
import csv
import copy



# helper functions
def name(geo):
  return geo.fileName.split('/')[-1].split('.')[0]


def farthest_pt(pts):
  dmax = 0
  for i in pts:
    for j in pts:
      if dist3(i,j) > dmax:
        dmax = dist3(i,j)
  return dmax


def farthest_pt2(pts):
  sumpts = [sum(p) for p in pts]
  minpt = pts.index(sumpts.index(min(sumpts)))
  maxpt = pts.index(sumpts.index(max(sumpts)))
  dmax = dist3(minpt, maxpt)
  return dmax*1.5
    


def checko(obj):
  unique_files, unique_items, unique_cells = None, None, None
  if type(obj) is not dict:
    print('Only works for dictionaries'); return 
  if len(obj['files']) != len(np.unique(obj['files'])):
    print('Duplicates found in files!')
  unique_files = len(np.unique(obj['files']))
  for k in obj.keys():
    if k != 'files' and k != 'cellTypes' and k != 'cellType':
      if len(obj[k]) != len(np.unique(obj[k])):
        print('Duplicates found in %s!' %k)
      unique_items = len(np.unique(obj[k]))
  try:
    unique_cells = len(np.unique(obj['cellTypes']))
  except:
    unique_cells = len(np.unique(obj['cellType']))
  print('Contents: %i unique files, %i unique items, %i cell types'
         %(unique_files, unique_items, unique_cells))
  return


def nodex(n):
  return [n.x, n.y, n.z]


def node3(n0, n1):
  return dist3(nodex(n0),nodex(n1))


def union(list1, list2, rall=False):
  """
  Returns the common segment between two lists. Union assumes there is only
  one common segment. First, it assumes the segment is not included in
  the lists, so it crawls through the neighbors of elements of both lists. If
  it finds more than 1 common segment, it assumes the segment is included
  in one of the lists, and then returns the common segment found by neighbors.
  rall - returns list of common segments (return all)
  """
  def union1(list1, list2): 
    for a in list1:
      for n in a.neighbors:
        if n in list2:
          return a
    return None
  def union2(list1, list2, rall=False):
    u = []
    for a in list1:
      for b in list2:
        for n in a.neighbors:
          if n in b.neighbors:
            if n not in u:
              u.append(n)
    if rall is True and len(u) != 1:
      # print('Found %i matches' %len(u))
      return u
    if len(u) > 1:
      return union1(list1, list2)
    elif len(u) == 1:
      return u[0]
    else:
      return None
  return union2(list1, list2, rall)



#######################################################################
# prepare for analysis -- load up hoc files
# default cellTypes (spreadsheet):
cellTypes =  ['LG', 'LG', 'LG', 'LG', 'LP', 'LP', 'LP', 'LP', 'PD', \
              'PD', 'PD', 'PD', 'GM', 'PD']

def get_hocfiles(directory='/home/alex/data/morphology/morphology-hoc-files/morphology/analyze/'):
  fils = os.listdir(directory)
  hocfiles = [i for i in fils if i.split('.')[-1]=='hoc']
  hocfiles = [directory+i for i in hocfiles]
  return hocfiles
  

def get_geofiles(directory='/home/alex/data/morphology/morphology-hoc-files/hocs/all/'):
  hocfiles = get_hocfiles(directory)
  geofiles = [demoReadsilent(g) for g in hocfiles]
  return geofiles, hocfiles


def match_celltype(names):
  celltype_dict = {'815_144_GM_scaled': 'GM',
                   '785_040_63xgly_new_LG': 'LG',
                   '798_007_63xgly_GM_knossos_scaled': 'GM',
                   '798_002_NoLoops_scaled': 'GM',
                   '836_047_63X_RH_199': 'LP',
                   '803_155_63x_JL_03': 'PD',
                   '836_095_63X_IM_34': 'LP',
                   '803_151_63x_IM_69': 'PD',
                   '791_076_63xgly_knossos_LG_soma-edit_scaled': 'LG',
                   '837_123_63x_IM_79_radius_2': 'LP',
                   '791_088_63xgly_knossos_LG_soma-edit_scaled': 'LG',
                   '786_062_63xgly_LG_ACS_10_scaled': 'LG',
                   '798_032_63x_NA_05': 'PD',
                   '837_103_63x_DK_38': 'LP',
                   '803_133_63x_NA_11': 'PD',
                   '791_127_63x_GM_ACS_soma-edit_scaled': 'GM'}
  return [celltype_dict[h] for h in names]
  





#######################################################################
# Branch angles

def dist3(pt0, pt1):
  if len(pt0) == len(pt1) and len(pt0) == 3:
    return math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
  else:
    print('dimension mismatch')
    print(pt0, pt1)




def get_angle(pt0, midpt, pt1):
  if pt0 in [midpt, pt1] or pt1 in [midpt, pt0] or midpt in [pt0,pt1]:
    print('Some points are the same!')
    print(pt0, midpt, pt1)
  PT0 = dist3(pt1, midpt)
  PT1 = dist3(pt0, midpt)
  MIDPT = dist3(pt0, pt1)
  try:
    ang = math.acos( (MIDPT**2 - PT1**2 - PT0**2) / (2*PT1*PT0) )
    ang = ang*180/math.pi
  except:
    ang = 'nan'
  return ang



def find_points(seg0, seg1):
  seg0list, seg1list = [], []
  pt0where, pt1where, midwhere = None, None, None
  switchdict = {0: -1, -1: 0}
  # make a list of the node locations
  for n in seg0.nodes:
    seg0list.append([n.x,n.y,n.z])
  for n in seg1.nodes:
    seg1list.append([n.x,n.y,n.z])
    # find the common node, then use that to find the distinct ones
  for n in seg0list:
    if n in seg1list:
      midpt = n
  if seg0list.index(midpt) != 0:
    pt0where = 0
    pt0 = seg0list[0]
  else:
    pt0where = -1
    pt0 = seg0list[-1]
  if seg1list.index(midpt) != 0:
    pt1where = 0
    pt1 = seg1list[0]
  else:
    pt1where = -1
    pt1 = seg1list[-1]
  
  f = True
  if pt0 == pt1 or pt0==midpt:
    f = False
    if pt0where == 0:
      try:
        pt0=seg0list[1]
        f = True
      except:
        pass
    elif pt0where == -1:
      try:
        pt0=seglist[-2]
        f = True
      except:
        pass
  if pt0 == pt1 or pt1==midpt:
    if pt1where == 0:
      try:
        pt1=seg1list[1]
        f = True
      except:
        pass
    elif pt1where == -1:
      try:
        pt1=seg1list[-2]
        f = True
      except:
        pass
  if f == False:
    print('Tried to find new coordinates, but failed. Skipping')
  if pt0 in [midpt, pt1] or pt1 in [midpt, pt0] or midpt in [pt0,pt1]:
    print(seg0list, seg1list)
  #print('pt0 at %i, pt1 at %i' %(pt0where, pt1where))
  if pt1 and pt0 and midpt:
    return pt0, midpt, pt1
  else:
    print('could not figure out segments %s and %s' %(seg0.name, seg1.name))
    print(seg0list, seg1list)
    return [False]


def branch_angles(geo, outfile=None, pairedAngles=False):
  f_angles, f_borders = [], []
  geo.calcForewardBranchOrder()
  for b in geo.branches:
    for n in b.neighbors:
      pts = find_points(n, b)
      if len(pts) == 3:
        pt0, midpt, pt1 = pts[0], pts[1], pts[2]
      f_angles.append(get_angle(pt0, midpt, pt1))
      f_borders.append(b.branchOrder)
  angles, borders = [], []
  # Purge the borders and angles of 'nan'
  for a in range(len(f_angles)):
    if str(f_angles[a]) == 'nan' or f_angles[a] == 'nan':
      pass
    else:
      angles.append(f_angles[a])
      borders.append(f_borders[a])
  if outfile is not None:
    with open('temp_angles.txt', 'w') as fOut:
      for a in angles:
        fOut.write('%.10f, \n' %a)
  if pairedAngles is False:
    return angles
  else:
    return angles, borders



##########################################################################
# Branch degree
# This gives two measures: (1) treating each branch as a node in
#   graph theory terms, gives the degree (# of total neighbors) of that
#   node; (2) how many branches arise from each branch point (not
#   including parent branch), so a traditional bifurcation would be 2.

def branch_degree(geo):
  """
  """
  bnodes = [len(b.neighbors) for b in geo.branches]
  return bnodes


def num_daughters(geo, branch=False, clean=True):
  """
  Only take daughters from the 1th end to be consistent (could crawl
  through whole tree structure but it's not clear that result would
  differ from this approach. clean=remove 0,1 nodes
  """
  branch_nodes = []
  
  if branch: # Use branches -- not such an excellent thingy
    for b in geo.branches:
      temp = []
      temp = [i for i in b.neighborLocations if i[0] == 1]
      branch_nodes.append(len(temp))
  
  else: # Use segments -- this is better (but much slower)
    pDF = PathDistanceFinder(geo, geo.soma)
    for seg in geo.segments:
      # Add the number of daughters at the farther neighbor
      if pDF.distanceTo(seg, 1) > pDF.distanceTo(seg, 0): # 1 is farther
        temp = len(seg.neighborsAt(1))
      else:
        temp = len(seg.neighborsAt(0))
      branch_nodes.append(temp)
  
  if clean:
    branch_nodes = [i for i in branch_nodes if i != 0 and i != 1]
    
  return branch_nodes



def degree_percent(deglist):
  """
  Return the % of branch degrees in each furcation,
  eg: [2,2,2,2,3,3] --> 2: 66%, 3: 34%
  """
  uniq = list(set(deglist))
  counts = {i: float(deglist.count(i))/len(deglist) for i in uniq}
  return counts





#######################################################################
# Path length and tortuosity

def path_lengths(geo):
  tips, tipinds = geo.getTipIndices()
  pDF = PathDistanceFinder(geo, geo.soma, 0.5)
  tipsegs = [geo.segments[i] for i in tips]
  path = [pDF.distanceTo(x,y) for x, y in zip(tipsegs, tipinds)]
  tort = [pDF.tortuosityTo(x,y) for x, y in zip(tipsegs, tipinds)]
  return path, tort


def path_lengths2(geo):
  # if FilamentIndex != geo.segments[index], use this: 
  tips, tipinds = geo.getTipIndices()
  tipsegs = [i for i in geo.segments if geo.getFilamentIndex(i) in tips]
  pDF = PathDistanceFinder(geo, geo.soma, 0.5)
  path, tort = [], []
  for x, y in zip(tipsegs, tipinds):          
    try:
      p, t = pDF.distanceTo(x,y), pDF.tortuosityTo(x,y)
      path.append(p)
      if t < 1: # Tortuosity cannot be < 1
        t = 1
      tort.append(t)
    except:
      continue
  return path, tort



def wen_tortuosity(geolist):
  """
  Returns pairs = [Euclidean, PathLengths] list, should be ~ 1, and 
  normalized tortuosity index (T) ~ normalized path length (l), 
  should be T = l/R -1. Inspired by figure 4E&F from Wen...Chklovskii, 2009
  """
  import scipy.stats as stats
  def seg_int(name):
    try:
      beh = int(name.split('[')[1].split(']')[0])
      return beh
    except:
      try:
        beh = int(name.split('_')[-1])
        return beh
      except:
        print(name)
    return 
  edists, plengths = [], []
  for g in geolist:
    tipsegs, tiplocs = g.getTipIndices()
    pDF = PathDistanceFinder(g, g.soma)
    tempedists, tempplengths = [], []
    for seg in g.segments:
      if seg_int(seg.name) in tipsegs:
        tempedists.append(dist3(g.soma.coordAt(0.5), seg.coordAt(0.5)))
        tempplengths.append(pDF.distanceTo(seg, 0.5))
    edists.append(np.mean(tempedists))
    plengths.append(np.mean(tempplengths))
  return edists, plengths
  
    
  
  




#######################################################################
# Sholl stuff

# Currently, only Hooser Sholl is used, which gives the most accurate results

def interpoint_dist(geo):
  # NOT CURRENTLY USED
  # determine the distances between successive points
  def nodex(node):
    return [node.x, node.y, node.z]
  dists = []
  for s in geo.segments:
    for n in range(len(s.nodes)-1):
      dists.append(dist3(nodex(s.nodes[n]), nodex(s.nodes[n+1])))
  print('Mean distance (%i points): %.5f +/- %.5f' 
         %(len(dists), np.mean(dists), np.std(dists)))
  return dists


def can_branch_be_added_again(geo, branch, sholl, key):
  # NOT CURRENTLY USED
  # profile branch
  soma = geo.soma.nodeAt(0)
  dD = [np.sign(node3(branch.nodes[n+1],soma)-node3(branch.nodes[n],soma))
        for n in range(len(branch.nodes)-1)]
  dLoc = [node3(i,soma) for i in branch.nodes[:-1]]
  changepts = []
  last_i = dD[0]
  for i in dD:
    if np.sign(i) != np.sign(last_i):
      changepts.append([dLoc[dD.index(i)], np.sign(i)])
      last_i = np.sign(i) # if there's a change, update the i
  # now figure out how many can be added
  if len(changepts) == 0:
    return False
  # if it only loops back once and hasn't been added twice, add it again
  elif len(changepts) == 1:
    if changepts[0][1] == -1: # if loops back
      if float(key) > changepts[0][0]:
        if sholl[key][1].count(branch) <= 1: # make 
          return True
    elif changepts[0][1] == 1: # if loops OUT
      if float(key) > changepts[0][0]:
        if sholl[key][1].count(branch) <= 1:
          return True
  else: # multiple change pts
    if sholl[key][1].count(branch) >= len(changepts) + 1:
      return False # already used all its changepts
    else:
      # assume forward order
      c = sholl[key][1].count(branch)
      change_sholls = [i[0] for i in changepts]
      if float(key) > min(change_sholls) and float(key) < max(change_sholls):
        return
  return


def interpolate_nodes(geo, return_nodes=False):
  # NOT CURRENTLY USED
  # find the most common distance betwixt successive nodes and then,
  # when successive nodes leave integer multiples of this distance
  # interpolate the difference to 'even' it out
  def nodex(node):
    return [node.x, node.y, node.z]
    
  def interp(pt1, pt2, ints):
    Xs = np.linspace(pt1[0], pt2[0], ints)
    Ys = np.linspace(pt1[1], pt2[1], ints)
    Zs = np.linspace(pt1[2], pt2[2], ints)
    return [[Xs[i],Ys[i],Zs[i]] for i in range(len(Xs))]
    
  dist = np.median(interpoint_dist(geo))
  pts = []
  segcount = -1
  for s in geo.segments:
    segcount = segcount + 1
    if segcount % 100 == 0:
      print('Completed %i/%i segments ' 
             %(segcount,len(geo.segments)))
    for n in range(len(s.nodes)-1):
      # if too far between nodes, interpolate
      if dist3(nodex(s.nodes[n]),nodex(s.nodes[n+1])) > 2 * dist:
        integer_interpolate = int((dist3(nodex(s.nodes[n]),
                                         nodex(s.nodes[n+1])))
                                   /dist)
        new_pts = interp(nodex(s.nodes[n]),nodex(s.nodes[n+1]),
                         integer_interpolate)
      # else just add the regular node pts
      else:
        new_pts = [nodex(s.nodes[n]), nodex(s.nodes[n+1])]
      # add the points as long as they don't already exist in pts
      for p in new_pts:
        if p not in pts:
          pts.append(p)
  # now should have all the points
  soma = geo.soma.coordAt(0.5)
  distances = []
  for p in pts:
    distances.append(dist3(soma, p))
  if return_nodes:
    return distances, pts
  else:
    return distances


def hooser_sholl(geo, sholl_lines=1000):
  """
  Sholl analysis without repeats or missed counts. Sholl_lines can be
  integer, whereupon the program creates that many evenly-spaced radii,
  or can be a vector whose values are used for sholl radii.
  """
  def get_neighbors(geo, branch, neb_segs):
    # For each possible neighbor, if it shares more than 1 node in common
    # with the branch, don't add it (it's probably the branch itself)
    possible_nodes = branch.nodes
    possible_segs, possible_locs = [], []
    for s in neb_segs:
      count = 0
      for n in s.nodes:
        if n in possible_nodes:
          count = count + 1
          if s.nodes.index(n) == len(s.nodes):
            possible_locs.append(-1)
          else:
            possible_locs.append(s.nodes.index(n))
          possible_segs.append(s)
      if count > 1:
        for j in range(count):
          possible_segs.pop()
          possible_locs.pop()
    # get the actual location of
    next_locs = []
    for i in possible_locs:
      if i == 0:
        next_locs.append(1)
      elif i == -1:
        next_locs.append(-2)
    possible_segs, possible_locs = possible_segs[0], possible_locs[0]
    return possible_segs.nodes[possible_locs]
    
  # helper functions
  def cross_node(geo, branch, nodeNum, sholl):
    # Determine if a sholl line is crossed by this node->next node
    soma = geo.soma.nodeAt(0)
    nebs = None
    #try: # if there is a next node, add it
    next_dist = node3(branch.nodes[nodeNum+1], soma)
    #except: # if not, find a ('the'?) neighboring node
    #  print('error!')
    #  try:
    #    nebs = [i[2] for i in branch.neighborLocations if i[0]==1]
    #  except:
    #    nebs = None
    # if this is the last node and there and neighbors with this node,
    # get the next node from the neighbors
    #if nebs:
    #  next_dist = get_neighbors(geo, branch, nebs)
    #else:
    #  return sholl
    node_dist = node3(branch.nodes[nodeNum], soma)
    for i in [float(k) for k in sholl.keys()]:
      if i < max([next_dist, node_dist]) and \
        i > min([next_dist, node_dist]):
          sholl[str(i)][0] = sholl[str(i)][0] + 1
          sholl[str(i)][1].append(branch)
    return sholl
      
  soma = geo.soma.nodeAt(0)
  sholl = {} # dictionary of crossings as keys, numcrossings as [0] (int)
             # and branches that cross (as objects? names?) as [1] (list)
  # integer mode
  if type(sholl_lines) is int:
    dists = [node3(n,soma) for n in geo.nodes]
    dists.sort()
    d99 = dists[int(len(dists)*.99)] # get 99% of the nodes
    lines = np.linspace(min(dists), d99, sholl_lines)
    for l in lines:
      sholl[str(l)] = [0,[]]
  # list mode
  elif type(sholl_lines) is list or type(sholl_lines) is np.ndarray:
    for l in sholl_lines:
      sholl[str(l)] = [0, []]
  else:
    print('sholl_lines input must be int or list or ndarray!')
    print('instead got %s' %str(type(sholl_list)))
    return None
  # go through branches and nodes and tabulate crossings
  for b in geo.branches:
    for nodeNum in range(len(b.nodes)-1): # here was change -1
      sholl = cross_node(geo, b, nodeNum, sholl)
  sholl_keys = list(sholl.keys())
  float_keys = [float(i) for i in sholl_keys]
  float_keys.sort()
  sholl_count = [sholl[str(i)][0] for i in float_keys]
  
  return [float_keys, sholl_count], sholl



######################################################################    
# Sub-tree analysis

# Note that sub-tree analysis works in tandem with Matlab for easier
#   3-D plotting and is somewhat interactive (it requires actions from
#   the user). The necessary Matlab scripts are ShollColor and ShollAxons.

# ALL SUBTREE ANALYSIS HAS BEEN MOVED TO neuron_Subtrees.py SINCE
# IT WAS GETTING VERY LONG.




### Subtree angles
"""def subtree_angles(subtrees, path, keep=10):
  # Returns the first N('keep') successive branch angles for each subtree.
  angles = []
  for t in subtrees:
    root = union(t, path) # Find the first non-mainpath segment
    count = 0
    prev_segs = [p for p in path]
    prev_segs.append(root)
    currsegs = [root]
    tree_angles = []
    while count < keep: # Keep going until 'keep' or run out of segs
      curr_angles = [] # Keep a list of angles for the current round
      next_segs = []
      for curr in currsegs:
        for n in curr.neighbors:
          if n not in prev_segs:
            pts = find_points(n, curr)
            curr_angles.append(get_angle(pts[0], pts[1], pts[2]))
            next_segs.append(n)
        prev_segs.append(curr)
      curr_segs = next_segs
      try:
        tree_angles.append(np.mean(curr_angles))
      except:
        count = 10
      count = count + 1
    # curr_angles should have <= 10 lists of angles, now these are averaged
    # angles.append([np.mean(o) for o in 
  return




## Still work in progress
def dendrogram_subtrees(geo, paths, subtrees):
  def plot_next(segs, prev_segs, x_range, prev_pts):
    if len(segs) != len(prev_pts):
      print('segs (%i) and prev_pts (%i) should be same length'
            %(len(segs), len(prev_pts)))
      return
    
  def subtree_pts(sub, x_range, paths):
    for s in sub:
      if s in np.array(paths).flatten():
        root = s
    
    return
    
  if len(paths) > 10:
    # There is only one path (one axon)
    m = len(subtrees)
  
  return
"""
  



######################################################################
# Partition asymmetry (symmetry)


# This version is CURRENTLY used.
def path_asymmetry(geo, bOrder=True):
  """
  This version uses paths and sets to calculate the length asymmetry.
  """
  geo.calcBranchOrder()
  pDF = PathDistanceFinder(geo, geo.soma, 0) # initialize pDF
  tipSegs, tipLocs = geo.getTipIndices()
  paths = [] # get all the paths
  for s in range(len(tipSegs)):
    try:
      paths.append(pDF.pathTo(geo.segments[tipSegs[s]],tipLocs[s]))
    except:
      continue
  length_asymmetry, orders = [], []
  for p1 in paths:
    if paths.index(p1) % 100 == 0:
      print('completed %i (of %i) paths' %(paths.index(p1),len(paths)))
    for p2 in paths:
      if p1 != p2: # make sure not the same
        p_root = [i for i in p1 if i in p2]
        # make sure they share common segments, then do analysis
        if len(p_root) > 0:
          l_p1 = sum([s.length for s in p1])
          l_p2 = sum([s.length for s in p2])
          l_root = sum([s.length for s in p_root])
          #
          if l_p1 > l_p2: # always put greater in denominator (customary)
            if l_p2 == 0:
              print('missed one')
              pass
            else:
              length_asymmetry.append((l_p2-l_root)/(l_p1-l_root))
              if bOrder is True:
                orders.append(path_branch_order(p1, p2))
          else:
            if l_p1 == 0:
              print('missed one')
            else:
              length_asymmetry.append((l_p1-l_root)/(l_p2-l_root))
              if bOrder is True:
                orders.append(path_branch_order(p1, p2))
  #if len(length_asymmetry) > 1000:
  #  length_asymmetry = length_asymmetry[::int(len(length_asymmetry)/1000)]
  if bOrder is True:
    return length_asymmetry, orders
  else:
    return length_asymmetry



def get_segment(geo, segname):
  for s in geo.segments:
    if s.name == segname:
      return s


def add_all_downstream(geo, seg, prev_seg, tips):
  # print('called add_all_downstream')
  newsegs = [seg]
  prevsegs = [prev_seg]
  for n in prev_seg.neighbors:
    if n != seg:
      prevsegs.append(n)
  go = True
  same = 0
  while go:
    old_len = len(newsegs)
    for n in newsegs:
      for neb in n.neighbors:
        if neb not in prevsegs:
          newsegs.append(neb)
      prevsegs.append(n)
    if len(newsegs) == old_len: # if no change in newsegs, increment
      same = same + 1
    else:
      same = 0
    if same > 10: # after 10 no-changes, stop
      go = False
  # get tips and lengths
  pDF = PathDistanceFinder(geo, seg, 0)
  numtips, length = 0, []
  for t in tips:
    try:
      if geo.segments[t] in newsegs:
        numtips = numtips + 1
        add_segs = pDF.pathTo(t,1)
        for a in add_segs:
          if a not in length:
            length.append(a)
    except:
      continue
  length = sum([l.length for l in length])
  return length, numtips



def path_branch_order(path1, path2):
  # Find the order of the last common branch
  # root = [i for i in path1 if i in path2]
  if path1[0] == path2[0]: # Should be soma
    # Get the first segment not in both (should be where they diverge)
    for p in path1:
      if p not in path2:
        return p.branchOrder



# OLD
def tips_asymmetry(geo):
  """
  Get the tip asymmetry of the neuron. Follow the soma's neighbors
  until there are more than 1, then start there.
  seg_lengths: dict with a section_name for keys and float as values
  seg_tips: dict with sec_name as key and list of segment objects as values
  """
  def get_bif_info(geo, seg, prev_seg, tips):
    # Given a branch and the previous branch (so we know the direction
    # of movement), make a dictionary with the non-previous neighbors as
    # keys and add all subsequent lengths and tips to the respective keys
    # print('called get_bif_info')
    forward_nebs = [n for n in seg.neighbors if n != prev_seg]
    neb_dict = {}
    lengths, tips = [], []
    for f in forward_nebs:
      length, tip = add_all_downstream(geo, f, seg, tips)
      # now it should be part, for future reference, too
      lengths.append(length)
      tips.append(tips)
      # neb_dict[f.name] = {'prevbranch': branch.name, 'length': 0, 'tips': 0}
    length_asym = []
    for l in lengths:
      try:
        length_asym.append(l/(sum(lengths)-l))
      except:
        continue
    #length_asym = [i/(sum(lengths)-i) for i in lengths]
    tip_asym = [] # [i/(sum(tips)-i) for i in tips]
    for t in tips:
      try:
        tip_asym.append(t/(sum(tips)-t))
      except:
        continue
    return length_asym, tip_asym
  
  # go through all branches in order
  tips, _ = geo.getTipIndices()
  #print('got tips')
  master_lengths, master_tips = [], []
  prevsegs = [geo.soma]
  newsegs = [i for i in geo.soma.neighbors if i not in prevsegs]
  print(newsegs)
  go, same = True, 0
  while go:
    old_len = len(newsegs)
    for n in newsegs:
      nebs = n.neighbors
      for k in nebs: # remove repeated segs
        if k in prevsegs:
          nebs.pop(nebs.index(k))
      if len(nebs) > 1: # now if there are > 1 neighbors (a branch)
        print('found a bifurcation!')
        for k in nebs:
          lengths, tips = get_bif_info(geo, n, k, tips)
          for i in lengths:
            master_lengths.append(i)
          for j in tips:
            master_tips.append(j)
          if k not in newsegs:
            newsegs.append(k) # done with k-th neighbor
      elif len(nebs) == 1:
        # if it's not a bifurcation, add it to prevsegs
        prevsegs.append(nebs[0])
      if n not in prevsegs:
        prevsegs.append(n) # done with n-th newseg
    if len(newsegs) == old_len:
      same = same + 1
      #print(same)
    else:
      same = 0
    if same > 1000:
      go = False
      return master_lengths, master_tips, prevsegs
    old_len = len(newsegs)
  # it should never get to this part, but in case it does
  return master_lengths, master_tips, prevsegs


def tip_coords(geo, seg_tips):
  # return x-y-z tuples for each tip; just use the (1) position of each tip seg
  tip_coords = {}
  for k in seg_tips.keys():
    tip_coords[k] = []
    for t in seg_tips[k]:
      tip_coords[k].append(t.coordAt(1))
  return tip_coords


def simplify_asymmetry(geo):
  # simplification of asymmetry data
  seg_lengths, seg_tips = tips_asymmetry(geo)
  sumlengths = sum([seg_lengths[k] for k in seg_lengths.keys()])
  sumtips = sum([len(seg_tips[k]) for k in seg_tips.keys()])
  lengths = [seg_lengths[k]/(sumlengths-seg_lengths[k]) for k in seg_lengths.keys()]
  tips = [float(len(seg_tips[k]))/float((sumtips-len(seg_tips[k]))) for k in seg_tips.keys()]
  return lengths, tips



def somatofugal_length(geo):
  print('Building lengths and tips...')
  prev_segs = [geo.soma]
  next_segs = [n for n in geo.soma.neighbors if n not in prev_segs]
  seg_lengths, seg_tips = {}, {}
  for n in next_segs:
    nlen, ntip = add_all_downstream(geo, n, geo.soma)
    seg_lengths[n.name] = nlen
    seg_tips[n.name] = ntip
  go, cnt = True, 0
  while go:
    old_len, new_next_segs = len(prev_segs), []
    for n in next_segs:
      for neb in n.neighbors:
        if neb not in prev_segs and neb not in next_segs:
          # now that have all the new segs, get their len & tips
          nlen, ntip = add_all_downstream(geo, neb, n)
          if neb.name in seg_lengths.keys() or neb.name in seg_tips.keys():
            print('Tried to add %s but is already there' %neb.name)
          else:
            seg_lengths[neb.name], seg_tips[neb.name] = nlen, ntip
          # add nebs to next_segs and n to prev_segs
          new_next_segs.append(neb)
      # put n into prev_segs
      prev_segs.append(next_segs.pop(next_segs.index(n)))
    # next_segs should be empty now
    next_segs = [i for i in new_next_segs]
    if len(prev_segs) == old_len:
      cnt = cnt + 1
    else:
      cnt = 0
    if cnt == 50 or len(prev_segs)==len(geo.segments):
      go = False
  print('Got %i (of %i) lengths and %i (of %i) tips)'
        %(len(seg_lengths), len(geo.segments),
          len(seg_tips), len(geo.segments)))
  return seg_lengths, seg_tips  



def tips_asymmetry_old(geo):
  """
  ############## OLD VERSION!!! #####################
  Get the tip asymmetry of the neuron. Follow the soma's neighbors
  until there are more than 1, then start there.
  seg_lengths: dict with a section_name for keys and float as values
  seg_tips: dict with sec_name as key and list of segment objects as values
  """
  seg_lengths, seg_tips = somatofugal_length(geo) # this will get a lot of them
  def get_bif_info(geo, seg, prev_seg, seg_lengths, seg_tips):
    # Given a branch and the previous branch (so we know the direction
    # of movement), make a dictionary with the non-previous neighbors as
    # keys and add all subsequent lengths and tips to the respective keys
    # print('called get_bif_info')
    forward_nebs = [n for n in seg.neighbors if n != prev_seg]
    neb_dict = {}
    lengths, tips = [], []
    for f in forward_nebs:
      # if this seg isn't already part of seg_lengths/tips, add it
      if f.name not in seg_lengths.keys():
        length, tip = add_all_downstream(geo, f, seg)
        seg_lengths[f.name] = length
      if f.name not in seg_tips.keys():
        seg_tips[f.name] = tip
      # now it should be part, for future reference, too
      lengths.append(seg_lengths[f.name])
      tips.append(seg_tips[f.name])
      # neb_dict[f.name] = {'prevbranch': branch.name, 'length': 0, 'tips': 0}
    length_asym = []
    for l in lengths:
      try:
        length_asym.append(l/(sum(lengths)-l))
      except:
        continue
    #length_asym = [i/(sum(lengths)-i) for i in lengths]
    tip_asym = [] # [i/(sum(tips)-i) for i in tips]
    for t in tips:
      try:
        tip_asym.append(t/(sum(tips)-t))
      except:
        continue
    return length_asym, tip_asym, seg_lengths, seg_tips
  
  # go through all branches in order
  master_lengths, master_tips = [], []
  prevsegs = [geo.soma]
  newsegs = [i for i in geo.soma.neighbors if i not in prevsegs]
  # print(newsegs)
  go, same = True, 0
  while go:
    old_len = len(newsegs)
    for n in newsegs:
      nebs = n.neighbors
      for k in nebs: # remove repeated segs
        if k in prevsegs:
          nebs.pop(nebs.index(k))
      if len(nebs) > 1: # now if there are > 1 neighbors (a branch)
        print('found a bifurcation!')
        for k in nebs:
          lengths, tips, seg_lengths, seg_tips = get_bif_info(geo, n, k, seg_lengths, seg_tips)
          for i in lengths:
            master_lengths.append(i)
          for j in tips:
            master_tips.append(j)
          if k not in newsegs:
            newsegs.append(k) # done with k-th neighbor
      if n not in prevsegs:
        prevsegs.append(n) # done with n-th newseg
    if len(newsegs) == old_len:
      same = same + 1
    else:
      same = 0
    if same > 10:
      go = False
      return master_lengths, master_tips, prevsegs
    old_len = len(newsegs)
  # it should never get to this part, but in case it does
  return master_lengths, master_tips, prevsegs
  


######################################################################
# Torques

def getNormVector(points):
  #print(points, np.shape(points))
  v1 = [points[1][0][i] - points[0][0][i] for i in range(3)]
  v2 = [points[2][0][i] - points[0][0][i] for i in range(3)]
  normVec = np.cross(v1,v2)
  return normVec


def angleBetween(plane1,plane2,planCoords):
  # get normal vectors
  n1, n2 = getNormVector(planCoords[plane1]), \
           getNormVector(planCoords[plane2])
  angle = np.arccos( (abs(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2])) /
                     ( np.sqrt(n1[0]**2+n1[1]**2+n1[2]**2) *
                       np.sqrt(n2[0]**2+n2[1]**2+n2[2]**2) ) )
  return angle*180/np.pi


def get_torques(geo):
  # return bifurcation torques
  Cons =  geo.connections
  Seg1s, Seg2s = [], []
  for c in Cons:
    Seg1s.append(c['filament1']) # here, location1 is always 0
    Seg2s.append(c['filament2']) # here, location2 is always 1
    #geometry.c['filament1'].coordAt(c['location1'])
  #
  tsegs = np.array([Seg1s,Seg2s]).T
  tsegs = tsegs.reshape(len(tsegs)*2)
  segs = set(tsegs)
  planCoords = {}
  #
  count = 0
  for seg in segs:
    friends, friendcoords = [], []
    for s in geo.segments:
      if s.name == seg:
        friends.append(s.name)
        if s.name in Seg1s:
          friends.append(Seg2s[Seg1s.index(s.name)])
        if s.name in Seg2s:
          friends.append(Seg1s[Seg2s.index(s.name)])
    #print('friends compiled')
   #   
    for s in geo.segments:
      if s.name in friends:
        friendcoords.append([s.coordAt(1)])
    count = count + 1
    #if count%100 == 0:
    #  print('%i of %i segments done' %(count, len(segs)))
    if len(friendcoords) > 2: # need 3 points to define plane
      planCoords[seg]=friendcoords
  #
  planCoordskeys = []
  for s in geo.segments: # loop through segments to find plane-neighbors
    if s.name in planCoords.keys():
      for n in s.neighbors:
        if n.name in planCoords.keys(): # if the neighbor is also a bifurcation
          planCoordskeys.append([s.name, n.name]) # add it
        else: # otherwise, keep looking for a neighbor that is
          for nn in n.neighbors:
            if nn.name in planCoords.keys():
              planCoordskeys.append([s.name, nn.name])
  #
  # get torques
  torques = []
  for P in planCoordskeys:
    torques.append(angleBetween(P[0],P[1],planCoords))
  return torques




###############################################################################
# Ellipse fitting, distance to nearest point stuff

# Neuropil fitting has replaced ellipse fitting, included here for fun
def neuropil_fit(geo):
  """
  Coarsely find which cubes the neuron contacts, then make finer cubes
  to see how well it fills those finer cubes. The first method generally
  finds where the neuron is, the second part tests how well it fills
  that space.
  """
  nodes = geo.nodes[::100]
  minx, maxx, miny, maxy, minz, maxz = np.inf, 0, np.inf, 0, np.inf, 0
  print('Retrieving the bounds for %s' %geo.name)
  def check_pt(pt, refmin, refmax):
    if pt < refmin:
      refmin = pt
    if pt > refmax:
      refmax = pt
    return refmin, refmax
  for n in nodes:
    minx, maxx = check_pt(n.x, minx, maxx)
    miny, maxy = check_pt(n.y, miny, maxy)
    minz, maxz = check_pt(n.z, minz, maxz)
  minmax = [[minx, maxx],[miny, maxy],[minz, maxz]]
  # Create cubic grid, 10 um^3
  biggrid = [np.linspace(h[0],h[1],int((h[1]-h[0])/10.)) for h in minmax]
  # Determine whether a cube contains any nodes
  print('Creating the large grid')
  n_cubes = len(biggrid[0])*len(biggrid[1])*len(biggrid[2])
  cube_idx = {}
  for i in biggrid[0]:
    for j in biggrid[1]:
      for k in biggrid[2]:
        cube_idx[i,j,k] = 0
  def in_cube(cube_idx, n, thresh=5.):
    for k in cube_idx.keys():
      if dist3(k, [n.x, n.y, n.z]) <= thresh: # 5 um
        cube_idx[k] = 1
    return cube_idx
  # Populate the cubes
  # nodes = geo.nodes[::10]
  for n in nodes:
    cube_idx = in_cube(cube_idx, n)
  # For each filled cube, divide it into smaller cubes!
  def three_randos(center):
    pts = []
    for p in range(10):
      pts.append([np.random.random(1)*10. + center[0]-5., # x
                  np.random.random(1)*10. + center[1]-5., # x
                  np.random.random(1)*10. + center[2]-5.]) # z
    return pts # Return the 10 random pts
  print('Picking the random points')
  micron_pts = []
  for k in cube_idx.keys():
    if cube_idx[k] == 1:
      temp_pts = three_randos(k)
      for t in temp_pts:
        micron_pts.append(t)
  # Now have all the points, get a random smattering of 1000 pts and find their distances
  r_pts = list(np.unique([int(np.random.random(1)*len(micron_pts)) for i in range(1000)]))
  r_pts = [micron_pts[u] for u in r_pts]
  pt_dists, count = [], 0
  for p in r_pts:
    pt_dists.append(min([dist3(p,[n.x,n.y,n.z]) for n in geo.nodes]))
    if count%100 == 0:
      print('%i/%i points processed so far' %(count, len(r_pts)))
    count = count + 1
  return pt_dists

#
#
# All ellipse stuff below this is OLD!
def quad_fit(geo):
  """
  Get the bounds of the neuropil (from skeleton) and fit a simple quadratic.
  Fits points at 20 um intervals.
  """
  minx, maxx, miny, maxy, minz, maxz = np.inf, 0, np.inf, 0, np.inf, 0
  def check_pt(pt, refmin, refmax):
    if pt < refmin:
      refmin = pt
    if pt > refmax:
      refmax = pt
    return refmin, refmax
  for n in geo.nodes:
    minx, maxx = check_pt(n.x, minx, maxx)
    miny, maxy = check_pt(n.y, miny, maxy)
    minz, maxz = check_pt(n.z, minz, maxz)
  mid = [np.mean([n.x for n in geo.nodes]),
         np.mean([n.y for n in geo.nodes]),
         np.mean([n.z for n in geo.nodes])]
  # Use mean-centered for quadratic fit (?)
  minx, maxx = minx - mid[0], maxx - mid[0]
  miny, maxy = miny - mid[1], maxy - mid[1]
  minz, maxz = minz - mid[2], maxz - mid[2]
  #
  def get_XYgrid(xmin, xmax, ymin, ymax, zval):
    gridpts = [] # Populate this with the grid 
    # Fit XY (minx, 0), (0, maxy), (maxx, 0); & reflect for bottom part
    coefXY = np.polyfit([xmin, 0, xmax], [0, ymax, 0], 2) # Force quad
    numxsteps = int((xmax - xmin)/20.) # max - (-min) = max+abs(min)
    xsteps = linspace(xmin, xmax, numxsteps)
    ylims = [sum([coefXY[0]*x**2, coefXY[1]*x, coefXY[2]]) for x in xsteps]
    for i in range(len(ylims)):
      yvals = np.linspace(0, ylims[i], int(ylims[0]/20.)) # positive y
      for j in yvals:
        gridpts.append([xstep[i], j, zval])
        gridpts.append([xstep[i], -j, zval])
        gridpts.append([-xstep[i], j, zval])
        gridpts.append([-xstep[i], -j, zval])
    return gridpts # All grid points for this level of z
  # Fit ZY and ZX to get the 
  return
  
  


# old ellipse stuff; nothing below this line is currently used for ellipse

def getNoSomaPoints(geo):
  # get the downsampled nodes, but not the soma
  somaPos = geo.soma.coordAt\
            (geo.soma.centroidPosition(mandateTag='Soma'))
  print('Soma position: %.5f, %.5f, %.5f' %(somaPos[0],somaPos[1],somaPos[2])) # works
  nodes = []
  for seg in geo.segments:
    if 'Soma' not in seg.tags:
      nodes.append(seg.coordAt(0))
      nodes.append(seg.coordAt(0.5))
      nodes.append(seg.coordAt(1))
  print('Sampled %i nodes' %len(nodes))
  
  return nodes



def findBounds(nodelist):
  # return the x,y,z bounds of the node list
  xs, ys, zs = [], [], []
  
  for n in range(len(nodelist)):
    xs.append(nodelist[n][0])
    ys.append(nodelist[n][1])
    zs.append(nodelist[n][2])

  bounds = {'xmin': min(xs), 'xmax': max(xs), 
            'ymin': min(ys), 'ymax': max(ys),
            'zmin': min(zs), 'zmax': max(zs)}
  
  return bounds



def getGridPoints(nodelist, pplot=False):
  # create a grid around the neuropil and use linspace to fill the volume
  bounds = findBounds(nodelist)
  gridpoints = []
  xs = np.linspace(bounds['xmin'], bounds['xmax'], 10)
  ys = np.linspace(bounds['ymin'], bounds['ymax'], 10)
  zs = np.linspace(bounds['zmin'], bounds['zmax'], 10)
  spacing = xs[1]-xs[0]
  
  # 1000 grid volume points
  for i in range(len(xs)-1):
    for j in range(len(ys)-1):
      for k in range(len(zs)-1):
        gridpoints.append([(xs[i]+xs[i+1])/2,
                           (ys[j]+ys[j+1])/2,
                           (zs[k]+zs[k+1])/2])
  print('gridpoints is length %i' %len(gridpoints))
  
  boxx, boxy, boxz = [], [], []
  for g in range(len(gridpoints)):
    boxx.append(gridpoints[g][0])
    boxy.append(gridpoints[g][1])
    boxz.append(gridpoints[g][2])
  
  nodex, nodey, nodez = [], [], []
  for n in range(len(nodelist)):
    nodex.append(nodelist[n][0])
    nodey.append(nodelist[n][1])
    nodez.append(nodelist[n][2])
  
  if pplot:
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
  #ax.plot(boxx, boxy)
    ax1.scatter(boxx, boxy, boxz, color='r', marker='.', alpha=0.5)
    ax1.scatter(nodex, nodey, nodez, color='k', marker='.', alpha=1)
  # ax.set_xlabel('')
  # plt.show()
    
  return gridpoints, spacing
  


def closestPoint(rectpoint, nodes):
  # find the closest neuron node to a rectangle point
  ptmin = np.inf
  ptind, pt = None, None
  for n in range(len(nodes)):
    dist = dist3(rectpoint, nodes[n])
    if dist < ptmin:
      ptmin = dist
      ptind = n
      pt = nodes[n]
  return ptind, ptmin


def closestPointPool(things):
  # find the closest neuron node to a rectangle point
  # things[0] = rect point, things[1] = all nodes
  things[0] = rectpoint
  things
  ptmin = np.inf
  ptind, pt = None, None
  for n in range(len(nodes)):
    dist = dist3(rectpoint, nodes[n])
    if dist < ptmin:
      ptmin = dist
      ptind = n
      pt = nodes[n]
  return ptind, ptmin


def getSurfacePoints(gridpoints, nodes, spacing, pplot=False):
  # given volume points and neuropil nodes, create downsampled
  # volume of the neuropil (if a neuron point is in a given cube, 
  # the cube is a 1, else 0
  ellipsePoints = []
  if type(gridpoints) is not np.ndarray:
    gridpoints = np.array(gridpoints)
  if type(nodes) is not np.ndarray:
    nodes = np.array(nodes)
  
  for b in range(len(gridpoints)):
    _, dist = closestPoint(gridpoints[b], nodes)
    if dist <= spacing/8.:
      ellipsePoints.append(gridpoints[b])
    if b % 100 == 0 and b != 0:
      print('%i/%i points examined' %(b, len(gridpoints)))
      
  print('Now have %i neuropil points' %len(ellipsePoints))
  
  surfx, surfy, surfz = [], [], []
  for s in ellipsePoints:
    surfx.append(s[0])
    surfy.append(s[1])
    surfz.append(s[2])
  if pplot:
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.scatter(surfx, surfy, surfz, color='g', marker='.')
    ax2.set_xlabel('X axis')
    ax2.set_ylabel('Y axis')
    ax2.set_zlabel('Z axis')
    plt.show()
  
  return ellipsePoints


def writeFile(points, outfile):
  # write points to a ascii; this is generally not necessary
  if outfile is None:
    outfile = 'neuropil_surfpoints.txt'  
  with open(outfile, 'w') as fOut:
    for p in range(len(points)):
      # print(points[p])
      ptstring = [str(points[p][0]), str(points[p][1]), str(points[p][2])]
      ptstr = ' '.join(ptstring)
      fOut.write(ptstr)
      fOut.write('\n')
      #print
  fOut.close()
  print('%s file written.' %outfile)
  return


# Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz = 1

def give_ellipse(axes, shrink, translate):
  """
  axes: [1x3], shrink: scalar (ratio), translate: [1x3]
  Returns a 2-D ellipse of points when given the 3 axes ([maj, min, 3rd])
  and where on the 3rd axis the current slice is
  --> axes = original evals ### scale omitted here
  --> 'shrink' is the ratio that determines 
      how large how the ellipse should be stretched in 2-D
  --> axes[2] not used in this version
  """
  norm_ax = [i/max(axes) for i in axes]
  xs = np.linspace(-norm_ax[0],norm_ax[0],1000)
  ys = [np.sqrt( (1 - (i**2/norm_ax[0])) * norm_ax[1] ) for i in xs]
  # get rid of the nans
  opts = [[x,y] for x,y in zip(xs,ys) if np.isfinite(y)]
  # need to get the negative part of the y half of the graph
  pts = []
  for p in opts:
    pts.append([p[0],-p[1]])
    pts.append(p)
  # pts are currently the 'largest' possible, need to shrink by 'where'
  pts = np.array(pts)
  pts = pts * shrink
  newpts = []
  for p in pts:
    _pt = [axes[0] * p[0] + translate[0],  \
           axes[1] * p[1] + translate[1],  \
           translate[2]]
    if _pt not in newpts:
      newpts.append(_pt)
  
  return newpts


def get_reduced_points(geo, outfile=None):
  # only pre-req is to run getNoSomaPoints first
  nodes = getNoSomaPoints(geo)
  gridpoints, spacing = getGridPoints(nodes)
  ellipsePoints = getSurfacePoints(gridpoints, nodes, spacing)
  #writeFile(ellipsePoints, outfile)
  
  return ellipsePoints


def check_eigen(s_vals, s_vecs, pts):
  """
  For singular value decomposition, check the orientations of vectors
  vs. the points they're supposed to represent
  """
  # Get zero-centered points first
  #means = [pts[i] for i in range(len(pts)) if i%100==0] # downsample
  means = pts
  _m = [np.mean([j[0] for j in means]), np.mean([j[1] for j in means]),
        np.mean([j[2] for j in means])]
  # subtract the mean but keep the shape
  newmeans = []
  for m in means:
    newmeans.append([m[0]-_m[0],m[1]-_m[1],m[2]-_m[2]])
  dmax = farthest_pt(pts)
  # get eigenvectors normalized by distance from farthest pts
  scales = [i/max(s_vals)*dmax for i in s_vals]
  print(scales)
  v1 = [[0,0,0],[scales[0]*s_vecs[0][0], scales[0]*s_vecs[1][0], 
                 scales[0]*s_vecs[2][0]]]
  v2 = [[0,0,0],[scales[1]*s_vecs[0][1], scales[1]*s_vecs[1][1], 
                 scales[1]*s_vecs[2][1]]]
  v3 = [[0,0,0],[scales[2]*s_vecs[0][2], scales[2]*s_vecs[1][2], 
                 scales[2]*s_vecs[2][2]]]
  print(v1,v2,v3)
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  
  for m in newmeans:
    ax.scatter(m[0],m[1],m[2], c='b', edgecolor='b', alpha=0.2)
  ax.plot([0,v1[1][0]], [0,v1[1][1]], [0,v1[1][2]], c='r')
  ax.plot([0,v2[1][0]], [0,v2[1][1]], [0,v2[1][2]], c='g')
  ax.plot([0,v3[1][0]], [0,v3[1][1]], [0,v3[1][2]], c='k')
  plt.show()
  return newmeans


def build_ellipse(geo):
  """
  Uses singular values from a uniformly resampled neuron grid to get
  major/minor axes to create an ellipsoid; scales and translates the
  ellipsoid back to neuron space.
  """
  gpts = get_reduced_points(geo)
  gmean = [np.mean([i[0] for i in gpts]),
           np.mean([i[1] for i in gpts]),
           np.mean([i[2] for i in gpts])]
  # get singular values
  _, s_vals, s_vecs = np.linalg.svd(gpts)
  s = np.array([i/max(s_vals) for i in s_vals])
  # scale singular values by longest distance
  dmax = farthest_pt(gpts)
  s = s * dmax
  # hyperbolic scaling reference for taper of top/bottom
  _x = np.linspace(0,10,50)
  _y = -_x**2 + 100
  y = [i/max(_y) for i in _y]
  y.reverse()
  zscale = [i for i in y]
  y.reverse()
  for i in y:
    zscale.append(i)
  eig_pts = []
  # make 100 layers of v3
  zlayers = np.linspace(-s[2],s[2],100)
  for v in zlayers:
    newpts = give_ellipse(s, zscale[list(zlayers).index(v)], 
                          [0,0,0])
    for p in newpts:
      eig_pts.append(p)
  eig_pts = np.array(eig_pts)
  # now have all eigen points, need to re-orient axes
  pts = eig_pts.dot(np.linalg.inv(s_vecs))
  # now translate:
  pts = [[p[0]+gmean[0], p[1]+gmean[1], p[2]+gmean[2]] for p in pts]
  return pts, gpts, eig_pts
  
    
def get_distances(geo, multi=None):
  """
  Return the "distances", the distance from each ellipse point to the
  closest point of the neuron's skeleton. Can be parallelized by multi=int.
  """
  if multi is None:
    ellipse_pts, _, _ = build_ellipse(geo)
    nodes = getNoSomaPoints(geo)
    distances = []
    # make sure this isn't ridiculously long
    if len(ellipse_pts) > 100000:
      ellipse_pts = ellipse_pts[::10]
    for e in ellipse_pts:
      _, d = closestPoint(e, nodes)
      distances.append(d)
      if ellipse_pts.index(e)%1000==0:
        print('%i (of %i) sampled so far' %(ellipse_pts.index(e), len(ellipse_pts)))
    return distances
  elif type(multi) is int:
    from multiprocessing import Pool
    p = Pool(multi)
    #distances = pool.map(closestPointPool, 
  return distances
  


#######################################################################
# simple branch stuff

def branch_lengths(geo, locations=False):
  lengths = [b.length for b in geo.branches]
  locations = [b.coordAt(0.5) for b in geo.branches]
  if locations:
    return lengths, locations
  else:
    return lengths


def branch_order(geo):
  geo.calcForewardBranchOrder()
  return [b.branchOrder for b in geo.branches]
    

def length_vs_dist(geo=None, lengths=None, locations=None):
  # Simple calculation of branch lengths
  if lengths is None and locations is None and geo is not None:
    lengths, locations = branch_lengths(geo, True)
  midpt = [np.mean([i[0] for i in locations]),
           np.mean([i[1] for i in locations]),
           np.mean([i[2] for i in locations])]
  print('Calculating distances....')
  distances = [dist3(pt, midpt) for pt in locations]
  return lengths, distances
    
  
# Multiple-regress plot
def multi_regress(lol):
  # lol must be a list of lists (each item is a list of 2 lists)
  return




#######################################################################
# Tip-to-tip distances


# This uses paths
def tip_to_tip(geo):
  """
  Who knows -- this might be important some day.
  """  
  tips, tipInds = geo.getTipIndices()
  tip_dists = []
  Tips = list(zip(tips, tipInds))
  randtip = int(np.random.random(1)*len(tips))
  pDF = PathDistanceFinder(geo, geo.segments[tips[randtip]], tipInds[randtip])
  for t in range(len(Tips)):
    try:
      tip_dists.append(pDF.distanceTo(geo.segments[Tips[t][0]], Tips[t][1]))
    except:
      print('missed one')
  return tip_dists



# tip dist matrix
def tip_matrix(geo, show=True):
  """
  Returns a matrix of tip distances (i,i along diagonal).
  """
  filamentInds, tipLocs = geo.getTipIndices()
  tipSegs = []
  for f in filamentInds:
    for seg in geo.segments:
      if seg.filamentIndex == f:
        tipSegs.append(seg)
  dists = np.zeros([len(tipSegs), len(tipSegs)])
  for t in range(len(tipSegs)):
    pDF = PathDistanceFinder(geo, tipSegs[t])
    for i in range(len(tipSegs)):
      dists[t,i] = pDF.distanceTo(tipSegs[i])
  if show is True:
    print('Plotting ...')
    showntips = int(len(tipSegs)/10)
    mask = np.zeros_like(dists)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style('white'):
      ax = sns.heatmap(dists, mask=mask, xticklabels=showntips, 
                       yticklabels=showntips, cmap='RdYlGn_r')
                       # cmap: 'YlGnBu'
    plt.show()
  return dists



def where_tips(geo, returnCoords=False):
  """
  This returns how far each tip is from the center of neuropil.
  """
  mid = [np.mean([n.x for n in geo.nodes]),
         np.mean([n.y for n in geo.nodes]),
         np.mean([n.z for n in geo.nodes])]
  tipFils, tipLocs = geo.getTipIndices()
  tipCoords = []
  for s in geo.segments:
    if s.filamentIndex in tipFils:
      tipCoords.append(s.coordAt(tipLocs[tipFils.index(s.filamentIndex)]))
  dists = [dist3(mid, t) for t in tipCoords]
  if returnCoords is True:
    return dists, tipCoords
  else:
    return dists


# Euclidean distance
def tip_to_tip_euclid(geo):
  """
  This tells us how far in euclidean space a random tip is from other tips.
  """
  filamentInds, tipLocs = geo.getTipIndices()
  tipSegs = []
  for f in filamentInds:
    for seg in geo.segments:
      if seg.filamentIndex == f:
        tipSegs.append(seg)
  if len(tipSegs) != len(tipLocs):
    print('Should be same number of tip segs (%i) and tip locs (%i):'
          %(len(tipSegs),len(tipLocs)))
    return None
  # Get the reference segment at random
  refNum = int(np.random.random(1)*len(tipSegs))
  refPt = tipSegs[refNum].coordAt(tipLocs[refNum])
  # Get the tip distances
  tip_dists = [dist3([refPt[0], refPt[1], refPt[2]],
                     [i for i in tipSegs[u].coordAt(tipLocs[u])]) 
                     for u in range(len(tipSegs))]
  return tip_dists
  




#######################################################################
# Fractal dimension                  (as per Caserta et al., 1995)

# For every point in the (resampled/interpolated) neuron, basically do
#   a center of mass calculation for 'radius of gyration'
#   (doesn't need to be interpolated -- sampled from every neuron point)
# If the result doesn't make sense, be sure to look at the log-log plot
#   and make sure the slope is being taken from the correct location.
#   Sometimes, especially near the higher end of the log x-axis, slope 
#   changes are small and so you get a value near 1.0 (when it should be 1.3-1.6)

# Helper functions
def element_histogram(data, bins):
  # Also returns the actual data, list (of len(bins)-1) of lists
  # This might be bottle neck -- omitted from function for now
  thing = [[] for i in bins]
  noplace = 0
  for d in data:
    for b in range(len(bins)):
      if d >= bins[-1]: # this should never really happen except maybe once
        thing[-1].append(d)
      elif d <= bins[0]:
        thing[0].append(d)
      elif d >= bins[b] and d < bins[b+1]:
        thing[b].append(d)
      else:
        noplace = noplace + 1
  if len(thing[-1]) > 0:
    print('lost %i points from last bin' %len(thing[-1]))
  print('Could not place %i (of %i) points' %(noplace, sum([len(i) for 
                                              i in thing])))
  return thing[:-1]



def pt8_slope(bins, masses, where=False):
  # Get the average constant 8-pt slope
  if len(bins) > len(masses):
    bins = bins[:len(masses)]
  if len(masses) > len(bins):
    masses = masses[:len(bins)]
  bins, masses = [np.log(i) for i in bins], [np.log(i) for i in masses]
  slopes = [(masses[i+1]-masses[i])/(bins[i+1]-bins[i]) for i in
                                                        range(len(bins)-1)]
  delta = [abs(slopes[i+1]-slopes[i]) for i in range(len(slopes)-1)]
  # find 8-pt min for delta in slope
  mins = [sum(delta[i:int(i)+7]) for i in range(len(delta)-8)]
  possible_pts = [i for i in mins]
  possible_pts.sort()
  for p in possible_pts:
    check = np.mean(slopes[mins.index(p):mins.index(p)+7])
    if check > 0:
      start_pt = mins.index(p)
      break
  frac_dim = np.mean(slopes[start_pt:start_pt+7])
  if where == True:
    return frac_dim, start_pt+3
  return frac_dim



# Fractal dimension
def fractal_dimension(geo, where=False):
  """
  Calculate fractal dimension.
  """
  pts = [nodex(n) for n in geo.nodes]
  if len(pts) > 10000:
    div = int(len(pts)/10000)
    pts = pts[::div] # downsample to ~ 10,000 pts for time
  maxm = farthest_pt(pts)
  bin_e = np.linspace(0., maxm, 1001)
  dists = []
  for i in pts:
    if pts.index(i)%1000==0:
      print('%i (of %i) centers analyzed' %(pts.index(i), len(pts)))
    for j in pts:
      dists.append(dist3(i,j))
  # now all dists are calculated
  hist, bin_e = np.histogram(dists, bin_e)
  #spheres = element_histogram(dists, bin_e)
  bins = [(bin_e[i]+bin_e[i+1])/2 for i in range(len(bin_e)-1)]
  masses = [np.sqrt((bins[i]**2 * hist[i]) / np.sqrt(hist[i])) 
            for i in range(len(bins))] # radius of gyration
  # if hist[i] == 0 this gives nan, clean these date
  bad_inds = [i for i in range(len(masses)) if str(masses[i])=='nan']
  masses = [masses[i] for i in range(len(masses)) if i not in bad_inds]
  bins = [bins[i] for i in range(len(bins)) if i not in bad_inds]
  # should now be a mass for every bin radius
  frac_dim = pt8_slope(bins, masses, where)
  if type(frac_dim) is list:
    return frac_dim[0], bins, masses, frac_dim[1]
  return [frac_dim, bins, masses]
  


# plot fit
def showFractalDimension(frac_d, bins, masses, fname=None):
  # Show a plot of the fit to the constant slope
  frac_dim, loc = pt8_slope(bins, masses, where=True)
  #if bins2 != bins or masses2 != masses:
  #  print('Discrepency in new elements of bins or masses; using old values')
  ptx, pty = [np.log(bins[0]), np.log(bins[loc]), np.log(bins[-1])], []
  pty.append(np.log(masses[loc])-(np.log(bins[loc])-np.log(bins[0]))*frac_dim)
  pty.append(np.log(masses[loc]))
  pty.append(np.log(masses[loc])+(np.log(bins[-1])-np.log(bins[loc]))*frac_dim)
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.scatter(np.log(bins), np.log(masses), color='b', edgecolor='b',
             alpha=0.2)
  ax.plot(ptx, pty, color='r', linewidth=3, alpha=0.5)
  ax.set_xlabel('Log radius um^x')
  ax.set_ylabel('Log mass')
  if fname is not None:
    ax.set_title('Fractal dimension for %s (%.2f)' %(fname, frac_dim))
  else:
    ax.set_title('Fractal dimension (%.2f)' %frac_dim)
  plt.show()
  return 
  


#######################################################################
# Soma position

# This is kind of a pain, only way I can think is to show the user
# and let them decide; could automate and just take opposite of axon ,
# but the program may find the axon wrong and then don't want to compound
# mistakes; also, GM projects to aln and this won't work for that

def display_simple_neuron(geo):
  # default is to show lots of soma and axon but little of everything else
  pts = [] # populate segments
  for s in geo.segments:
    pts.append(s.coordAt(0.5))
  axons = geo.findAxons() # populate axon
  axs = []
  for a in axons:
    for n in a.nodes:
      axs.append([n.x,n.y,n.z])
  sms = [] # populate soma
  for n in geo.soma.nodes:
    sms.append([n.x,n.y,n.z])
  # plotting shit
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for p in pts:
    ax.scatter(p[0],p[1],p[2], c='b',edgecolor='b', alpha=0.2)
  for p in axs:
    ax.scatter(p[0],p[1],p[2], c='r',edgecolor='r', alpha=0.5)
  for p in sms:
    ax.scatter(p[0],p[1],p[2], c='k', edgecolor='k', alpha=0.5)
  ax.set_xlabel('x axis')
  ax.set_ylabel('y axis')
  ax.set_zlabel('z axis')
  ax.set_aspect('equal')
  plt.show()
  return



def place_soma(geo, stn_val):
  """
  Places the stn in 3-D to find the soma's position. stn_val should be
  a 3-tuple where any None value becomes the mean for that coordinate;
  i.e.: for 836_047 the y and z are none but the x is -100 to place the 
  stn at the end of the neuropil, so stn_val = [-100, None, None]
  This function also assumes an X-Y- major projection! Important for phi.
  """
  # These copied from branch_angles
  def dist3(pt0, pt1):
    if len(pt0) == len(pt1) and len(pt0) == 3:
      return math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
    else:
      print('dimension mismatch')
      print(pt0, pt1)
  #
  def get_angle(pt0, midpt, pt1):
    if pt0 in [midpt, pt1] or pt1 in [midpt, pt0] or midpt in [pt0,pt1]:
      print('Some points are the same!')
      print(pt0, midpt, pt1)
    PT0 = dist3(pt1, midpt)
    PT1 = dist3(pt0, midpt)
    MIDPT = dist3(pt0, pt1)
    try:
      ang = math.acos( (MIDPT**2 - PT1**2 - PT0**2) / (2*PT1*PT0) )
      ang = ang*180/math.pi
    except:
      ang = 'nan'
    return ang
  #
  pts = [s.coordAt(0) for s in geo.segments]
  means = [np.mean([m[0] for m in pts]),
           np.mean([m[1] for m in pts]),
           np.mean([m[2] for m in pts])]
  for s in range(len(stn_val)): # make sure stn_val is kosher
    if stn_val[s] == None or stn_val[s] == 0:
      stn_val[s] = means[s]
  theta_val = stn_val
  phi_val = [means[0],means[1],100] # doesn't really matter for z, just needs be positive
  soma_val = geo.soma.coordAt(0)
  theta = get_angle(theta_val, means, soma_val)
  phi = abs(get_angle(phi_val, means, soma_val))-90. # from X-Y plane
  r = dist3(means, soma_val)
  return 180-theta, phi, r



def plot_soma_positions(arr, types=None):
  """
  Input is a list of 3-tuples [theta (from stn in X-Y), phi (elevation
  from X-Y plane at stn and center of neuropil), r (distance from center 
  of neuropil].
  """
  def polar_to_rect(p):
    x = p[2]*np.cos(p[0]/180*np.pi)
    y = p[2]*np.sin(p[0]/180*np.pi)
    z = p[2]*np.sin(p[1]/180*np.pi)
    return [x,y,z]
  pts = [polar_to_rect(p) for p in arr]
  stn_length = max([p[2] for p in arr])
  # get colors
  colors = ['darkkhaki','royalblue','forestgreen','tomato']
  if types:
    if type(types[0]) is not int:
      names = list(set(types))
      types = [names.index(i) for i in types]
    cols = [colors[i] for i in types]
    patches = []
    for n in range(len(names)): # should not exceed 3 (4 cell types)
      patches.append( mpatches.Patch(color=colors[n], label=names[n]) )
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  # set up the axes
  ax.plot([0,stn_length], [0,0],[0,0], linewidth=5, c='k')
  ax.plot([0,-stn_length*.5],[0,0], [0,0], linewidth=1, c='k')
  ax.plot([0,0],[stn_length*.5,-(stn_length*.5)],[0,0], linewidth=1,c='k')
  ax.text(1.1*stn_length, 0, -20, r'stn', style='italic', fontsize=20)
  stnlen = '%.0f um' %(stn_length)
  ax.text(1.1*stn_length, 0, -50, stnlen, fontsize=15)
  # plot the soma
  for p in pts:
    if types:
      ax.scatter(p[0],p[1],p[2],s=100,c=cols[pts.index(p)], 
                 edgecolor=cols[pts.index(p)])
    else:
      ax.scatter(p[0],p[1],p[2],s=100)
  if types:
    plt.legend(handles=patches, loc='best', fontsize=15)
  ax.set_zlim([-stn_length*.5, stn_length*0.5])
  ax.set_xlim([-stn_length*.25, stn_length])
  ax.set_ylim([-stn_length*.25, stn_length])
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_zticks([])
  plt.axis('off')
  plt.show()
  return

###########################################################################
# How 'branched' is the neuron -- get carrier points

# Could just do tips, but that assumes the trajectory is insignificant


def carrier_spacing(geo, step=20.):
  """
  Create carrier points along every path at step-um intervals.
  """
  # Get the paths
  pDF = PathDistanceFinder(geo, geo.soma)
  tipInds, tipLocs = geo.getTipIndices()
  pathsegs = [seg for seg in geo.segments if seg.filamentIndex in tipInds]
  paths = [pDF.pathTo(pseg, 0.5) for pseg in pathsegs]
  
  # Find out may points and assign them
  pathlengths = [sum([s.length for s in p]) for p in paths]
  num_sp_pts = [int(s/float(step)) for s in pathlengths]
  pts = []
  for p in range(len(paths)):
    if p%100 == 0:
      print('%i / %i paths analyzed' %(p, len(paths)))
    path = paths[p]
    dist = list(np.linspace(0,pathlengths[p], num_sp_pts[p])) # pop(0) this as you go
    dist.pop(0)
    so_far = 0
    
    for seg in path:
      try: 
        if dist[0] > so_far and dist[0] < so_far+seg.length: # Put it in here
          remain = dist[0] - so_far
          this_node = seg.nodes[int(remain/seg.length*len(seg.nodes))]
          pts.append([this_node.x, this_node.y, this_node.z]) # Undershoots, but okay
          dist.pop(0)
        else: # Increment so_far because dist[0] wasn't within this seg
          so_far = so_far + seg.length
      except:
        pass
    
    # Got all the dists for this path
  # Got all the paths
  # Now trim the pts array so it's not redundant
  # return pts
  print('Initially %i points. Removing redundant points...' %len(pts))
  cnt = 0
   
  new_pts = [] # Remove redundant points
  for p in pts:
    if p not in new_pts:
      new_pts.append(p)
  pts = new_pts
  
  while cnt < len(pts): # Remove nearby points
    p = pts[cnt]
    dists = [math.sqrt(sum([(p[i]-pt[i])**2 for i in range(3)]))
                                        for pt in pts]
    popit = [d for d in range(len(dists)) if dists[d] < step/2 and dists[d] != 0]
    new_pts = [pts[i] for i in range(len(pts)) if i not in popit]
    pts, cnt = new_pts, cnt + 1

  print('Finally left with %i points.' %len(pts))
  return pts



def longest_shortest(geo, pts=None):
  """
  Given a geo file +/- carrier points, find the range of total wiring.
  """
  if pts is None:
    pts = carrier_spacing(geo)
  
  gpt = geo.soma.coordAt(0)
  # Longest distance
  maxdist = sum([ math.sqrt(sum([(gpt[i]-pt[i])**2 for i in range(3)])) for pt in pts ])
  
  # Get shortest distance, use networkx minimum spanning tree
  c_pts = {}
  for p in range(len(pts)): # Make the dictionary 
    c_pts[p] = pts[p]
  G = nx.Graph()
  G.add_nodes_from(c_pts)
  
  # Add the edges
  edgelist = []
  for n1 in range(len(pts)):
    for n2 in [i for i in range(len(pts)) if i != n1]:
      edgelist.append([n1, n2])
  G.add_edges_from(edgelist)
  
  # Describe MST
  MST = nx.minimum_spanning_tree(G)
  mindist = 0.
  for edge in MST.edges():
    try:
      mindist = mindist + math.sqrt(sum([(c_pts[edge[0]][i]-c_pts[edge[1]][i])**2 
                                         for i in range(3)]))
    except:
      print(edge)
  
  # Get the current geo wiring
  pDF = PathDistanceFinder(geo, geo.soma)
  tipInds, tipLocs = geo.getTipIndices()
  paths = [s for s in pDF.pathTo(tipInds[l], tipLocs[l]) for l in range(len(tipInds))]
  wiring = sum([b.length for b in geo.branches])
  print('Neuron %s is %.2f long, potential length ranges from %.2f - %.2f'
        %(geo.name, wiring, mindist, maxdist))
  
  return wiring, mindist, maxdist

  

##########################################################################
# geo stuff : tips are closer than expected by chance

def ddist3(pt1, pt2, node1=False, node2=False):
  # Euclidean distance to a point
  if node1:
    pt1 = [pt1.x, pt1.y, pt1.z]
  if node2:
    pt2 = [pt2.x, pt2.y, pt2.z]
  return math.sqrt(sum([(pt2[i]-pt1[i])**2 for i in range(3)]))



def getRank(tipInfo, nodePaths, distBounds):
  """
  This determines the rank of the actual tip compared to nodes of 
  similar euclidean distance. Rank 1 = tip path is longer than all node paths,
  rank 0 = tip path is super short
  """
  print('Calculating ranks ...')
  t_ranks = []
  
  for t in tipInfo:
    # Find the key for this tip and the node-paths; skips the longest tip
    try:
      # loc identifies which nodes to access (nodes with euc dist similar to tip value)
      loc = len([i for i in distBounds if i < t[0]])-1
      # Compare the tip path to the node-paths for that euc distance
      if len(nodePaths[loc]) == 0:
        print('No Nodes match euc tip dist of %.2f' %t[0])
        print('Bounds are: %.4f - %.4f' %(distBounds[loc], distBounds[loc+1]))
      else:
        t_ranks.append(float(len([i for i in nodePaths[loc] if i < t[1]])) /
                     float(len(nodePaths[loc])))
    except:
      print('Looking for %.4f, but no hits: range is %.4f - %.4f'
            %(t[0], min(distBounds), max(distBounds)))
  
  print(' ... done.')
  return t_ranks

  

def tips_path_points(geo):
  """
  This funciton compares tip (path) distances to path distances of points
  at similar Euclidean distances from the soma. 
  If tip paths are closer than test paths, tips are closer than expected
  by chance. tolerance = 5um by default. (needed??)
  """
  # Get path lengths and the euc distances of the relevant tips
  pDF = PathDistanceFinder(geo, geo.soma)
  tipInds, tipLocs = geo.getTipIndices()
  distBounds = []
  # [ [eucDist, pathLength], [], ... ]
  tipInfo = [[ddist3(geo.soma.nodeAt(0), 
              seg.nodeAt(tipLocs[tipInds.index(seg.filamentIndex)]),True,True),
              pDF.distanceTo(seg, tipLocs[tipInds.index(seg.filamentIndex)])]
             for seg in geo.segments if seg.filamentIndex in tipInds]
  
  # Find nodes that match this euc distance and get their path lengths
  _, distBounds = np.histogram([t[0] for t in tipInfo], bins=100)
  nodePaths = {i: [] for i in range(len(distBounds)-1)}
  print('Populating nodePaths dictionary...')
  
  for nod in geo.nodes:
    if geo.nodes.index(nod)%1000==0:
      print('%i / %i nodes examined' %(geo.nodes.index(nod), len(geo.nodes)))
    D = ddist3(geo.soma.nodeAt(0), nod, True, True)
    
    if D > min(distBounds) and D < max(distBounds): # If node is within range, else skip
      loc = len([i for i in distBounds if i < D])-1
      # Append the 'path length' to relevant euc distance list 
      nodePaths[loc].append(pDF.distanceTo(nod.segments[0], 0))
  
  # Now should have all the distances, calculate interval-wise non-parametric p-value
  temp_ranks = getRank(tipInfo, nodePaths, distBounds)
  print('Have %i / %i ranks.' %(len(temp_ranks), len(tipInfo)))
  return temp_ranks
  


def show_tips(geofils, labelsin, plusone=0, switch=True):
  """
  Show the locations of the tips of each geofile.
  """
  if switch:
    for i in range(len(geofils)-1):
      geofils.append(geofils.pop(0))
      labelsin.append(labelsin.pop(0))
  
  # Get the tip coords for each file first
  tipCoords = []
  for geo in geofils:
    tipInds, tipLocs = geo.getTipIndices()
    tNodes = [seg.nodeAt(tipLocs[tipInds.index(seg.filamentIndex)])
              for seg in geo.segments if seg.filamentIndex in tipInds]
    tCoords = [ [n.x, n.y] for n in tNodes ]
    tipCoords.append(tCoords)
  
  # Condition the data
  tipCoords = [t[::int(len(t)/min([len(i) for i in tipCoords])*5)]
               for t in tipCoords]
  
  # Now plot these mofos
  fig = plt.figure()
  sq = int(np.sqrt(len(geofils)+plusone))
  plots = [fig.add_subplot(sq, sq, i) for i in range(len(geofils))]
  colors = ['darkkhaki', 'royalblue', 'forestgreen','tomato']
  L = list(np.unique(labelsin))
  C = [L.index(i) for i in labelsin]
  
  for p in range(len(plots)):
    plots[p].scatter([i[0] for i in tipCoords[p]], 
                     [i[1] for i in tipCoords[p]], color=colors[C[p]],
                     s=10, alpha=0.6, edgecolor=colors[C[p]])
  
  plt.show()
  return







############################################################################
# Radius stuff (Figure 11)

def single_ratios(rlist, skip=2):
  """
  Given a list, this assumes the format is parent1 daugh1a daught1b ...
  parentn daughtna daughtnb and returns a list daught1a/parent1 ...
  daughtnb/parentn. (Exactly 2/3 the size of input.)
  skip = len(str fields) that being each object (file, celltype)
  """
  outlist, count = [], 0
  if type(rlist[3]) is str: # First few items are usually str
    nlist = []
    for r in rlist:
      try:
        nlist.append(float(r))
      except:
        if str(r) == 'x':
          pass
        else:
          nlist.append(r)
    rlist = nlist
  rlist = [n for n in rlist if n != 0]
  for i in range(int((len(rlist)-skip)/3)):
    outlist.append(rlist[i*3+1+skip]/rlist[i*3+skip])
    outlist.append(rlist[i*3+2+skip]/rlist[i*3+skip])
  return outlist
  
  

def combined_ratios(rlist, skip=2):
  """
  Same as single_ratios but sums daughters; returns (daughters1a+b)/parent1.
  (Exactly 1/3 the size of input.)
  """
  outlist = []
  if type(rlist[3]) is str:
    nlist = []
    for r in rlist:
      try:
        nlist.append(float(r))
      except:
        if str(r) == 'x':
          pass
        else:
          nlist.append(r)
    rlist = nlist
  rlist = [n for n in rlist if n != 0]
  for i in range(int((len(rlist)-skip)/3)):
    outlist.append((rlist[i*3+1+skip]+rlist[i*3+2+skip])/rlist[i*3+skip])
  return outlist



def rall_analysis(lol):
  """
  Each list is [filename, celltype, p1, d1a, d1b, p2, d2a, d2b ... ]
  """
  ralls = []
  for l in lol:
    print('Calculating rall exponents for %s' %l[0])
    ralls.append([l[0]]) # Add filename
    idx = lol.index(l) # Current cell index
    ralls[idx].append(l[1]) # Add celltype
    startlist = [i*3+2 for i in list(range(int((len(l)-2)/3)))]
  
    for s in startlist:
      for rad in [1,2]: # Two 
        # Send the par, [d1, d2] pairs to fit_rall_exp
        ralls[idx].append(fit_rall_exp(l[s], l[s+rad])) 
  ralls = [[i for i in k if i is not None] for k in ralls]
  print('All rall exponents returned')
  return ralls



def rall_analysis_dict(aDict, whichKeys, newKeys=None, norm=False):
  """
  This acts on an entire dict with the given keys; new keys can be passed,
  otherwise they default to +'_rall'.
  """
  bDict = copy.deepcopy(aDict)
  if newKeys is None:
    newKeys = [s+'_rall' for s in whichKeys]
  elif len(newKeys) != len(whichKeys):
    print('whichKeys (len %i) and newKeys (%i) should be same length' 
          %(len(whichKeys), len(newKeys)))
    return
  
  # For each key, compute the ralls and add it
  for k in range(len(whichKeys)):
    bDict[newKeys[k]] = [[] for i in bDict[whichKeys[k]]]
    for b in range(len(bDict[newKeys[k]])): # For each cell in this key
      if len(bDict[whichKeys[k]][b]) % 3 == 0: # If divisible by 3 (par, d1, d2)...
        num_trips = int(len(bDict[whichKeys[k]][b])/3)
        for s in range(num_trips): # For each triplet
          
          bDict[newKeys[k]][b].append(fit_rall_exp(bDict[whichKeys[k]][b][s*3], # par
                                                   [bDict[whichKeys[k]][b][s*3+1], # d1
                                                   bDict[whichKeys[k]][b][s*3+2]],
                                                   norm)) # d2
  
  return bDict


  
# FIX RALL!!!!!!!!!!!!!!!!!!!!!!!!11

def fit_rall_exp(par, daugs, norm=False): # This is the Rall power calculator
  """
  (daughter1^x + daughter2^x + ...) = parent^x
  """
  if type(daugs) is not list:
    print('Daugters should be list, got %s instead' %str(type(daugs)))
    print(par, daugs)
    return
  
  if norm:
    mmax = max([par, daugs[0], daugs[1]])
    if mmax == 0:
      print('Parent: %.4f, Daugs: %.4f, %.4f; divide by zero' %(par, daugs[0], daugs[1]))
      return None
    par = par/mmax
    daugs = [d/mmax for d in daugs]
  
  try:
    exps = np.linspace(0.001, 50, 10000)
    remain = [(sum([d**x for d in daugs]) - par**x) for x in exps]
    exps = [exps[i] for i in range(len(remain)) if remain[i]>0]
    remain = [remain[i] for i in range(len(remain)) if remain[i]>0]
    return exps[remain.index(min(remain))]
    # return math.log(par, daug)
  
  except:
    print('Wrong type: par: %s, daug: %s' %(str(type(par)), str(type(daugs))))
    print(par, daugs)
    return


def plot_rall(par_daugs=None, rang=[0.001,50], steps=10000, bench=1.5,
              xlim=[0,5], ylim=[0,5], showleg=False):
  """
  Examine a rall fit if numbers seem funky. Can submit one value triplet
  ([par, d1, d2]) or several ([[p,d1,d2],[p,d1,d2],[p,d1,d2]]).
  
  """
  # Check the data first; if none are supplied, use defaults.
  if par_daugs is None:
    par_daugs = [ [1.0, 0.5, 0.5], [1.0, 1.0, 0.5], [1.0, 1.5, 0.5],
                  [1.0, 1.5, 1.5], [1.0, 0.6, 0.6], [1.0, 0.7, 0.7],
                  [1.0, 0.8, 0.8], [1.0, 0.9, 0.9] ]
  if type(par_daugs) is list:
    if type(par_daugs[0]) is not list: # Assume a single-triplet entry
      par_daugs = [par_daugs]
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  for entry in par_daugs:
    par, daugs = entry[0], entry[1:3]
    xs = np.linspace(0.001, 50, 10000)
    ys = [(sum([d**x for d in daugs]) - par**x) for x in xs]
    xs = [xs[i] for i in range(len(ys)) if ys[i]>0]
    ys = [ys[i] for i in range(len(ys)) if ys[i]>0]
    lab = 'p=%.1f, d=%.1f,%.1f' %(par, daugs[0], daugs[1])
    ax.plot(xs, ys, label=lab, c=np.random.rand(3))
  
  if bench is not None:
    try:
      ax.plot([1.5,1.5], [ylim[0], ylim[1]], '--', c='purple', linewidth=2, alpha=0.5)
    except:
      print('bench must be int, floar or None; got' + type(bench))
  
  ax.set_xlabel('Exponent')
  ax.set_ylabel('(D1^x+D2^x)-P^x')
  ax.set_xlim(xlim)
  ax.set_ylim(ylim)
  for pos in ['top', 'right']: # Also hide these borders for all plots
    ax.spines[pos].set_visible(False)
  
  if showleg:
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
  plt.show()
  return


#def fit_rall_exp_list(parlist, dauglist): # Deprecated #
#  """
#  Given a list of ratios this fits the rall exponents (3/2 is 'ideal').
#  """
#  if type(parlist) is not list and type(dauglist) is list:
#    return fit_rall_exp(parlist, dauglist)
#  assert len(parlist) == len(dauglist), "input1 (length: %i) and %input2 (length: %i)" \
#                                        %(len(parlist), len(dauglist))
#  return [math.log(parlist[i], dauglist[i]) for i in range(len(parlist))]




def pixel_div(x, y=None):
  # Get the 'average' pixel size to divide out.
  if y is None:
    y = x
  return 0.5*(np.mean([x,y]) + np.sqrt(2*x*y))



def div_radius(tips, divs, hand=False):
  # Divide non-string numbers by the divisor, len(tips)==len(divs)
  for t in range(len(tips)):
    for i in range(len(tips[t])):
      if hand is True:
        if type(tips[t][i]) is not str and tips[t][i] > 1: # Not binary
          tips[t][i] = tips[t][i]/divs[t]
      else: # Not hand
        if type(tips[t][i]) is not str:
          tips[t][i] = tips[t][i]/divs[t]
  return tips



def plot_all_rall(ralldict, bounds=[0,10]):
  """
  Given a dict of rall objects this makes all the plots for figure 3 (E-H).
  Keys of rall dict must include: prim_sec_raw_rall, sec_tert_raw_rall,
  tips_raw_rall init_rall, and cellTypes.
  """
  targets = ['prim_sec_raw_rall', 'sec_tert_raw_rall', 'tips_raw_rall']
  for att in targets:
    if att not in ralldict.keys():
      print('Dictionary is missing %s!' %att)
      return None
  
  # This function calls hori_scatter from pretty_plot.py
  figs = []
  for target in targets:
    labs, vals = condition_by_name(ralldict['cellTypes'], ralldict[target])
    figs.append(hori_scatter(vals, labs, axes=['','Rall power'], 
                             title=target, bounds=bounds, switch=True,
                             shade=False, fill=True, bench=1.5, retplot=True))
  
  rlabs, rvals = condition_by_name(ralldict['cellTypes'],ralldict['init_rall'])
  rlabs = [rlabs[i] for i in [0,4,8,12]]
  rvals = [[rvals[u*4+i] for i in range(4)] for u in range(4)]
  print(rvals)
  rfig = hori_scatter(rvals, rlabs, axes=['','Rall power'], size=25,
                      title='init_rall', bounds=bounds, showmean=False,
                      shade=False, fill=True, bench=1.5, retplot=True)
  
  plt.show()
  return

  


#########################################################################
# Quadratic taper by path for estimated tip 
# See bottom of quaddiameter.py


def collapse_list_to_dict(lols, keys):
  """
  Each lol in lols should be loaded as a csv (or similar), i.e.:
  lol1 = get_csv(), ...
  Keys should be the data name of each lol loaded by csv.
  newdict = collapse_list_to_dist([lol1, lol2, ...], ['tips', 'prim', ...])
  """
  newdict = {}
  newdict['files'] = [g[0] for g in lols[0]]
  newdict['cellTypes'] = [g[1] for g in lols[0]]
  
  for lol in lols: # For each list
    # Assume file and celltype are first, and len(keys)==len(lols)
    newdict[keys[lols.index(lol)]] = [[] for i in range(len(lol))]
    for l in lol: # For each cell, match its placement
      newdict[keys[lols.index(lol)]][newdict['files'].index(l[0])] = l[2:]
  
  return newdict
  

    

def get_csv(csvfile):
  # Return list of csv row elements
  arr = []
  with open(csvfile, 'r') as cfile: # Assumes text, not binary
    creader = csv.reader(cfile) # Assumes csv
    for row in creader:
      nrow = []
      for i in row:
        try:
          nrow.append(float(i))
        except:
          nrow.append(i)
      arr.append(nrow)
  return arr



def save_csv(tips, outfile, cols=None, rows=None):
  # Save a tips-like list of lists as a .csv
  if rows is not None:
    if cols is None:
      cols = range(max([len(t) for t in tips]))
    get = [[i] for i in rows]
    for k in cols:
      get[0].append(k)
    for t in range(len(tips)): # len(tips) = len(get)-1
      for c in tips[t]:
        get[t+1].append(c)
  else:
    get = tips
  with open(outfile, 'w') as fOut:
    for g in get:
      fOut.write(','.join([str(i) for i in g]))
      fOut.write('\n')
  print('%s written.' %outfile)
  return



############
if __name__ == "__main__":
  print("Module needs to be used interactively.")
































