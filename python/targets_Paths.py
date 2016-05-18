# targets_Paths.py

"""
Most of the path stuff, including finding targets, etc.
This mostly replaces adriane2.py.
"""



# Need some imports for pDF and plt and such
from pyramidal_readExportedGeometry import * 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import scipy
import scipy.stats as stats
import scipy.optimize as op
if sys.version_info != 2:
  print('Will not be able to use MCMC or maximum likelihood estimator')
else: 
  import pymc as pm
import seaborn as sns
sns.set_style('white')
from PIL import Image



def pyplot_defaults():
  """ After calling pandas or seaborn, the default plt can be messed up.
  If this causes inconsistencies, can call this to revert back to defaults. """
  import matplotlib as mpl
  mpl.rcParams.update(mpl.rcParamsDefault)
  return


######################################################################
# Helper functions

def dist3(pt0, pt1):
  """
  In theory this can handle node objects.
  """
  if type(pt0) is not list:
    try:
      pt0 = [pt0.x, pt0.y, pt0.z] # Try for a node
    except:
      print('what do I do with ' + pt0)
  if type(pt1) is not list:
    try:
      pt1 = [pt1.x, pt1.y, pt1.z] # Try for a node
    except:
      print('what do I do with ' + pt1)
  try:
    return math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
  except:
    print('dimension mismatch')
    print(pt0, pt1)
  return None


def xy_vs_z(pt0, pt1):
  if type(pt0) is not list:
    try:
      pt0 = [pt0.x, pt0.y, pt0.z] # Try for a node
    except:
      print('what do I do with ' + pt0)
  if type(pt1) is not list:
    try:
      pt1 = [pt1.x, pt1.y, pt1.z] # Try for a node
    except:
      print('what do I do with ' + pt1)
  xy = math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(2)]))
  z = math.sqrt((pt0[2]-pt1[2])**2)
  if xy > z:
    return 'xy'
  else:
    return 'z'


def get_colors():
  # Colors with max contrast
  contcols = ['lightskyblue', 'brown', 'orange', 'springgreen',
              'fuchsia', 'tomato', 'gold', 'indigo',
              'darkslateblue', 'black', 'darkgreen', 'aqua',
              'darkorchid', 'grey', 'salmon', 'plum',
              'coral', 'sienna', 'darkkhaki', 'yellowgreen',
              'deeppink', 'ivory', 'orchid', 'lightsteelblue']
  return concols



########################################################################
# Loading and reading neurons

def findBounds(geo):
  """
  Find the bounds of the neuron object.
  """
  bounds = {'minx': np.inf, 'maxx': 0,
            'miny': np.inf, 'maxy': 0,
            'minz': np.inf, 'maxz': 0}
  # Find the bounds of the geo object
  for n in geo.nodes:
    if n.x > bounds['maxx']: # Xs
      bounds['maxx'] = n.x
    if n.x < bounds['minx']:
      bounds['minx'] = n.x
    if n.x > bounds['maxy']: # Ys
      bounds['maxy'] = n.y
    if n.x < bounds['miny']:
      bounds['miny'] = n.y
    if n.x > bounds['maxz']: # Zs
      bounds['maxz'] = n.z
    if n.x < bounds['minz']:
      bounds['minz'] = n.z
  #print(bounds)
  print('Min X: %.2f, Max X: %.2f,\nMin Y: %.2f, Max Y: %2.f, \nMin Z: %.2f, Max Z: %.2f'
        %(bounds['minx'], bounds['maxx'], bounds['miny'],
         bounds['maxy'], bounds['minz'], bounds['maxz']))
  return



def loadTargets(infile, kind='xml', voxel=[.732,.732,.49], order=None,
                paired=False):
  """
  Load the targets from an xml file (as from knossos).
  """
  # Get the tuple from the line
  def get_xyz(line):
    splitLine = line.split('"')
    temp = []
    for u in [5, 7, 9]:
      try:
        temp.append(float(splitLine[u]))
        #print(float(splitLine[u]))
      except:
        print('problem converting to float' + splitLine[u])
    return temp
  
  targets = []
  with open(infile, 'r') as fIn:
    for line in fIn:
      splitLine = line.split('"')
      if len(splitLine) > 15:
        # Format is:
        # <node id="1" radius="3.8042" x="458" y="290" z="120" inVp="0" inMag="1" time="147586"/>
        if splitLine[0].split(None)[0].split('<')[1] == 'node':
          targets.append(get_xyz(line))
  
  for t in range(len(targets)):
    targets[t] = [targets[t][i]*voxel[i] for i in range(3)]
  if order is not None: # Sort these in "adriane's order"
    if len(order) == len(targets):
      new_targs = [targets[i] for i in order]
    return new_targs
  return targets




########################################################################
# Paths --- with and without targets

"""
This function is called extensively.
"""

def getPaths(geo, targets, show='solo'):
  """
  Return the distance from the soma to each of the uncaging targets.
  If the targets were ordered in loadTargets, their order is the key.
  show = 'solo' -- individual; show='all' -- overlapped
  """
  pDF = PathDistanceFinder(geo, geo.soma)
  targ_dict = {i: {'location': targets[i], 'pathLength': None,
                   'closestFil': None} for i in range(len(targets))}
  seglist = []
  
  for t in range(len(targets)):
    # Find the closest segment first
    mindist, minnode = np.inf, None
    for node in geo.nodes:
      if dist3(node, targets[t]) < mindist:
        mindist = dist3(node, targets[t])
        minnode = node
    
    # Found the min node, assign its segment to targ_dict
    seglist.append(minnode.segments[0])
    targ_dict[t]['closestFil'] = seglist[-1].filamentIndex
    targ_dict[t]['pathLength'] = pDF.distanceTo(seglist[-1]) # Just take the first one
  # Got all the targets assigned
  
  #### Plotting stuff! ####
  if show is 'solo': 
    branchpts = []
    # Plot by branches
    for b in geo.branches: 
      branchpts.append([[n.x, n.y] for n in b.nodes])
    if len(targets) <= 8:
      fig = plt.figure(figsize=(6,4), dpi=150)
    else:
      fig = plt.figure(figsize=(6,8), dpi=150)
    for t in targ_dict.keys():
      if len(targets) <= 8:
        ax = fig.add_subplot(2,4,t+1) # Works for <= 8 targets, assumes keys are ints
      elif len(targets) > 8 and len(targets) <= 16:
        ax = fig.add_subplot(4,4,t+1)
      else:
        ax = fig.add_subplot(5,5,t+1)
      for b in branchpts: # Plot the background skeleton first
        for s in range(len(b)-1):
          plt.plot([b[s][0], b[s+1][0]],
                   [b[s][1], b[s+1][1]], color='black', alpha=0.3)
      path = pDF.pathTo(seglist[t], 0)
      col = contcols[t-1]
      for p in range(len(path)-1):
        c0, c1 = path[p].coordAt(0), path[p+1].coordAt(0)
        plt.plot([c0[0], c1[0]], [c0[1], c1[1]], linewidth=3,
                 color=col, alpha=0.7)
      plt.plot(targets[t][0], targets[t][1], color=col, markeredgecolor='none',
               marker='o')
      ax.set_title('Target %i: %.2f um' %(t+1, targ_dict[t]['pathLength']))
      if t != 0:
        for pos in ['top', 'right', 'bottom', 'left']:
          ax.spines[pos].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    plt.show()
  
  if show is 'all':
    branchpts = []
    # Plot by branches
    for b in geo.branches: 
      branchpts.append([[n.x, n.y] for n in b.nodes])
    #
    plt.figure(figsize=(5,6), dpi=150) #ax2 = plt.subplot2grid((1,3), (0,2))
    ax1 = plt.subplot2grid((1,5), (0,0), colspan=4)
    ax2 = plt.subplot2grid((1,5), (0,4))
    for b in branchpts: # Plot the background skeleton first
      for s in range(len(b)-1):
        ax1.plot([b[s][0], b[s+1][0]],
                 [b[s][1], b[s+1][1]], color='black', alpha=0.3)
    for t in targ_dict.keys():
      path = pDF.pathTo(seglist[t], 0)
      col = contcols[t-1]
      #for p in range(len(path)-1): ### Show paths! ###
      #  c0, c1 = path[p].coordAt(0), path[p+1].coordAt(0)
      #  ax1.plot([c0[0], c1[0]], [c0[1], c1[1]], '--', linewidth=3,
      #           color=col, alpha=0.7)
      ax1.plot(targets[t][0], targets[t][1], color=col, markeredgecolor='none',
               marker='o')
      ax1.text(targets[t][0], targets[t][1]+20., t, color=col, fontweight='bold')
      ax2.text(0.,float(t)/len(targ_dict.keys())+0.03, 'Target %i: %.0f um' %(t, targ_dict[t]['pathLength']), 
               fontsize=5, color=col, fontweight='bold')
    plt.show()
  
  return targ_dict



def pathsByLength(geo, targ_dict, background=True):
  """
  Overlay paths all on one tree colored by path length.
  """
  lengths = [targ_dict[s]['pathLength'] for s in targ_dict.keys()]
  plt.figure(figsize=(4,6))
  
  if background:
    branchpts = []
    # Plot by branches
    for b in geo.branches: 
      branchpts.append([[n.x, n.y] for n in b.nodes])
    for b in branchpts: # Plot the background skeleton first
        for s in range(len(b)-1):
          plt.plot([b[s][0], b[s+1][0]],
                   [b[s][1], b[s+1][1]], color='black', alpha=0.3)
  
  # Plot each thing
  pDF = PathDistanceFinder(geo, geo.soma)
  for t in targ_dict.keys():
    targ = [seg for seg in geo.segments if 
              seg.filamentIndex==targ_dict[t]['closestFil']][0]
    path = pDF.pathTo(targ)    # This makes sure the colors range the whole spectrm
    pcolor = (targ_dict[t]['pathLength']-min(lengths))/(max(lengths)-min(lengths))
    for p in range(len(path)-1):
      c0, c1 = path[p].coordAt(0), path[p+1].coordAt(0)
      plt.plot([c0[0], c1[0]], [c0[1], c1[1]], linewidth=3,
               color=cm.hsv(pcolor), alpha=0.7)
  plt.show()
  



"""
*pretty_skeleton* from pretty_plot.py to plot the skeletons.
"""



def targetDistances(geo, targs, show=False, save=False):
  """
  Calculate filament distances between targets.
  """
  # Get the target nodes and segments first
  def segfil(geo, filnum):
    return [seg for seg in geo.segments if seg.filamentIndex == filnum][0]
  #
  targ_dict = getPaths(geo, targs, show=False) # Do not show paths
  #targ_list = [[k, targ_dict[k]['location']] for k in targ_dict.keys()]
  #print([i[0] for i in targ_list])
  targ_fils = {k: targ_dict[k]['closestFil'] for k in targ_dict.keys()}
  print(targ_fils)
  targ_segs = [segfil(geo, targ_dict[i]['closestFil']) for i in 
               range(len(targ_dict.keys()))]
  # Calculate filament distances between these
  dists = []
  for refseg in targ_segs:
    # Avg seg length is like 5 um, so doens't matter which end we select
    pDF = PathDistanceFinder(geo, refseg, 0)
    dists.append([pDF.distanceTo(seg, 0) for seg in targ_segs])
    #print('For target %i:' %targ_segs.index(refseg))
    #print('0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10')
    #print('\t'.join(['%.0f' %i for i in dists[targ_segs.index(refseg)]]))
  dists = np.array(dists)
  
  ### Plotting ###
  if show or save:
    plt.figure(figsize=(7,4), dpi=150)
    #
    # First heatmap
    mask = np.tri(dists.shape[0], k=-1) # Mask out lower triangle
    dists = np.ma.array(dists, mask=mask)
    ax1 = plt.subplot2grid((1,3), (0,0), colspan=2) # pos=(0,0), span 2 cols
    heatmap = ax1.pcolor(dists, cmap=plt.cm.get_cmap('jet'))
    labels = range(dists.shape[0])
    #dists = dists[:-1,1:]
    #ax1.imshow(dists, interpolation='nearest', cmap=plt.cm.get_cmap('jet', 10))
    for y in list(range(dists.shape[0])):
      for x in list(range(dists.shape[1])):
        if dists[y,x] != 'nan' and dists[y,x] != 0:
          ax1.text(x + 0.5, y + 0.5, '%.0f' %dists[y,x],
                   horizontalalignment='center', # right
                   verticalalignment='center',) # bottom
    #
    plt.colorbar(heatmap)
    ax1.set_title('%s' %geo.name)
    for axis in [ax1.xaxis, ax1.yaxis]: # Center axis ticks
      axis.set(ticks=np.arange(0.5, len(labels)), ticklabels=labels)
    ax1.yaxis.tick_right() # Move ticks to right
    for pos in ['left', 'top']: # Hide some axes
      ax1.spines[pos].set_visible(False)
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_xticklabels(range(dists.shape[0]), ha='center', minor=False)
    ax1.set_yticklabels(range(dists.shape[1]), ha='center', minor=False)
    ax1.patch.set_visible(False) # Remove ugly whitespace
    #
    # Then skeleton
    ax2 = plt.subplot2grid((1,3), (0,2))
    branchpts = []
    for b in geo.branches: 
      branchpts.append([[n.x, n.y] for n in b.nodes])
    #
    for b in branchpts: # Plot the background skeleton first
      for s in range(len(b)-1):
        ax2.plot([b[s][0], b[s+1][0]],
                 [b[s][1], b[s+1][1]], color='black', alpha=0.3)
    #
    # Label the points, remove the labels
    for t in targ_dict.keys():
      ax2.text(targ_dict[t]['location'][0], targ_dict[t]['location'][1],
               '%i' %t, fontsize=10, color='red')
    for pos in ['bottom', 'left', 'top', 'right']:
      ax2.spines[pos].set_visible(False)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
  if show:
    plt.show()
  if save:
    plt.savefig('%s-tips.eps' %geo.name, bbox_inches=0)
    plt.close()
  #
  # Print these things
  print(targ_dict)
  return dists



#
def showTipPaths(geo, targs, tips=None):
  """
  Show a few paths between targets(?).
  """
  def segfil(geo, filnum):
    return [seg for seg in geo.segments if seg.filamentIndex == filnum][0]
  #
  targ_dict = getPaths(geo, targs, show=False) # Do not show paths
  #targ_list = [[k, targ_dict[k]['location']] for k in targ_dict.keys()]
  #print([i[0] for i in targ_list])
  targ_fils = {k: targ_dict[k]['closestFil'] for k in targ_dict.keys()}
  print(targ_fils)
  targ_segs = [segfil(geo, targ_dict[i]['closestFil']) for i in 
               range(len(targ_dict.keys()))] #seg for seg in geo.segments if seg.filamentIndex in targ_fils]
  # Calculate filament distances between these
  dists, paths = [], []
  for refseg in targ_segs:
    # Avg seg length is like 5 um, so doens't matter which end we select
    pDF = PathDistanceFinder(geo, refseg, 0)
    dists.append([pDF.distanceTo(seg, 0) for seg in targ_segs])
    paths.append([pDF.pathTo(seg, 0) for seg in targ_segs])
    #print('For target %i:' %targ_segs.index(refseg))
    #print('0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10')
    #print('\t'.join(['%.0f' %i for i in dists[targ_segs.index(refseg)]]))
  dists = np.array(dists)  
  # Plot theses paths
  plot_these, plot_dists = [], []
  if tips is None:
    for y in range(3):
      temp = [i if i != 0 else np.mean(dists[y]) for i in dists[y]]
      plot_these.append(paths[y][temp.index(min(temp))])
      plot_dists.append(temp[temp.index(min(temp))])
      plot_these.append(paths[y][temp.index(max(temp))])
      plot_dists.append(temp[temp.index(max(temp))])
  elif type(tips) is list and len(tips[0]) == 2: # [ [t1,t2], [t2,t3] ]
    plot_these = [paths[t[0]][t[1]] for t in tips]
    plot_dists = [dists[t[0]][t[1]] for t in tips]
  
  fig = plt.figure()
  ax2 = fig.add_subplot(111)  
  #ax2 = plt.subplot2grid((1,1), (0,0), colspan=2)
  branchpts = []
  for b in geo.branches: 
    branchpts.append([[n.x, n.y] for n in b.nodes])
  #
  for b in branchpts: # Plot the background skeleton first
    for s in range(len(b)-1):
      ax2.plot([b[s][0], b[s+1][0]],
               [b[s][1], b[s+1][1]], color='black', alpha=0.3)
  # Add some paths
  patches = []
  for path in plot_these:
    col = contcols[plot_these.index(path)]
    for p in range(len(path)-1): ### Show paths! ###
      c0, c1 = path[p].coordAt(0), path[p+1].coordAt(0)
      ax2.plot([c0[0], c1[0]], [c0[1], c1[1]], linewidth=3,
               color=col, alpha=0.7,)
    print('%s: %.0f' %(contcols[plot_these.index(path)],
                       plot_dists[plot_these.index(path)]))
    # Add the patch for this color
    patches.append(mpatches.Patch(color=col,
                   label='%.0f' %plot_dists[plot_these.index(path)]))
  # Label the points, remove the labels
  for t in targ_dict.keys():
    ax2.text(targ_dict[t]['location'][0], targ_dict[t]['location'][1],
             '%i' %t, fontsize=10, color='red')
  for pos in ['bottom', 'left', 'top', 'right']:
    ax2.spines[pos].set_visible(False)
  ax2.get_xaxis().set_visible(False)
  ax2.get_yaxis().set_visible(False)
  plt.legend(handles=patches)
  ax2.set_title('%s' %geo.name)
  plt.show()
  return



"""
This function is used a lot.
"""

def getDiameters(geo, diamfile, voxel=[0.732, 0.732, 0.488], report=False,
                 tofile=None, show='all'):
  """
  Return the length of each pair of nodes and give the midpoint (to
  help assign it to an uncaging position). Returns targ_dict.
  """
  if type(geo) is str:
    geo = demoReadsilent(geo)
  
  nodes, edges = {}, []
  with open(diamfile, 'r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        qsplit = line.split('"')
        if splitLine[0] == '<node': # Found a new node
              #0      1   2      3   4   5   6   7   8  9  10     11...
          #<node id="55"radius="1.5" x="493" y="160" z="43" inVp="0" inMag="1" time="11006907"/>
          nodes[int(qsplit[1])] = [float(qsplit[u]) for u in [5, 7, 9]]
        
        elif splitLine[0] == '<edge':
          #         0       1    2     3   4
          # <edge source="44" target="45"/>
          edges.append([int(qsplit[u]) for u in [1, 3]])
        else:
          pass
  
  # Re-scale the coordinates, find the distances and midpts
  for k in nodes.keys():
    nodes[k] = [nodes[k][i]*voxel[i] for i in range(3)]
  dists, midpts = [], []
  for e in edges:
    t_d = dist3(nodes[e[0]], nodes[e[1]])
    dists.append([t_d, xy_vs_z(nodes[e[0]], nodes[e[1]])])
    t_pt = [(nodes[e[0]][i]+nodes[e[1]][i])/2. for i in range(3)]
    midpts.append(t_pt)
    #print(t_pt)
    #print('    diameter = %.3f' %t_d)
  
  print(len(edges), edges)
  print(dists)
  # Assume two targets, link together closest ones
  targets = {}
  used = []
  for mid0 in range(len(midpts)):
    mind, idx = np.inf, None
    for mid1 in range(len(midpts)):
      D = dist3(midpts[mid0], midpts[mid1])
      if D < mind and mid0 != mid1 and mid1 not in used:
        mind = D
        loc = mid1
    if mid0 not in used and loc not in used:
      if dists[mid0][1] == dists[loc][1]: # would re-write xy or z
        print('Tried to overwrite %s with %s' 
              %(dists[mid0][1], dists[loc][1]))
        print(dists[mid0][0], dists[loc][0])
        if dists[mid0][0] >= dists[loc][0]: # Larger (or same)
          dists[loc][1] = ['xy', 'z', 'xy'][['xy', 'z', 'xy'].index(dists[loc][1])+1]
          print('Replaced(v1) %s with %s' %(dists[mid0][1], dists[loc][1]))
        else:
          repeated = ['xy', 'z', 'xy'][['xy', 'z', 'xy'].index(dists[loc][1])]
          replace = ['xy', 'z', 'xy'][['xy', 'z', 'xy'].index(dists[loc][1])+1]
          print(repeated, replace)
          dists[mid0][1] = replace
          print('Replaced(v2) %s with %s' %(dists[loc][1], dists[mid0][1]))
      targets[len(targets.keys())] = {'midpt': midpts[mid0], dists[mid0][1]: dists[mid0][0],
                                     dists[loc][1]: dists[loc][0]}
      used.append(mid0)
      used.append(loc) # Add to used
  
  if report:
    print(targets)
  pops = [] # Assign as targets
  for k in targets.keys():
    for op in ['z', 'xy', 'midpt']: # Make sure nothing is missing
      if op not in targets[k].keys():
        print('%s missing!' %op)
        pops.append(k)
  for p in pops:
    targets.pop(p)
  targets_ = [targets[k]['midpt'] for k in range(len(targets.keys()))]
  
  targ_dict = getPaths(geo, targets_, show=show) # Do not show paths
  for t in targ_dict.keys():
    targ_dict[t]['branchOrder'] = [seg.branchOrder for seg in geo.segments
                                   if seg.filamentIndex == targ_dict[t]['closestFil']][0]
    targ_dict[t]['xy_diam'] = targets[t]['xy']
    targ_dict[t]['z_diam'] = targets[t]['z']
  if report: ## Report findings
    for k in targets.keys():
      print('Target %i at (%.3f, %.3f, %.3f)' %(k, targ_dict[k]['location'][0],
             targ_dict[k]['location'][2], targ_dict[k]['location'][2]))
      print('Distance from soma: %.3f um, diameter: %.3f x %.3f um' 
             %(targ_dict[k]['pathLength'], targets[k]['xy'], targets[k]['z']))
      print('   Cross sectional area: %.3f um' %(np.pi*targets[k]['xy']*targets[k]['z']))
  
  # Save to file
  if tofile is not None:
    with open(tofile, 'w') as fOut:
      fOut.write(','.join(['target','x','y','z','xy','z','area','branchOrder']))
      fOut.write('\n')
      for k in targ_dict.keys():
        temp_ = [k, targ_dict[k]['location'][0], targ_dict[k]['location'][1],
                 targ_dict[k]['location'][2], targ_dict[k]['pathLength'],
                 targets[k]['xy'], targets[k]['z'],
                 np.pi*targets[k]['xy']*targets[k]['z'],
                 targ_dict[k]['branchOrder']]
        temp_ = [str(i) for i in temp_]
        fOut.write(','.join(temp_))
        fOut.write('\n')
    print('%s written' %tofile)
  
  return targ_dict



########################################################################
# Axon identification, some analysis and some axon plotting


def idAxons(geo, retax=False, show=True):
  """
  """
  btags = []
  for b in geo.branches:
    for k in b.tags:
      if k.split('_')[0] != 'filament':
        btags.append([geo.branches.index(b),k])
  if show:
    for b in btags:
      if b[1] == 'Axon':
        c0, c1 = geo.branches[b[0]].coordAt(0), geo.branches[b[0]].coordAt(1)
        thiscol = color=np.random.random(3)
        plt.plot([c0[0], c1[0]], [c0[1], c1[1]], color=thiscol, lw=3, alpha=0.6)
        plt.text(c1[0], c1[1], b[0], fontsize=10, color=thiscol)
    for s in geo.segments:
      c_ = s.coordAt(0)
      plt.plot(c_[0], c_[1], '.', color='gray', alpha=0.4)
    plt.show()
  if retax:
    return [j[0] for j in btags]
  return




def idMoreAxons(geo, N=2, show=True):
  """
  Helps identify axons and the branches leading up to axons. N is the number
  of upstream neighbors to be co-identified.
  """
  def upstream(geo, ref_ind, test_ind): # Is test closer to soma than ref? T/F
    pDF, dists = PathDistanceFinder(geo, geo.soma), []
    for bra in [ref_ind, test_ind]: # Get distances for random segments
      segs = [sg for sg in geo.segments if sg.name in geo.branches[bra].tags]
      dists.append(pDF.distanceTo(segs[np.random.randint(0,len(segs))], 0.0))
    if dists[0] < dists[1]: # ref < test = downstream
      return False
    elif dists[0] > dists[1]: # ref > test = upstream (add)
      return True
    elif dists[0] == dists[1]:
      print('Equal test and reference distances: %.2f and %.2f' %(dists[1], dists[0]))
      return True
    else:
      print('Something screwy happened'); return None
    
  # Add the upstream axons
  axs = idAxons(geo, retax=True, show=False) # Axon indices
  for n in range(N): # Add N layers
    for ax in axs:
      for r in [0,1]:
        nebs = geo.branches[ax].neighborsAt(r)
        for n in nebs:
          ix = geo.branches.index(n)
          if upstream(geo, ax, ix) and ix not in axs:
            axs.append(ix)
  
  # Plotting
  if show:
    for ax in axs:
      c0, c1 = geo.branches[ax].coordAt(0), geo.branches[ax].coordAt(1)
      thiscol = color=np.random.random(3)
      plt.plot([c0[0], c1[0]], [c0[1], c1[1]], color=thiscol, lw=3, alpha=0.6)
      plt.text(c1[0], c1[1], ax, fontsize=10, color=thiscol)
    for s in geo.segments:
      c_ = s.coordAt(0)
      plt.plot(c_[0], c_[1], '.', color='gray', alpha=0.4)
    plt.title(geo.name)
    plt.show()
  return




def axonBranchpoint(geo, axInds):
  """
  Find the last branch points of the axon. Should be used interactively
  with idMoreAxons or idAxons. axInds are branch indices.
  """
  pDF = PathDistanceFinder(geo, geo.soma)
  axons = {ind: {'soma': None, 'seg': None} 
           for ind in axInds}
  
  # Get the axon distances and find the specific segment
  for ax in axInds:
    hdists = []
    for r in [0,1]:
      segs = [sg for sg in geo.segments if sg.name in geo.branches[ax].tags]
      for seg in segs:
        hdists.append(pDF.distanceTo(seg,r))
    mInd = hdists.index(min(hdists))
    axons[ax]['soma'] = hdists[mInd]
    if mInd > len(segs)-1:
      mInd = mInd - len(segs)
    axons[ax]['seg'] = segs[mInd]
  
  # Make the path distances grid
  for ax in axons.keys():
    pDF = PathDistanceFinder(geo, axons[ax]['seg'])
    axons[ax]['path_dists'] = {ax1: pDF.distanceTo(axons[ax1]['seg'], 0) 
                               for ax1 in axons.keys()}
    axons[ax]['euc_dists'] = {ax2: dist3(axons[ax]['seg'].nodeAt(0),
                                    axons[ax2]['seg'].nodeAt(0))
                              for ax2 in axons.keys()}
    axons[ax]['length'] = geo.branches[ax].length
  
  return axons




def axonBranchAnalysis(geo, axInds=None, axons=None):
  """
  Will plot the traces on a neuron, heatmap and fit
  """
  if axInds is None:
    if axons is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
  if axons is None:
    if axInds is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
    else:
      axons = axonBranchpoint(geo, axInds)
  
  collist = ['blue', 'red', 'forestgreen', 'goldenrod', 'purple', 'yellowgreen',
             'skyblue', 'tomato', 'darkgray']
  klist = list(axons.keys())
  # Show some stuff
  for s in geo.segments:
    c0 = s.coordAt(0)
    c1 = s.coordAt(1)
    plt.plot([c0[0], c1[0]], [c0[1], c1[1]], color='gray', lw=1, alpha=0.4)
  for ax in klist:
    pDF = PathDistanceFinder(geo, axons[ax]['seg'])
    coll = collist[klist.index(ax)]
    for ax1 in axons.keys():
      path = pDF.pathTo(axons[ax1]['seg'], 0)
      for s in path:
        c0 = s.coordAt(0)
        c1 = s.coordAt(1)
        plt.plot([c0[0], c1[0]], [c0[1], c1[1]], color=coll, 
                 lw=2.5, alpha=0.4)
    
    # Segment thing
    axcoord = axons[ax]['seg'].coordAt(0)
    plt.text(axcoord[0], axcoord[1], ax, color=coll, fontsize=15)
  plt.title(geo.name)
  
  # Now do the points and the fit
  plt.figure()
  xs, ys = [], []
  for ax in klist:
    coll = collist[klist.index(ax)]
    for u in axons[ax]['path_dists'].keys():
      plt.plot(axons[ax]['path_dists'][u], axons[ax]['euc_dists'][u], 'o',
               color=coll, markeredgecolor='none', alpha=0.6)
      xs.append(axons[ax]['path_dists'][u])
      ys.append(axons[ax]['euc_dists'][u])
  
  nxs, nys = [], []
  for i in range(len(xs)):
    if xs[i] not in nxs:
      nxs.append(xs[i])
      nys.append(ys[i])
  beta, alpha, r, p, yerr = stats.linregress(nxs, nys)
  plt.plot(nxs, [i*beta+alpha for i in nxs], '--', color='black',
           alpha=0.6, lw=2)
  plt.fill_between(nxs, [i*beta+alpha+yerr for i in nxs],
                   [i*beta+alpha-yerr for i in nxs], lw=0., facecolor='gray',
                   alpha=0.4)
  plt.text(min(plt.xlim()), max(plt.ylim())-10, 
           'beta : %.2f,  alpha : %.2f,  $\mathrm{r^2}$ : %.2f'
           %(beta, alpha, r**2))
  plt.xlabel('Path distance (um)')
  plt.ylabel('Euclidean distance (um)')
  plt.title(geo.name)
  plt.show()
  return




def branchesBetween(geo, axInds=None, axons=None):
  """
  Find the avg number of branches between two axon roots
  """
  if axInds is None:
    if axons is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
  if axons is None:
    if axInds is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
    else:
      axons = axonBranchpoint(geo, axInds)
  
  # Do the calcuations
  segs = {}
  for ax in axons.keys():
    pDF = PathDistanceFinder(geo, axons[ax]['seg'])
    for ax2 in axons.keys():
      if ax != ax2:
        segs['%i-%i' %(ax, ax2)] = [i.name for i in pDF.pathTo(axons[ax2]['seg'])]
  
  # Now figure out how many branches compose these segments
  blist = {}
  for b in geo.branches:
    for t in b.tags:
      blist[t] = geo.branches.index(b)
  branches = []
  for segkey in segs.keys():
    temp_branches = [blist[i] for i in segs[segkey]]
    branches.append(len(list(set(temp_branches))))
    print('Branches between %s: %i' %(segkey, branches[-1]))
  
  print('Average inter-branch number: %.2f +/- %.2f' %(np.mean(branches), np.std(branches)))
  return np.mean(branches)

#


def branchScramble(geo, dists, axInds=None, axons=None, N=1000, show=True):
  """
  Gets the inter-branch intervals and then makes many combinations of
  branches (to get their lengths) to get a distribution on which the
  true inter-branch distances (dists) can be plotted. (N trials.)
  """
  if axInds is None:
    if axons is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
  if axons is None:
    if axInds is None:
      print('Need either axInds or axons (dict from axonBranchpoints)')
      return None
    else:
      axons = axonBranchpoint(geo, axInds)
  
  # Get the avg inter-branch number
  interbr = branchesBetween(geo, axons=axons)
  # scrams = []
  trials = [np.random.randint(0, len(geo.branches), int(interbr)+1)
            for u in range(N)]
  scrams = [sum([geo.branches[ix].length for ix in tri[:-1]]) 
            + geo.branches[tri[-1]].length*(interbr%1)
            for tri in trials]
 # scrams = [i + geo.branches[tri[-1]].length*(interbr%1) for tri in trials]
  
  # Calculate the p-values
  print(scrams[:10], np.mean(scrams))
  pvals = [float(len([s for s in scrams if s < di]))/
           float(len(scrams)) for di in dists]
  print('Chance of occurring randomly: %.5f +/- %.5f '
        %(np.mean(pvals), np.std(pvals)))
  
  # Plotting
  if show:
    plt.figure()
    plt.hist(scrams, bins=50, color='blue', edgecolor='white',
             alpha=0.7)
    plt.axvline(np.mean(pvals), 0, plt.ylim()[-1], linestyle='--', color='red', lw=1.5)
    plt.axvspan(np.mean(pvals)-np.std(pvals), np.mean(pvals)+np.std(pvals),
                color='red', alpha=0.4)
    plt.text(np.mean(pvals), plt.ylim()[-1]-10, '$\mathrm{p =}$ %.5f' %np.mean(pvals))
    plt.show()
  
  return



########################################################################
# Other plotting functions


def pathStack(geo, diamfile, voxel=[.732,.732,.488], colors=None):
  """
  This writes the paths to an image stack so the paths can be visualized
  in 3D.
  """
  targ_dict = getDiameters(geo, diamfile, voxel=voxel, report=False, show=False)
  pDF = PathDistanceFinder(geo, geo.soma)
  paths = {t:[] for t in targ_dict.keys()}
  
  # Each path should be a collection of x,y,z points, as ints, 
  #   to write to an img stack
  for t in targ_dict.keys():
    t_seg = [seg for seg in geo.segments 
             if seg.filamentIndex == targ_dict[t]['closestFil']][0]
    try:
      segs = pDF.pathTo(t_seg)
    except:
      print('No segments match filament index %i' %targ_dict[t]['closestFil'])
    # Add the nodes as [x,y,z] tuples
    for s in segs:
      for n in s.nodes:
        if [n.x, n.y, n.z] not in paths[t]:
          paths[t].append([n.x, n.y, n.z])
    
  # Replace the float paths with ints
  for t in paths.keys():
    for c in range(len(paths[t])):
      paths[t][c] = [int(paths[t][c][i]/voxel[i]) 
                         for i in range(3)]
  return paths



def pathsIn3D(geo, diamfile=None, paths=None, voxel=[.732,.732,.488],
              colors=None):
  """
  """
  if diamfile is None and paths is None:
    print('Need either diamfile or paths!'); return
  if paths is None:
    paths = pathStack(geo, diamfile, voxel=voxel, colors=colors)
  
  # Interpolate pts
  newpaths = {}
  for p in paths.keys():
    print('Interpolating %i / %i paths ...' 
          %(list(paths.keys()).index(p)+1, len(paths.keys())))
    temp_ = []
    for u in range(len(paths[p])-1):
      internode = dist3(paths[p][u], paths[p][u+1])
      if internode > 1.5: # Add nodes in between
        axes = [np.linspace(paths[p][u][ix], paths[p][u+1][ix], int(internode)) 
                for ix in range(3)]
        for a in range(len(axes)):
          temp_.append([int(axes[0][a]), int(axes[1][a]), int(axes[2][a])])
      else:
        temp_.append(paths[p][u])
        temp_.append(paths[p][u+1])
    newpaths[p] = temp_
  
  dims = [ max([max([pt[i] for pt in newpaths[p]]) for p in newpaths.keys()])
           for i in range(3)]
  if colors is None:
    colors = [np.random.random(3) for i in newpaths.keys()]
  
  # Create the image array and add the path points to it
  imarr = [ [ [ [0.,0.,0.] for j in range(dims[1]+1)] # y-range
             for i in range(dims[0]+1)] # x-range
           for k in range(dims[2]+1)] # z is first for slice indexing
  for p in newpaths.keys():
    print('Writing path %i / %i to image ...' 
          %(list(paths.keys()).index(p)+1, len(paths.keys())))
    for pt in newpaths[p]: # For each point
      #print(pt)
      imarr[pt[2]] [pt[0]] [pt[1]] = colors[list(newpaths.keys()).index(p)]
  
  return imarr
  


def arrStack(imarr, imdir):
  """
  Save an np.arr(dtype=np.uint8) as a tiff stack.
  """
  return



def convertToRGB(imdir):
  # Convert a 1-int B&W image to an RGB 3x1 uint8 array
  # Load up the image stack as black and white
  print('Loading the image from %s ... ' %imdir)
  imArr = []
  imlist = os.listdir(imdir)
  imlist = [f for f in imlist if f.split('.')[-1] == 'tif']
  imlist = [imdir + f for f in imlist]
  for im in imlist:
    if imlist.index(im) % 20 == 0:
      print('   %i/%i images processed' %(imlist.index(im), len(imlist)))
    temp_im = np.array(Image.open(im))
    add_im = [ [ [] for i in range(temp_im.shape[0])] 
              for j in range(temp_im.shape[1])]
    for i in range(temp_im.shape[0]):
      for j in range(temp_im.shape[1]):
        add_im[i][j] = [np.uint8(u) for u in [temp_im[i][j] for g in range(3)]]
    imArr.append(np.array(add_im))
  return imArr



def addStacks(geo, diamfile, imdir, tostack, 
              colors=None, voxel=[.732,.732,.488], width=6):
  """
  Adds two stacks together. A path stack is created from pathStack
  and is added to an image stack (from imdir). 
  The paths are padded with *width* voxels on all sides to make them more visible.
  Make sure tostack is an existing directory! imdir can also be the array itself.
  """
  try:
    p_ = os.listdir(tostack)
  except:
    print('Directory %s does not exist!' %tostack)
  
  paths = pathStack(geo, diamfile, voxel, colors)
  # Set the colors
  if colors is None:
    colors = [np.random.random(3) for c in range(len(paths.keys()))]
  padwidth = [] # Create the padding offsets (+/- width voxels)
  for p in range(1,width+1):
    padwidth.append(-p)
    padwidth.append(p)
  print('Padding:'); print(padwidth)
  # Get the image array into an RGB format
  if type(imdir) is str:
    imArr = convertToRGB(imdir) 
  else:
    print('Assuming imdir is the array...')
  imArr = np.array(imdir)
  
  print(imArr.shape)
  # Add the paths to this array
  print('Padding the directory with width=%i' %width)
  skip = 0
  for t in paths.keys(): # For each path
    print('Path %i of %i ....' %(list(paths.keys()).index(t), len(paths.keys())))
    for c in range(len(paths[t])): # For each node/coordinate
      for i in padwidth: # i is X-dim (actually y) -- "rows"
        for j in padwidth: # j is Y-dim (actually x) -- "cols"
          for k in padwidth: # k in Z-dim
            coords = []
            for test in [[2,k], [0,i], [1,j]]:
              try:
                coords.append(paths[t][c][test[0]]+test[1])
              except:
                coords.append(0)
            z, x, y = coords
            if x == y == z == 0:
              skip = skip + 1
            else:
              imArr[z][x][y] = colors[list(paths.keys()).index(t)]
  
  # Make sure all images have RGB values, not just B&W
  print('Skipped %i values.\nChecking all values are RGB-compatible...' %skip)
  return imArr
  
  for z in range(len(imArr)):
    
    for x in range(len(imArr[z])):
      for y in range(len(imArr[z][x])):
        try:
          L = len(imArr[z][x][y])
        except:
          L = 1
        if L == 1:
          imArr[z][x][y] = [np.uint8(i) for i in [imArr[z][x][y] for qq in range(3)]]
        elif L == 3:
          imArr[z][x][y] = [np.uint8(i) for i in imArr[z][x][y]]
        else:
          print('Shit! Only found %i values for %i, %i, %i (z,x,y)'
                %(len(imArr[z][x][y]), z, x, y))
          imArr[z][x][y] = [np.uint8(i) for i in [0,0,0]] # Paint it black
  
  # Now can write the whole image stack
  cnt = 0
  for zslice in imArr:
    zim = Image.fromarray(zslice)
    zim.save(tostack+geo.name+'_%i.tif' %cnt)
    cnt = cnt + 1
  print('Saved %i images into stack %s' %(cnt, tostack))
  
  return



########################################################################
# Tools

def simple_csv(ddict, outfile, idx=False):
  """ Write a simple dictionary to a csv file. """
  nrows = len(ddict.keys())
  
  # Turn 'em into a string and write 'em
  with open(outfile, 'w') as fOut:
    try:
      sd = list(ddict[list(ddict.keys())[0]].keys())
    except:
      sd = list(ddict.keys())
    fOut.write(','.join([str(u) for u in sd]))
    fOut.write('\n')
    
    for k in ddict.keys():
      if idx:
        ddict[k]['idx'] = k
      try:
        fOut.write(','.join([str(j) for j in ddict[k].values()]))
      except:
        fOut.write(','.join([str(j) for j in ddict[k]]))
      fOut.write('\n')
  print('%s written' %outfile)
  return

#


def row_entries(listofrows, rowheaders, outfile):
  """Easily writes a list of data with row titles to a data file, 
  i.e.: row_entries(data, [g.name for g in geofiles], outfile)
  """
  if len(listofrows) != len(rowheaders):
    print('Need same number of rowheaders as rows of data!')
    return None
  with open(outfile,'w') as fOut:
    for row in range(len(listofrows)):
      temp = [str(rowheaders[row])]
      for item in listofrows[row]:
        temp.append(str(item))
      fOut.write(','.join(temp))
      fOut.write('\n')
  print('%s written' %outfile)
  return










########################################################################

if __name__ == "__main__":
  print('Module is used interactively.')






































