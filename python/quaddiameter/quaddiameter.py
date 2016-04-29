# quaddiameter.py - Using the quadratic tapering constants dervied
#                   by H Cuntz et al (TREES Toolbox), fit a quadratic
#                   taper to a given skeleton.
# usage: python quaddiameter.py hocFile fittedFile
# Here, the radii from hocFile are ignored; fittedFile is a hoc file
# with the same nodes but with different radii (quad-fit).


import os, sys
sys.path.append('/home/alex/code/morphology/python/functions/')
import numpy as np
import matplotlib.pyplot as plt
from neuron_readExportedGeometry import *

  

###########################################################
class Quadfit():
  """
  Fit quadratic taper to the given neuron. primNeurDiam by default
  is the mean of initial diameters.
  """
  def __init__(self, hocfile, primNeurDiam=15.4, newHoc=True):
    self.hocfile = hocfile
    self.newhocfile = hocfile.split('.')[0]+'_quaddiam.hoc'
    self.pts, self.paths = [], []
    self.geo = demoReadsilent(self.hocfile)
    self.P, self.ldend = self.load_params()
    self.primNeurDiam = primNeurDiam
    self.properties = {} # Properties
    self.usedSegs, self.missed = [], []
    self.tipInds, self.tipLocs = self.geo.getTipIndices()
    self.pDF = PathDistanceFinder(self.geo, self.geo.soma) # Same pDF is used
    self.filamentInds = {s.filamentIndex: None for s in self.geo.segments}
    
    # Populate the seg dict by filamentIndex
    self.segDict()
    self.axonDiameter()
    self.assignAll()
    self.zeroRads()
    self.getProps()
    if newHoc:
      self.as_hoc()
    return
  
  
  
  def load_params(self):
    """
    Load the parameters for the fits. These are derived in Cuntz... Segev 
    (2007). For each segment of normalized length ldend[i], P[i] are the
    parameters for  y = P[i,0]x^2 + P[i,1]x + P[i,2] = 
    """
    P = []
    with open('P.txt','r') as fIn:
      for line in fIn:
        if line:
          splitLine = line.split(None)
          P.append([float(i) for i in splitLine])
    
    # All ldend values are tab or space-separated on one line
    ldend = []
    with open('ldend.txt', 'r') as fIn:
      for line in fIn:
        if line:
          splitLine = line.split(None)
          for s in splitLine:
            ldend.append(float(s))
    return P, ldend
  
  
  
  def solve_path(self, Pars, path, ldend, masterrad):
    """
    Fit diamaters by path lengths. Given current params (len=3), the
    path (list of seg instances) and the length to fit to. Masterrad
    scales the radii 
    """
    plength = sum([s.length for s in path])
    nnodes = sum([len(s.nodes) for s in path])
    xs = np.linspace(0,ldend, nnodes)
    # Rads has same number of nodes as path
    rads = [Pars[0]*i**2 + Pars[1]*i + Pars[2] for i in xs]
    rads = [r/max(rads)*masterrad for r in rads]
    cnt = 0
    
    # Assign the radii for the segs in the path
    for s in path:
      avg_over_n = len(s.nodes)
      self.filamentQuads[s.filamentIndex]['radius'] = mean(rads[cnt:cnt+avg_over_n])
      self.usedSegs.append(s)
      cnt = cnt + avg_over_n
    return self
  
  
  
  def uniquePath(self, potpath):
    """
    Given a potential path, return only the unused segments.
    """
    for seg in potpath:
      if seg in self.usedSegs:
        potpath.pop(potpath.index(seg))
    self.potentialPath = potpath
    return self
    
    
    
  def segDict(self):
    """
    Populate the seg dict by filamentIndex
    """
    # Get the paths
    for seg in self.geo.segments:
      if seg.filamentIndex in self.filamentInds.keys():
        self.filamentInds[seg.filamentIndex] = seg
    for t in range(len(self.tipInds)):
      self.paths.append(self.pDF.pathTo(self.filamentInds[self.tipInds[t]], self.tipLocs[t]))
      
    # Map path lengths
    self.path_lengths = [sum([p.length for p in path]) for path in self.paths] 
    # Initialize filamentQuads
    self.filamentQuads = {s.filamentIndex: {'length': s.length} for s in self.geo.segments}
    # For longest path, make the initial segment diam match the given diam
    # Then taper the others as needed
    self.long_path = max(self.path_lengths)
    
    # Set the longest path first
    idx = self.path_lengths.index(self.long_path)
    self.solve_path(self.P[-1], self.paths[idx], self.ldend[-1], self.primNeurDiam)
    self.paths.pop(idx) # Remove already used path
    
    # Now do it for all other paths
    gotone = False
    for p in self.paths:
      prim_cnt = 1 # Keep track of these to make sure they equal the number of paths
      self.uniquePath(p)
      for seg in self.potentialPath: # Get the radius of the 'root' of the path
        for neb in seg.neighbors:
          if neb in self.usedSegs: # This is the root
            newPrim = self.filamentQuads[neb.filamentIndex]['radius']
            gotone = True # Did it get at least one new root?
      
      if gotone:
        prim_cnt = prim_cnt + 1 # Increment the new prim count
            
      # Find the lengths and params
      p_idx = int(self.ldend[int(sum([s.length for s in self.potentialPath])
                                 /self.long_path)])
      self.solve_path(self.P[p_idx], self.potentialPath, self.ldend[p_idx], newPrim)
    # Should now have radii for all the segments in all the paths
    print('Found %i/%i paths' %(prim_cnt, len(self.paths)+1))
    return self
  
  
  
  def axonDiameter(self):
    """
    Set the diameters for the axons segments.
    """
    unused = [seg for seg in self.geo.segments if seg not in self.usedSegs]
    farthest = max([self.pDF.distanceTo(seg, 0.5) for seg in unused])
    axCnt = 0
    # Scale the axon by length of the longest path 
    #  (axon cannot get thinner than the thinnest tip)
    if farthest > self.long_path:
      get = 1
    else: 
      get = farthest/self.long_path
    
    # For each path, get the Param values as a fraction of farthest
    for seg in unused:
      idx = int(self.pDF.distanceTo(seg, 0.5)/farthest * get * len(self.P))
      # Set the new diameter
      x = self.pDF.distanceTo(seg, 0.5)/farthest*self.ldend[int(get)]
      diam = self.P[int(get)][0]*x**2 + \
             self.P[int(get)][1]*x + \
             self.P[int(get)][2]
      # Scale the radius by the soma rad
      diam = diam/self.P[int(get)][2] * self.primNeurDiam
    
      # If this already has a radius, don't replace it
      radthing = None
      try:
        radthing = self.filamentQuads[seg.filamentIndex]['radius']
      except:
        pass
      if radthing is None:
        self.filamentQuads[seg.filamentIndex]['radius'] = diam
        self.usedSegs.append(seg)
        axCnt = axCnt + 1
      
    print('Assigned %i radii from axon segments' %axCnt)
    return self
  
  
  
  def assignAll(self):
    """
    If there were any unassigned radii, assign them now.
    """
    for seg in self.geo.segments:
      filInd, filRad = None, None
      try:
        filInd = self.filamentQuads[seg.filamentIndex]
        filRad = self.filamentQuads[seg.filamentIndex]['radius']
      except:
        pass
      
      # If this filament doesn't have an entry, make it
      filCnt, radCnt = 0, 0
      if filInd is None:
        self.filamentQuads[seg.filamentIndex]['length'] = seg.length
        filCnt = filCnt + 1
      
      # Or just add the radius by averaging neighboring radii
      if filRad is None: # Will always be None if filInd is none
        nebrads = []
        for neb in seg.neighbors:
          try:
            nebrads.append(self.filamentQuads[neb.filamentIndex]['radius'])
          except:
            pass
        self.filamentQuads[seg.filamentIndex]['length']['radius'] = np.mean(nebrads)
        radCnt = radCnt + 1
    
    print('Created %i new filaments; assigned %i new radii' 
          %(filCnt, radCnt))
    return self
  
  
  
  def zeroRads(self):
    """
    Make sure there are no zero radii.
    """
    # Get the minimum radius that isn't 0.
    minrad = np.inf
    for segkey in self.filamentQuads.keys():
      currad = self.filamentQuads[segkey]['radius']
      if currad < np.inf and currad > 0.0001:
        minrad = currad
    
    # Assign all zero rads this new small, non-zero value
    for segkey in self.filamentQuads.keys():
      if self.filamentQuads[segkey]['radius'] <= 0.0001:
        self.filamentQuads[segkey]['radius'] = minrad
    
    return self
    
  
  
  def getProps(self):
    """
    Get some rudimentary properties of the quad object.
    """
    minrad, maxrad = np.inf, 0.
    tiprads, somarad = [], None
    
    for segkey in self.filamentQuads.keys():
      currad = self.filamentQuads[segkey]['radius']
      
      # Get soma rad
      if segkey == self.geo.soma.filamentIndex:
        somrad = currad
      
      # Get extreme rad stuff
      if currad > maxrad:
        maxrad = currad
      if currad < minrad:
        minrad = currad
      
      # It it's a tip, add it to the tips list
      if segkey in self.tipInds:
        tiprads.append(self.filamentQuads[segkey]['radius'])
    
    # Put all of these into a dict
    self.properties['min radius'] = minrad
    self.properties['max radius'] = maxrad
    self.properties['tip radii'] = tiprads
    self.properties['soma radius'] = somarad
    
    return self
  
  
  
  def seg_points(self, fil_name):
    """
    Return the x,y,z,quadrad for points in the given segment.
    """
    pts, fil_ind = [], 'no matching segment found!'
    winner_seg = None
    for seg in self.geo.segments:
      if seg.name == fil_name:
        winner_seg = seg
        try:
          fil_ind = seg.filamentIndex
        
          for n in seg.nodes: # For each 'node'
            pts.append([n.x, n.y, n.z,
                self.filamentQuads[fil_ind]['radius']])
        except:
          print(' -->  missed: %s' %fil_name)
          self.missed.append(fil_name)
          
    try:
      print('Found %i/%i points for segment %s ' %(len(self.pts),
                                                   len(winner_seg.nodes),
                                                   winner_seg.name))
      self.pts = pts
    except:
      print('No matching segment for %s' %fil_name)
      self.pts = None
    return self
    
  
  
  def as_hoc(self):
    """
    This saves the quad object as a hoc file with the new diameters.
    """
    print('Writing %s ... ' %self.newhocfile)
    with open(self.hocfile, 'r') as fIn:
      with open(self.newhocfile, 'w') as fOut:
        skip = 0
        l = 0 # Track the lines
        for line in fIn:
          if skip > 0: # Skip lines to accomodate for used pt3dadd lines
            skip = skip - 1
            l = l + 1
            pass
          
          # This should only be entered if skip == 0
          else:
            if line:
              if len(line.split('_')) > 0:
                if line.split('_')[0] == 'filament': # Found a seg, write it
                  # get the points for this segment and write them
                  self.seg_points(line.split(None)[0])
                  fOut.write(line)
                  fOut.write('  pt3dclear()\n')
                  for p in self.pts:
                    fOut.write('  pt3dadd(%.4f,%.4f,%.4f,%.4f)\n'
                               %(p[0], p[1], p[2], p[3]))
                  fOut.write('}\n')
                  skip = len(self.pts) + 2
                
                else:
                  fOut.write(line)
                  skip = 0
              else:
                fOut.write(line)
                skip = 0
            else:
              fOut.write(line)
              skip = 0
        
        # Now should have all lines written
    # Done with writing new file
    print('File %s written' %self.newhocfile)
    return
    

  
def compareQuad_Geo(quad):
  """
  quad object should have geo as an attribute.
  """
  # Compile the radii to compare
  compar = [] # Items: [ [filInd, quadRad, geoRad], ... ]
  klist = list(quad.filamentQuads.keys())
  
  for k in klist:
    try:
      temp = [k, quad.filamentQuads[k]['radius']]
      for seg in quad.geo.segments:
        if seg.filamentIndex == k:
          temp.append(seg.avgRadius)
      compar.append(temp)
    except:
      pass
      
  # Compare
  dif = [c[1]/c[2] for c in compar]
  print('%i segments differ by %.3f +/- %.3f\n   range: %.3f - %.3f'
         %(len(dif), np.mean(dif), np.std(dif), min(dif), max(dif)))
  return compar
  



############################################################

# Analysis

def compare_tips(quad, show=False):
  """
  Compare the radii at the tips between the quad and geo objects.
  """
  # Structure of data: [filInd, quad, geo]
  tiprads = {t: [quad.filamentQuads[t]['radius']] for t in quad.tipInds}
  for seg in quad.geo.segments:
    if seg.filamentIndex in tiprads.keys():
      tiprads[seg.filamentIndex].append(seg.nodes[int(len(seg.nodes) *
                                                      quad.tipLocs[quad.tipInds.index(seg.filamentIndex)])-1].r1)
  # Plotting
  if show:
    fig = plt.figure()
    ax1 = fig.add_subplot(211) # Paired-line plot
    ax2 = fig.add_subplot(212) # A few random cross-sections?
    
    # Paired-line plot
    for segkey in tiprads.keys():
      #print(tiprads[segkey])
      ax1.plot([0,1], tiprads[segkey], linewidth=1, color='gray',
               alpha=0.2)
      ax1.scatter(0, tiprads[segkey][0], marker='o', color='lightcoral', 
               edgecolor='lightcoral', alpha=0.2)
      ax1.scatter(1, tiprads[segkey][1], marker='o', color='palegreen',
               edgecolor='palegreen', alpha=0.2)
    ax1.set_xlim((-0.2, 1.25))
    # ax1.set_ylim((0,quad.properties['max radius']*1.1))
    qmean = np.mean([tiprads[k][0] for k in tiprads.keys()])
    gmean = np.mean([tiprads[k][1] for k in tiprads.keys()])
    ax1.plot([0,1], [qmean, gmean], linewidth=3, color='k', alpha=0.8)
    
    # Also make a line at the mean
    ax1.plot([-.1,.1], [qmean, qmean], linewidth=3, color='lightcoral', alpha=0.6)
    ax1.text(-.2, qmean-0.1, '%.2f' %qmean, fontsize=14)
    ax1.plot([.9,1.1], [gmean, gmean], linewidth=3, color='palegreen', alpha=0.6)
    ax1.text(1.1, gmean-0.1, '%.2f' %gmean,  fontsize=14)
    
    # More plot specifics
    ax1.text(-.1,-.15, 'Quad radius', fontsize=15)
    ax1.text(.9,-.15, 'Actual radius', fontsize=15)
    for pos in ['left', 'top', 'right', 'bottom']:
      ax1.spines[pos].set_visible(False) 
    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])
    
    # Random tip cross-sections
    r_idx = [list(tiprads.keys())[int(u)] for u in 
             np.random.random([4])*len(tiprads.keys())] # Just 4 samples
    print('Random tips are:', r_idx)
    # Get the relevant radii
    quadrads, georads = [tiprads[i][0] for i in r_idx], [tiprads[i][1] for i in r_idx]
    maxrad = max([max(quadrads), max(georads)])
    
    # Plot the circles
    centers = [[maxrad, 3*maxrad], [3*maxrad, 3*maxrad], 
               [maxrad, maxrad], [3*maxrad, maxrad]]
    cirs = []
    for c in range(4):
      cirs.append(plt.Circle((centers[c][0], centers[c][1]),
                             radius=georads[c], fc='palegreen', 
                             edgecolor='palegreen', alpha=0.5))
      cirs.append(plt.Circle((centers[c][0], centers[c][1]),
                             radius=quadrads[c], fc='lightcoral',
                             edgecolor='lightcoral', alpha=0.9))

    for t in range(4):
      ax2.text(centers[t][0]-.2, centers[t][1], 'Tip %i' %r_idx[t])
    for cir in cirs:
      plt.gca().add_artist(cir)
    
    # Clean it up
    ax2.set_xlim((0,4*maxrad))
    ax2.set_ylim((0,4*maxrad))
    for pos in ['left', 'top', 'right', 'bottom']:
      ax2.spines[pos].set_visible(False) 
    ax2.get_xaxis().set_ticks([])
    ax2.get_yaxis().set_ticks([])
    ax2.plot([0,1], [0,0], 'k', linewidth=3)
    ax2.text(0,0, '1 um')
    plt.show()
    
  return tiprads

  

def branch_ratios(quad):
  """
  Calculate branch ratios and compare between quad and geo object.
  """
  geoRads = {seg.filamentIndex: {'radius': seg.avgRadius,
                                 'daughters': []} for seg in quad.geo.segments}
  quadRads = {k: {'radius': quad.filamentQuads[k]['radius'],
                  'daughters': []} for k in list(quad.filamentQuads.keys())}
  
  # Now go through paths and find downstream neighbors
  # Parent segment cannot go through its daughters on its way to soma
  for path in quad.paths:
    for seg in path:
      if len(seg.neighbors) > 2:
        nebs = seg.neighbors
        
        for n in nebs:
          try:
            if n in path[::-path.index(seg)]: # If it's in path back to soma
              go = False
            else:
              go = True
          except: # This is the 0th seg in the path, must be parent
            if path.index(seg) == 0:
              go = True
            else:
              go = False
        
        # Else, get the radii
        if go:
          for n in nebs:
            geoRads[seg.filamentIndex]['daughters'].append(n.avgRadius)
          go = False # Reset this thing
  
  # Now have all the radii 
  return quadRads, geoRads



def get_ratios(quad, quadRads=None, geoRads=None):
  """
  Show the branch ratios; can calculate them fresh if no rads given.
  """
  def getPair(quadRads, geoRads, k, num):
    return [quadRads[k]['radius']/quadRads[k]['daughters'][num],
            geoRads[k]['radius']/geoRads[k]['daughters'][num]]
  
  if quadRads is None and geoRads is None:
    quadRads, geoRads = branch_ratios(quad) # Get the radii
  
  # Clean dictionaries first
  for k in quadRads.keys():
    if len(quadRads[k]['daughters']) == 0:
      quadRads.pop(k)
  for k in geoRads.keys():
    if len(geoRads[k]['daughters']) == 0:
      geoRads.pop(k)
  
  # Paired-line data
  paired, quadSingle, geoSingle = [], [], []
  for k in quadRads.keys():
    if k in geoRads.keys(): # Found a match
      
      # If same number of daughters, they are paired
      if range(len(quadRads[k]['daughters'])) == range(len(geoRads[k]['daughters'])):
        for d in range(len(quadRads[k]['daughters'])):
          paired.append(getPair(quadRads, geoRads, k, d))
      
      # Else, log them separately
      else:
        for d in quadRads[k]['daughters']:
          quadSingle.append(quadRads[k]['radius']/d)
        for d in geoRads[k]['daughters']:
          geoSingle.append(geoRads[k]['radius']/d)
  
  # Got all the branch ratios
  return paired, quadSingle, geoSingle



def load_params(addpath='/home/alex/code/morphology/python/build-morphology/'):
  """
  Load the parameters for the fits. These are derived in Cuntz... Segev 
  (2007). For each segment of normalized length ldend[i], P[i] are the
  parameters for  y = P[i,0]x^2 + P[i,1]x + P[i,2] = 
  """
  P = []
  with open(addpath+'P.txt','r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        P.append([float(i) for i in splitLine])
  
  # All ldend values are tab or space-separated on one line
  ldend = []
  with open(addpath+'ldend.txt', 'r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        for s in splitLine:
          ldend.append(float(s))
  return P, ldend



def get_quad_tips(geo=None, primNeurDiam=15.4, minTip=0.0001, paths=None):
  """
  For a simpler analysis, this just fits the paths to a quad taper
  and returns what the tip should be (cannot be less than minTip).
  """
  def solve(params, x):
    return params[0]*x**2 + params[1]*x + params[2]
  
  # Get paths first
  if paths is None:
    if geo is not None:
      pDF = PathDistanceFinder(geo, geo.soma)
      tipInds, tipLocs = geo.getTipIndices()
      # Get the segments
      tipSegs = []#[s for s in geo.segments if s.filamentIndex==ind] for ind in tipInds]
      for i in tipInds:
        tipSegs.append([s for s in geo.segments if s.filamentIndex==i][0])
      paths = [pDF.distanceTo(tipSegs[u], tipLocs[u]) for u in range(len(tipInds))]
  
  # Fit longest path to longest taper and find the other fits
  P, ldend = load_params() # Load parameters
  norm_pathlengths = [(p/max(paths))*max(ldend) for p in paths]
  inds = []
  for p in norm_pathlengths:
    for l in range(len(ldend)):
      if ldend[l] > p: # If it just became greater than the path, use it
        inds.append(l)
        break
  
  # Get the 'fitted' tip diameters
  tip_diams = [(solve(P[u], ldend[u])/solve(P[u],0))*primNeurDiam
               for u in inds]
  return tip_diams
  













  
#############################################################
if __name__ == "__main__":
  args = sys.argv
  if len(args) > 2:
    hocfile = args[1]
    Quadfit(hocfile)
























