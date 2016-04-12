#!/usr/bin/python



_usageStr=\
"""usage: neuron_readExportedGeometry.py geoFile
  get dictionary object describing neuron model geometry info by reading file
"""



import os, sys, re, math
from NeuronGeometry import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
import math
import pickle



class HocGeometry(Geometry):
  def __init__(self, _fileName=None):
    Geometry.__init__(self)
    self._openFilament = None
    self._connections = []
    self.connections = []
    self._filamentNames = []
    self._filaments = {}
    self._filamentNameType = None
    self._warnRepeatFilaments = True
    
    if _fileName is not None:
      self.setFileName(_fileName)
      self.readGeometry()
      
  
  def readGeometry(self):
    """
    get dictionary object describing neuron model geometry info by reading file
    """
     
    lineNum = 0
    with open(self.fileName, 'r') as fIn:
      # read the geometry file
      try:
        for line in fIn:
          # loop through each line in the file
          
          # inc the line number
          lineNum = lineNum + 1
          # parse the line in geometry file, adding info to geometryInfo
          self._parseHocGeometryLine(line)
    
      except IOError as err:
        sys.tracebacklimit = 0
        raise IOError('Error reading %s line %d: %s' % \
                      (self.fileName, lineNum, err.message))
    
    if self._openFilament:
      raise IOError('Error reading %s, filament %s open at end of file' %
                    (self.fileName, self._openFilament))
    
    # connect filaments and remove filaments and _connections, leaving segments
    # and nodes
    self._connectFilaments()
    
    
    # make compartments from hemispheres remaining at the end of unconnected
    # segments
    #self._addOneNodeCompartments()
    

  def getSomaIndex(self):
    """
    return (filamentIndex, position)
      filamentIndex indexes the .hoc file filament that contains the Soma
      position is a float between 0 and 1 that points to the soma centroid on
        the filament
    """
    # get the Soma
    soma = self.soma
    filamentIndex = self.getFilamentIndex(soma)
    
    # get the centroid of the Soma
    centroid = soma.centroidPosition(mandateTag='Soma')

    return (filamentIndex, centroid)

  
  def getTipIndices(self):
    """
    return (filamentInds, positions)
      filamentInds is a list of indices to .hoc file filaments that contain
        terminal segments
      positions is a list of floats (0 or 1) that point to the terminal end
        of each terminal segment
    NOTE: This will NOT contain Axon or Soma even if they are terminal segments
    """
    self.checkConnectivity(removeDisconnected=True)
    
    soma = self.soma
    axons = self.findAxons()
    
    def _terminalEnd(seg):
      n0, n1 = False, False
      for loc, nLoc, node in seg.neighborLocations:
        if loc == 0.0:
          if n1:
            return None
          else:
            n0 = True
        elif loc == 1.0:
          if n0:
            return None
          else:
            n1 = True
      if n0:
        return 1.0
      else:
        return 0.0
    
    ends = ((self.getFilamentIndex(s), _terminalEnd(s)) for s in self.segments
            if (s != soma and s not in axons))
    try:
      filamentInds, positions = zip( *((f, e) for f, e in ends
                                     if e is not None))
    except ValueError:
      raise ValueError('No tip indices found!?!')
    
    return filamentInds, positions


  def getSegLengths(self):
    """
    return (segLength)
    """
    cnt = 0
    for segs in len(range(self.segments)):
      cnt += 1
    
    return cnt

  
  def getTips(self):
    """
    return (filamentInds, positions)
      filamentInds is a list of indices to .hoc file filaments that contain
        terminal segments
      positions is a list of floats (0 or 1) that point to the terminal end
        of each terminal segment
    NOTE: This will NOT contain Axon or Soma even if they are terminal segments
    """
    self.checkConnectivity(removeDisconnected=True)
    
    soma = self.soma
    axons = self.findAxons()
    
    def _terminalEnd(seg):
      n0, n1 = False, False
      for loc, nLoc, node in seg.neighborLocations:
        if loc == 0.0:
          if n1:
            return None
          else:
            n0 = True
        elif loc == 1.0:
          if n0:
            return None
          else:
            n1 = True
      if n0:
        return 1.0
      else:
        return 0.0
    
    ends = ((s, _terminalEnd(s)) for s in self.segments
            if (s != soma and s not in axons))
    try:
      terminalSegs, positions = zip( *((f, e) for f, e in ends
                                     if e is not None))
    except ValueError:
      raise ValueError('No tip indices found!?!')
    
    return terminalSegs, positions
  
  
  def getAxonIndices(self):
    """
    return (filamentInds, positions)
      filamentInds is a list of indices to .hoc file filaments that contain
        terminating branches
      positions is a list of floats (0 or 1) that point to the terminal end
        of each terminating branch
    """
    self.checkConnectivity(removeDisconnected=True)
    axons = self.findAxons()
    if not axons:
      return [], []
    
    def _terminalEnd(seg):
      n0, n1 = False, False
      for loc, nLoc, node in seg.neighborLocations:
        if loc == 0.0:
          assert not n1, 'Axon is an isolated segment'
          n0 = True
        elif loc == 1.0:
          assert not n0, 'Axon is an isolated segment'
          n1 = True
      if n0:
        return 1.0
      else:
        return 0.0
    
    ends = ((self.getFilamentIndex(s), _terminalEnd(s)) for s in axons)
    filamentInds, positions = zip( *((f, e) for f, e in ends) )
    
    return filamentInds, positions


  def _parseHocGeometryLine(self, line):
    """
    Read a line from hoc file specifying geometry, and update geometryInfo
    appropriately.
    openFilament = name of filament if in a declaration block, otherwise = None
    """
    splitLine = line.split(None)
    if not splitLine:
      return
    
    if self._openFilament:
      self._parseDefineFilament(line)
    elif splitLine[0] == 'connect':
      self._addConnection(splitLine)
    elif splitLine[0] == 'create':
      self._createFilaments(splitLine)
    elif splitLine[0] == 'neuron_name':
      self.name = splitLine[-1]
    elif splitLine[0].lower() == "range":
      if len(splitLine) < 7:
        raise IOError(\
          'range should be of form "range minX maxX minY maxY minZ maxZ"')
      self.minRange = tuple([float(x) for x in splitLine[1:6:2]])
      self.maxRange = tuple([float(x) for x in splitLine[2:7:2]])
    elif splitLine[0] in self._filamentNames:
      self._openFilament = splitLine[0]
    elif splitLine[0]+'[0]' in self._filamentNames:
      self._openFilament = splitLine[0]+'[0]'
  
  def _parseDefineFilament(self, line):
    """
    Parse a line in a filament declaration block. Add node, clear nodes, or
    close block. If multiple nodes are added in one declaration block, connect
    them. Update geometryInfo['nodes'] and geometryInfo['filaments']
    appropriately.
    """
    splitLine = re.split(',|\)|\(', line.strip())
    
    openSegment = self.segments[self._filamentNames.index(self._openFilament)]
      
    if splitLine[0] == '}':
      self._openFilament = None
    elif splitLine[0] == 'pt3dclear':
      openSegment.clear()
    elif splitLine[0] == 'pt3dadd':
      if not (len(splitLine) == 6 or
              (len(splitLine) == 7 and splitLine[-2] == '0')):
        raise IOError('Unexpected form for pt3dadd')
      x,y,z,d = tuple(float(s) for s in splitLine[1:5])
      if d <= 0:
        if d == 0:
          raise ValueError('pt3dadd with diameter = 0.0')
        else:
          raise ValueError('pt3dadd with diameter < 0.0')
        
      self._addNode(openSegment, x, y, z, 0.5 * d)

      if len(openSegment.nodes) > 1:
        node0 = openSegment.nodes[-2]
        node1 = openSegment.nodes[-1]
        self._addCompartment(openSegment, node0, node1, append=True)
    else:
      raise IOError('Invalid filament command')


  def _addConnection(self, splitLine):
    """
    Return dict describing connection between two filaments
    """
    (name1, location1) = re.split('\(|\)', splitLine[1])[0:2]
    (name2, location2) = re.split('\(|\)', splitLine[2])[0:2]
    
    connection = { \
      'filament1' : name1, \
      'location1' : float(location1), \
      'filament2' : name2, \
      'location2' : float(location2) \
    }
    self._connections.append(connection)
    self.connections.append(connection)


  def _createFilaments(self, splitLine):
    """
    Add requested number of filaments to geometry, as segments
    """
    if '[' and ']' in splitLine[1]:
      # hoc produced by Imaris, requests variable number of filaments
      # create baseName[numFilaments]
      baseName, numFilamentsStr = re.split('\[|\]', splitLine[1])[0:2]
      numFilaments = int(numFilamentsStr)
      if self._filamentNameType in [None, 'Imaris']:
        self._filamentNameType = 'Imaris'
      else:
        self._filamentNameType = 'Mixed'
        warn('Filament index will not reliably match numbers in filament name')
      thisType = 'Imaris'
    else:
      # hoc produce by Amira, requests 1 filament
      # create baseName
      baseName, numFilaments = splitLine[1], 1
      if self._filamentNameType in [None, 'Amira']:
        self._filamentNameType = 'Amira'
      else:
        self._filamentNameType = 'Mixed'
        warn('Filament index will not reliably match numbers in filament name')
      thisType = 'Amira'
    
    for n in range(numFilaments):
      if thisType == 'Imaris':
        name = '%s[%d]' % (baseName, n)
      else:
        name = baseName
      if name in self._filamentNames:
        raise IOError('%s already created' % name)
      newSeg = self._addSegment(name)
      newSeg.filamentIndex = len(self._filamentNames)
      self._filamentNames.append(name)
      self._filaments[newSeg.filamentIndex] = newSeg


  def _connectFilaments(self):
    """
    Loop through requested filament _connections.
    For each connection connect two filaments together by joining the nodes at
      their ends. Note that this removes a node for each connection
    """
    def _getSegmentFromFilament(_filament):
      _segment = self._filaments[self._filamentNames.index(_filament)]
      return _segment      
    
    while self._connections:
      connection = self._connections.pop()
      
      # get the filaments and locations
      location0 = connection['location2']
      filament0 = connection['filament2']
      segment0 = _getSegmentFromFilament(filament0)

      location1 = connection['location1']
      filament1 = connection['filament1']
      segment1 = _getSegmentFromFilament(filament1)
      
      self._connectSegments(segment0, location0, segment1, location1)
    
    self.nodes = [n for n in self.nodes if n not in self._removeNodes]
    self._removeNodes = set()

  
  def getFilamentIndex(self, seg):
    """
    Return index to filament from original .hoc file
    """
    #filamentIndex = int(seg.name.split('[')[1].split(']')[0])
    return seg.filamentIndex
    
  
  def getFilament(self, index):
    """
    return a segment based upon filament number
    """
    return self._filaments[index]
      


###############################################################################
def _parseArguments():
  arguments = sys.argv
  
  if len(arguments) != 2:
    print(_usageStr)
    raise TypeError('Incorrect number of arguments.')
  
  geoFile = arguments[1]
  return geoFile


def suggestProps(geometry, tau_m=215.2, tau_1=29.0, R_0=2.1, R_in=9.7):
  # R in MOhm, Tau in ms
  Cm = 1.0 # uF/cm^2, gospel
  Rm = 1.0e-3 * tau_m / Cm  # MOhm cm^2
  g1 = 1.0e-6 / Rm # S/cm^2
  RTau1 = 1.0e-3 * tau_1 / Cm
  g2 = 1.0e-6 / RTau1
  
  somaArea = geometry._soma.surfaceArea * 1.0e-2
  cellArea = geometry.surfaceArea * 1.0e-2
  
  g3 = 1.0e-6 / R_0 / somaArea
  g4 = 1.0e-6 / R_0 / cellArea
  g5 = 1.0e-6 / R_in / somaArea
  g6 = 1.0e-6 / R_in / cellArea
  
  print('Potential leak conductances:')
  print(g1, g2, g3, g4, g5, g6)
  
  """
  L = math.pi / math.sqrt(tau_m / tau_1 - 1) #unitless ratio of length / lambda
  RInf = R_in * math.tanh(L)
  print('Rm = %g MOhm cm^2' % Rm)
  print('L = %g' % L)
  print('RInf = %g MOhm' % RInf)
  # get d (equivalent diameter in cm) somehow...
  # maybe guess it's the smallest diameter in the primary neurite?
  somaIndex, somaPos = geometry.getSomaIndex()
  d = 2 * geometry.segments[somaIndex].minRadius * 1.0e-4 # cm
  
  Ri = (RInf * math.pi / 2 * d**1.5)**2 / Rm  # MOhm cm
  Ra = 1.0e6 * Ri # Ohm cm
  print('cm = %g uF/cm^2, gLeak = %g S/cm^2, Ra = %g Ohm cm'
        % (Cm, g, Ra))
  """
  

###############################################################################

def get_slopes(pts):
  # return the fractal dimension cumulative mass-fit slopes for zipped masses
  slopes, r_vals = [], []
  # use 2-point slopes
  for p in range(len(pts)-1):
    slopes.append((pts[p][1]-pts[p+1][1])/(pts[p][0]-pts[p+1][0]))
  fig3 = plt.figure()
  ax3 = fig3.add_subplot(111)
  ax3.scatter([pts[i][0] for i in range(len(pts)-1)],slopes)
  ax3.set_ylabel('Slope')
  ax3.set_xlabel('log Radius')
  plt.show()



def dist(pt1, pt2):
  #print('pt2 passed to dist ', type(pt2), pt2)
  return sqrt((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)

def cumulative_mass(geometry, trees):
  """
  Calculate the fractal dimension of each hemisphere
  """
  treecoords = dict.fromkeys(trees.keys(),[])
  listcoords = dict.fromkeys(trees.keys(),[])
  print('Populating tree coordinates for fractal dimension...')
  for s in geometry.segments:
    for k in trees.keys():
      if s.name in trees[k]:
        listcoords[k].append([tuple([n.x,n.y,n.z]) for n in s.nodes]) # or list
  for k in listcoords.keys():
    for s in listcoords[k]:
      for n in s:
        treecoords[k].append(n)
    #print('shape of %s is ' %k)
    #print(treecoords[k][:10])
  
  print('Populating the distance list (repeated calls to dist func)')
  fracdim = dict.fromkeys(treecoords.keys(),[])
  # populate the distances
  count = 0
  for k in treecoords.keys(): # for each tree
    treecoords[k] = np.array(treecoords[k]) # convert to np array for norm
    for n1 in treecoords[k]: # for each node in each tree
      for n2 in treecoords[k]:
        #print('n1 ' + n1 + ' n2 ' + n2)
        if (n1 != n2).all():
          # print(n1, n2)
          #fracdim[k].append(np.linalg.norm(treecoords[k][n1]-
          #                                 treecoords[k][n2]))
          fracdim[k].append(dist(n1,n2))
      if count%100==0:
        print('%i / %i' %(count, len(treecoords[k])))
      count = count + 1
  
  # save progress as pickle dictionary
  pickle.dump(fracdim, open('current_fracdim.p', 'wb'))
  
  print('Fracdim dictionary populated. Creating histogram')
  # populate the hists
  hists = dict.fromkeys(fracdim.keys(),[])
  bin_edges = None
  for k in fracdim.keys():
    if not bin_edges:
      print('k is %s and contains' %k)
      print(fracdim[k])
      hists[k], bin_edges = np.histogram(fracdim[k], bins=100)
    else:
      hists[k], _ = np.histogram(fracdim[k], bin_edges)
  # make hists cumulative
  masses = dict.fromkeys(treecoords.keys(),[])
  for k in hists.keys():
    masses[k] = [sum(hists[k][:n]) for n in range(len(hists[k])+1)] # could use cumsum here instead
  
  # slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
  
  print('Plotting binned cumulative masses')
  # plot findings
  bin_s = bin_edges # [(bin_edges[i]+bin_edges[i+1])/2 for i in range(len(bin_edges)-1)]
  Cs = ['b','r','g','y','k'] # color list
  fig3 = plt.figure()
  ax3 = fig3.add_subplot(111)
  count = 0
  for k in masses.keys():
    print(np.shape(bin_s), np.shape(masses[k]))
    ax3.scatter(np.log(bin_s), np.log(masses[k]), c=Cs[count])
    count = count+1
  ax3.set_ylabel('log Masses')
  ax3.set_xlabel('log Radius')
  plt.show()
  
  # save progress as pickle dictionary
  pickle.dump(masses, open('current_masses.p','wb')) 
  
  """go = input('Print masses and bins? (y/n)')
  if go == 'y' or go == 'Y':
    print('bins:')
    print(np.log(bin_s))
    print('masses:')
    for k in masses.keys():
      print(np.log(masses[k]))"""
  
  for k in masses.keys():
    pts = list(zip(np.log(bin_s), np.log(masses[k])))
    get_slopes(pts)
  
  return



def hemispheres(geometry, show=False):
  """
  This function runs neighbor-generation analysis to assign segments 
  to one sub-graph. Disputes are settled by either (1) lower generation
  wins, (2) similar radius wins, or (3) chance
  """
  trees, treecoords = {}, {}
  for a in geometry._axons:
    trees[a.name] = [a.name]
    treecoords[a.name] = [(a.coordAt(0)[0],a.coordAt(0)[1],a.coordAt(0)[2]),
                          (a.coordAt(0.5)[0],a.coordAt(0.5)[1],a.coordAt(0.5)[2]),
                          (a.coordAt(1)[0],a.coordAt(1)[1],a.coordAt(1)[2])]
    print('Axon neighbors are: %s' %[l.name for l in a.neighbors])
  # unassigned = [s.name for s in geometry.segments]
  
  # create a neighbors dict for easier usage
  nebs, coords = {}, {}
  for s in geometry.segments:
    nebs[s.name] = [l.name for l in s.neighbors]
    coords[s.name] = [(s.coordAt(0)[0],s.coordAt(0)[1],s.coordAt(0)[2]),
                      (s.coordAt(0.5)[0],s.coordAt(0.5)[1],s.coordAt(0.5)[2]),
                      (s.coordAt(1)[0],s.coordAt(1)[1],s.coordAt(1)[2])]
    if s.name == 'filament_100000028[1075]':
      print('Adding axon neighbors and coords')
  print('tree keys: %s' %(trees.keys()))
  
  #print(nebs.keys())
  
  stop=0
  while len(nebs.keys()) > 1 and stop < 10000000: # as long as there are unassigned elements
    for k in trees.keys(): # add the neighbors one generation at a time
                           # for each tree,
      trees[k] = list(set(trees[k]))
      for seg in trees[k]: # find neighbors for each seg
        # print('seg is %s' %str(seg))
        if seg in nebs.keys():
          for n in range(len(nebs[seg])):
            trees[k].append(nebs[seg][n])
          for n in range(len(coords[seg])):
            treecoords[k].append(coords[seg][n])

          # print('added %s' %seg)
          del coords[seg] # delete the used segs
          del nebs[seg]
  
          break
          
        else:
          #print('seg %s not found in nebs %s.' %(str(seg), k))
          stop = stop +1
          #break
    # print('%i segments remain' %len(nebs.keys()))
  
  for t in trees.keys():
    trees[t] = set(trees[t])
    treecoords[t] = set(treecoords[t])
    print('Tree %s has %i segments' %(str(t),len(trees[t])))
  
  # plotting stuff -- this takes super long, best to skip for debugging
  if show==True:
    count = 0
    Cs = ['b','r','y','k','g']
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    for t in treecoords.keys():
      prog = 0
      #treecoords[t] = signal.decimate(treecoords[t],10)
      for p in treecoords[t]:
        ax2.scatter(p[0],p[1],p[2], c=Cs[count], marker='.',
                    edgecolor=Cs[count])
        prog = prog+1
        if prog%10000 == 0:
          print('Progress: %i / %i' %(prog, len(treecoords[t])))
      count=count+1
    print('Displaying graph...')
    plt.show()
  
  return geometry, trees



def targetAxons(geometry, show=False):
  """
  This function asks the user to identify the axons and relies
  on Ted's axon-finding procedure to produce guesses.
  """
  # self.findAxon is a good place to start
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111, projection='3d')
  axCoords, axnames = [], []
  
  if len(geometry._axons) == 1: # if only one axon, just keep it
    show = False
    
  else:
  
    for a in geometry._axons:
      axCoords.append([a.coordAt(0), a.coordAt(0.5), a.coordAt(1)])
      axnames.append(a.name)
      if show == True:
        ax1.scatter(a.coordAt(0)[0], a.coordAt(0)[1], a.coordAt(0)[2],
                    color = 'r', marker = '^', s=50)
        ax1.scatter(a.coordAt(0.5)[0], a.coordAt(0.5)[1], a.coordAt(0.5)[2],
                    color = 'r', marker = '^', s=50)
        ax1.text(a.coordAt(0.5)[0], a.coordAt(0.5)[1], a.coordAt(0.5)[2], 
                     axnames.index(a.name), 'x')
        ax1.scatter(a.coordAt(1)[0], a.coordAt(1)[1], a.coordAt(1)[2],
                    color = 'r', marker = '^', s=50)
    if show == True:
      for s in geometry.segments:
        ax1.scatter(s.coordAt(0)[0], s.coordAt(0)[1], s.coordAt(0)[2],
                    color='k', marker='.')
        ax1.scatter(s.coordAt(0.5)[0], s.coordAt(0.5)[1], s.coordAt(0.5)[2],
                    color='k', marker='.')
        ax1.scatter(s.coordAt(1)[0], s.coordAt(1)[1], s.coordAt(1)[2],
                    color='k', marker='.')
      
      ax1.scatter(geometry.soma.coordAt(0.5)[0], geometry.soma.coordAt(0.5)[1],
                  geometry.soma.coordAt(0.5)[2], c='b',marker='*', s=50)
      ax1.text(geometry.soma.coordAt(0.5)[0], geometry.soma.coordAt(0.5)[1],
              geometry.soma.coordAt(0.5)[2], 'soma', 'x')
      print('Soma is: %s' %str(geometry.soma.name))
      
      plt.show()
  
    print('Found %i axons:' %len(geometry._axons))
    print(axnames)
    
    which_ax = input('Which axons to keep? Separate with spaces. Return to exit. ')
    which_ax = [int(i) for i in which_ax.split(None)]
    if len(which_ax) > 0:
      geometry._axons = [geometry._axons[i] for i in which_ax]
      geometry, trees = hemispheres(geometry)
    else:
      print('Done')
  geometry, trees = hemispheres(geometry)
  
  return geometry, trees
  



def makeConnectionsNames(geometry):
  # make filament blah blah into names to refer to segments.name
  name_connections = []
  for con in geometry.connections:
    name_connections.append([con['filament2'], con['filament1']])
  geometry.name_connections = name_connections
  
  return geometry



def paintPlots(geometry):
  # graphical representation of hemisphere analysis
  # get the nodes for each branchpoint and make the edges
  edge_from = [] # two nodes indices define the edge
  nodes = {} # nodes are hashed with integers in a dict
  geometry = makeConnectionsNames(geometry)
  for s in geometry.segments:
    coord1 = s.coordAt(0)
    coord2 = s.coordAt(1)
    if coord1 not in nodes.values():
      nodes[len(nodes.keys())] = coord1
    if coord2 not in nodes.values():
      nodes[len(nodes.keys())] = coord2
      
  for con in geometry.name_connections:
    #print(con)
    c1, c2 = con[0], con[1]
    for s in geometry.segments:
      if s.name == c1:
        coord2 = s.coordAt(0)
      if s.name == c2:
        coord1 = s.coordAt(1)
    # found the relevant coords (nodes)
    n1, n2 = None, None
    for k in nodes.keys():
      if nodes[k] == coord1:
        n1 = k
      if nodes[k] == coord2:
        n2 = k
    # found indexed nodes
    edge_from.append([n1, n2])
  # finished indexing nodes and connections
  
  G = nx.Graph()
  G.add_edges_from(edge_from)
  print('Found %i edges in the neuron ' %len(G.edges()))
  #nx.draw_shell(G, node_size=0.1)
  #plt.show()
  
  
  return
  
  
  
  """
  # start @ soma
  for n in geometry.soma.nodes:
    print(n.x, n.y, n.z)
  
  if len(geometry.soma.neighbors) > 1
  for n in geometry.soma.neighbors:
  """
  
  #print(type(geometry.soma.nodes))
  


def demoRead(geoFile):
  ### Read in geometry file and pre-compute various quantities
  geometry = HocGeometry(geoFile)
  
  tips, tipPositions = geometry.getTips()
  pDF = PathDistanceFinder(geometry, geometry.soma)
  #tortuosities = [geometry.pathTortuosity(pDF.pathTo(tip, pos))
  #                for tip, pos in zip(tips, tipPositions)]
  tortuosities = [pDF.tortuosityTo(tip, pos)
                  for tip, pos in zip(tips, tipPositions)]
  meanTort, stdTort = mean(tortuosities), std(tortuosities)
  
  print('From soma to tips, tortuosity is %.1f +- %.1f'
        % (meanTort, stdTort))
        
  Cons =  geometry.connections
  Seg1s, Seg2s = [], []
  for c in Cons:
    Seg1s.append(c['filament1']) # here, location1 is always 0
    Seg2s.append(c['filament2']) # here, location2 is always 1
    #geometry.c['filament1'].coordAt(c['location1'])
  
  paintPlots(geometry)
  geometry, trees = targetAxons(geometry)
  cumulative_mass(geometry, trees)
  
  
    
    


###############################################################################
if __name__ == "__main__":
  # get the geometry file
  geoFile = _parseArguments()
  # run a demo of capabilities
  demoRead(geoFile)
  #exit
  sys.exit(0)
