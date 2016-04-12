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
import math



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
def demoReadOld(geoFile):
  ### Read in geometry file and pre-compute various quantities
  geometry = HocGeometry(geoFile)

  ### This section gets indices to the Axon, Soma, and neurite tips,
  ### then measures distance from Axon to Soma, and from Axon to a random
  ### neurite tip.
  axons = geometry.findAxons()

  if axons:
    axonInds, axonTipPos = geometry.getAxonIndices()
    print('Axon tip length = %g' % axons[0].length)
    print('Axon tip pos = %f' % axonTipPos[0])
    pDF = PathDistanceFinder(geometry, axons[0], axonTipPos[0])
    soma = geometry.soma
    somaInd, somaPos = geometry.getSomaIndex()
    print('Path distance from Axon tip to Soma = %g' %
          pDF.distanceTo(soma, somaPos))

    import random
    random.seed
    
    neuriteTipInd, neuriteTipPos = geometry.getTipIndices()
    whichTip = random.randint(0, len(neuriteTipInd) - 1)
    randomTip, randomPos = neuriteTipInd[whichTip], neuriteTipPos[whichTip]
    print('Path distance from Axon tip to random neurite tip = %g' %
          pDF.distanceTo(geometry.getFilament(randomTip), randomPos))

  else:
    print('No axons found')
  
  
  if os.access('steady_voltages.pickle', os.R_OK):
    import cPickle
    from math import sqrt, log
    from matplotlib import pyplot

    pDF = PathDistanceFinder(geometry, somaInd)
    dSeg = [pDF.distanceTo(s) for s in geometry.segments]
    with open('steady_voltages.pickle', 'r') as fIn:
      vSteady = cPickle.loads(fIn.read())
    eLengths = pDF.getElectrotonicLengths(vSteady)
    logELengths = [log(eL) for eL in eLengths]
    dist = [pDF.distanceTo(s) for s in range(len(vSteady))]
    diam = [sqrt(s.avgRadius) for s in geometry.segments]
    """
    pyplot.figure()
    pyplot.plot(dist, vSteady, 'r.')
    pyplot.figure()
    pyplot.plot(dist, logELengths, 'b.')
    pyplot.figure()
    pyplot.plot(diam, logELengths, 'b.')
    pyplot.show()
    """
  
    suggestProps(geometry)
  
  ### Display summary info
  geometry.displaySummary()



def dist(pt1, pt2):
  return math.sqrt((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2)


def taper(geometry):
  # Do the taper thing
  tapers, coords = [], []
  first, last = None, None
  for s in geometry.segments:
    for n in s.nodes:
      if not first:
        first = n.r1
      else:
        last = n.r1
    tapers.append((first-last)/s.length)
    coords.append((s.coordAt(0.5)[0], s.coordAt(0.5)[1], s.coordAt(0.5)[2]))
  # tapers compiled
  
  mean_pt = ( sum([i[0] for i in coords])/len(coords), 
              sum([i[1] for i in coords])/len(coords),
              sum([i[2] for i in coords])/len(coords) )
  dists = [dist(i,mean_pt) for i in coords]
  
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111, projection='3d')
  ax1.scatter([i[0] for i in coords], [i[1] for i in coords],
              [i[2] for i in coords], c=dists, edgecolor=None)
  plt.show()
  
  return



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
  
  taper(geometry)
  
    
    


###############################################################################
if __name__ == "__main__":
  # get the geometry file
  geoFile = _parseArguments()
  # run a demo of capabilities
  demoRead(geoFile)
  #exit
  sys.exit(0)
