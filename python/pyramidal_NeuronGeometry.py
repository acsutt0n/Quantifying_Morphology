#!/usr/bin/python



_usageStr=\
"""usage: neuron_readExportedGeometry.py geoFile
  get dictionary object describing neuron model geometry info by reading file
"""


import os
from scipy import special, mean, std
from collections import deque
import matplotlib.pyplot as pyplot
from math import log, sqrt, atan, isnan, pi, acos
from bisect import bisect_left
import numpy as np

"""
Geometry class public methods: (self is always first argument)
 setFileName(_fileName)
 numCompartments()
 readGeometry()  pure virtual
 displaySummary()
 findBranches()
 checkConnectivity()
 shollAnalysis()
"""

terminalColors = {
  'endColor'   : '\033[0m',

  'black'      : '\033[0;30m',
  'red'        : '\033[0;31m',
  'green'      : '\033[0;32m',
  'yellow'     : '\033[0;33m',
  'blue'       : '\033[0;34m',
  'purple'     : '\033[0;35m',
  'cyan'       : '\033[0;36m',
  'lightGray'  : '\033[0;37m',

  'darkGray'   : '\033[1;30m',
  'boldRed'    : '\033[1;31m',
  'boldGreen'  : '\033[1;32m',
  'boldYellow' : '\033[1;33m',
  'boldBlue'   : '\033[1;34m',
  'boldPurple' : '\033[1;35m',
  'boldCyan'   : '\033[1;36m',
  'white'      : '\033[1;37m'
}


def warn(warnStr, extraInfo='', color='boldRed'):
  if extraInfo:
    print(terminalColors[color] + warnStr + terminalColors['endColor']
          + ': ' + extraInfo)
  else:
    print(terminalColors[color] + warnStr + terminalColors['endColor'])


def cumsum(values, start=0.0, yieldStart=True):
  """
  Return generator to running cumulative sum of values
  """
  if yieldStart:
    yield start
  
  for v in values:
    start += v
    yield start


"""
class PathDistanceFinder
Finds network distances from one segment/branch at specified position
Can report distance to another segment/branch at specified position via
  .distanceTo()
Can report optimal path to another segment/branch at specified position va
  .pathTo()
Can report tortuosity of optimal path via
  .tortuosityTo()
Can compute electrotonic lengths from a list of voltages via
  .getElectrotonicLengths()
"""
class PathDistanceFinder(object):
  def __init__(self, geometry, segment, pos=0.5, warnLoops=False):
    self.geometry = geometry
    if type(segment) == int:
      segment = geometry.segments[segment]
      self.network = geometry.segments
    elif segment in geometry.segments:
      self.network = geometry.segments
    elif segment in geometry.branches:
      self.network = geometry.branches
    else:
      raise TypeError('Start segment must be an index to geometry.segments, or'
                      ' an object from geometry.segments or geometry.branches')
    
    self.warnLoops = warnLoops
    self.startSegment = segment
    self.startPos = pos
    self.startCoord = segment.coordAt(pos)
    self.distances, self.branchOrders = self._computeDistances()
    
  
  def distanceTo(self, segment, pos=0.5):
    # return distance to specified segment at the specified location
    if type(segment) == int:
      segment = self.network[segment]
    
    if segment not in self.distances:
      raise KeyError('%s is not reachable from the network' % segment.name)

    return min(baseD + segment.length * abs(pos - startPos) for
               baseD, startPos, pathDesc, path in self.distances[segment])
  
  def pathTo(self, segment, pos=0.5):
    # return optimal path to specified segment at specified location
    return min(((baseD + segment.length * abs(pos - startPos), path) for
               baseD, startPos, pathDesc, path in self.distances[segment]),
               key=lambda x:x[0])[1]
  
  def pathDescriptionTo(self, segment, pos=0.5):
    # return optimal path to specified segment at specified location
    return min(((baseD + segment.length * abs(pos - startPos), pathDesc) for
               baseD, startPos, pathDesc, path in self.distances[segment]),
               key=lambda x:x[0])[1]
  
  def tortuosityTo(self, segment, pos=0.5):
    stopCoord = segment.coordAt(pos)
    euclideanD = sqrt(sum((r0-r1)**2 for r0,r1 in zip(self.startCoord,
                                                      stopCoord)))
    pathD = self.distanceTo(segment, pos)
    tortuosity = pathD / euclideanD
    if tortuosity < 1.0:
      warn('Tortuosity < 1',
           self.pathDescriptionTo(segment, pos))
      raise RuntimeError('Path to %s at %g has tortuosity < 1'
                         % (segment.name, pos))
    return tortuosity
  
  
  def branchOrder(self, segment):
    return self.branchOrders[segment]
  
  
  def getElectrotonicLengths(self, steadyVoltages):
    # given a list of steady-state voltages, return electrotonic lengths of all
    #  segments
    return [self._getElectrotonicLength(ind, seg, steadyVoltages)
            for ind, seg in enumerate(self.network)]
  
  
  def _getElectrotonicLength(self, index, segment, steadyVoltages):
    # given list of steady-state voltages, compute the electrotonic length of a
    #  specific segment
    dSeg = self.distanceTo(segment)
    vSeg = steadyVoltages[index]
    segOrder = self.branchOrders[segment]

    # use only one other neighbor, pick the closest one that is of a different
    # branch order
    
    neighbor = min((s for s in segment.neighbors
                    if self.branchOrders[s] != segOrder),
                   key=lambda n: abs(self.distanceTo(n) - dSeg))
    dNeighbor = self.distanceTo(neighbor)
    vNeighbor = steadyVoltages[self.network.index(neighbor)]
    eLength = (dNeighbor - dSeg) / log(vSeg / vNeighbor)
    if eLength < 0:
      nIndex = self.network.index(neighbor)
      nOrder = self.branchOrders[neighbor]
      print('eLength=%g: ind=%d, nInd=%d, d=%g, nd=%g, v=%g, nV=%g, '
            % (eLength, index, nIndex, dSeg, dNeighbor, vSeg, vNeighbor)
            + 'order=%d, nOrder=%d' % (segOrder, nOrder))
    return eLength
      
  
  def _computeDistances(self):
    # use Dijkstra's algorithm to find path distance from start to rest of
    # neuron.
    # Keep track of effect of startPos (starting position in startSegment)
    # Also keep track of effect of pos of each final segment
    segment, startPos = self.startSegment, self.startPos
    # distances is a dict object, with segments as keys
    # the values are a list of paths, each path described by a tuple
    #   (pathDistance, connecting location in the segment, path description,
    #    list of segments in path )
    distances = { segment : [(0.0, startPos,
                              segment.name + '(%.1f)' % startPos, [segment])] }
    branchOrders = { segment : 0 }
    openSegments = { segment, }
    while openSegments:
      segment = openSegments.pop()
      for currentD, startPos, segPathDesc, segPath in distances[segment]:
        branchOrderInc = 1
        if branchOrders[segment] > 0 and len(segment.neighbors) <= 2:
          branchOrderInc = 0
        #branchOrderInc = int(len(segment.neighbors) > 2)
        for neighbor, (connectLoc, nConnectLoc, node) \
            in zip(segment.neighbors, segment.neighborLocations):
          pathD = currentD + segment.length * abs(startPos - connectLoc)
          # check if neighbor in distances?
          if neighbor not in distances:
            # found path to new segment
            nPath = segPath + [neighbor]
            nPathDesc = segPathDesc + '->(%.1f)->' % connectLoc + \
                  neighbor.name + '(%.1f)' % nConnectLoc
            distances[neighbor] = [(pathD, nConnectLoc, nPathDesc, nPath)]
            openSegments.add(neighbor)
            branchOrders[neighbor] = branchOrders[segment] + branchOrderInc
          else:
            # Either there is a loop involving this segment, or the path is
            #  backtracking
            # Check if the current path is an efficient route to the loop
            efficient = True
            insertInd = None
            loopDistances = distances[neighbor]
            for ind, (loopD, loopPos, loopPathDesc, loopPath) in \
                enumerate(loopDistances):
              traverse = neighbor.length * abs(loopPos - nConnectLoc)
              if pathD >= loopD + traverse:
                # the new path to neighbor is too slow to ever be useful
                efficient = False
                break
              elif pathD + traverse < loopD:
                # the new path to neighbor renders an old one obsolete
                loopDistances.pop(ind)
              if insertInd is None and pathD < loopD:
                insertInd = ind
            if efficient:
              nPath = segPath + [neighbor]
              nPathDesc = segPathDesc + '->(%.1f)->' % connectLoc + \
                    neighbor.name + '(%.1f)' % nConnectLoc
              pathInfo = (pathD, nConnectLoc, nPathDesc, nPath)
              if insertInd is None:
                loopDistances.append(pathInfo)
              else:
                loopDistances.insert(insertInd, pathInfo)
              distances[neighbor] = loopDistances
              openSegments.add(neighbor)
              branchOrders[neighbor] = branchOrders[segment] + branchOrderInc
              if self.warnLoops and len(loopDistances) > 1:
                warn('%d efficient paths to %s.' % (len(loopDistances), neighbor.name))
                for ind, loopPos, loopPath in loopDistances:
                  print(loopPath)
    
    return distances, branchOrders



class Geometry:
  def __init__(self, _fileName = None):
    # who knows, do something?
    self.path = None
    self.fileName = None
    self.name = None
    
    # store the geometry info here
    self.nodes = []
    self.segments = []
    self.branches = []
    self.compartments = []
    self.branchOrders = None
    self.tags = {'*' : 0}
    
    # helper sets for efficient deleting
    self._removeNodes = set()
    self._removeSegments = set()
    # keep track of which objects have had connectivity checked
    self._connectivityChecked = set()
    
    self._soma = None
    self._somaBranch = None
    self._axons = None
    self._axonsBranch = None
    
    self.surfaceArea = 0.0   # mm^2
    self.volume = 0.0        # mm^3

    self.minRange = [float('nan'), float('nan'), float('nan')]
    self.maxRange = [float('nan'), float('nan'), float('nan')]
    
    if _fileName is not None:
      self.setFileName(_fileName)
      self.readGeometry()
  
  
  def setFileName(self, fileName):
    self.path = os.path.dirname(os.path.abspath(fileName))
    self.fileName = fileName
    self.name = os.path.basename(fileName).split('.')[0]
  
  def numCompartments(self):
    return len(self.compartments)
  
  
  def readGeometry(self):
    raise RuntimeError( \
      'Geometry must be a subclass that knows how to read files')  
  
  
  def displaySummary(self):
    """
    Display summary statistics of neuron geometry
    """
    print("total number of nodes: %d" % len(self.nodes))
    print("total number of compartments: %d" % len(self.compartments))
    print("total number of segments: %d" % len(self.segments))

    subGraphs = self.checkConnectivity(removeDisconnected=False, debugInfo=True)
    
    print("number of connected nodes: %d" % len(self.nodes))
    print("number of connected compartments: %d" % len(self.compartments))
    print("number of connected segments: %d" % len(self.segments))
    
    self.findBranches()
    print("number of branches: %d" % len(self.branches))
    
    soma = self.soma
    somaArea = sum(c.surfaceArea for c in soma.compartments \
                   if 'Soma' in c.tags)
    print('Soma Area = %g mm^2' % somaArea)
    print('Found %d axon%s' % (len(self._axons), 's'*(len(self._axons)!=1)))
    print("volume: %g mm^3" % self.volume)
    print("surface area: %s mm^2" % self.surfaceArea)
    self.calcBranchOrder(doPlot=False)
    self.shollAnalysis(straightenNeurites=True)
    self.mergeBranchesByDistanceToEdge()
    pyplot.show() ##### SHOW PLOTS!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  #############################################################################
  def getProperties(self, passiveFile="", display=False, # CHANGED FALSE
                    makePlots=False):
    def _dispListStats(L, confidence = 0.05, display=True, printName=""):
      # return median, lowBound, highBound
      sortedL = sorted(L)
      numL = len(L)
      if numL % 2:
        # odd number of elements
        medianInd = int((numL - 1) / 2)
        median = sortedL[medianInd]
      else:
        medianInd = numL / 2
        median = (sortedL[int(medianInd)] + sortedL[int(medianInd - 1)]) / 2.0
      lowInd = int(round( 0.5 * confidence * numL ))
      highInd = int(round( (1.0 - 0.5 * confidence) * numL ))
      
      low = sortedL[lowInd] ; high = sortedL[highInd]
      
      if display and printName:
        print('%s = %.2f +%.2f -%.2f'
              % (printName, median, high - median, median - low))
    
    def _plotTraces(timeTrace, vTraces):
      from scipy import array
      fig = pyplot.figure()
      axes = fig.add_subplot(1,1,1)
      y = array(vTraces.values()).transpose()
      axes.plot(timeTrace, y)
      pyplot.ylabel('Membrane Potential (mV)')
      pyplot.xlabel('Time (ms)')
      pyplot.title('Model Response to Step Current')
      pyplot.tight_layout()
    
    from scipy.optimize import brentq, fmin
    def _rallLaw(p, *ratios):
      try:
        return sum(r**p for r in ratios) - 1.0
      except OverflowError:
        return float('inf')
    def _rallLawTrouble(p, *ratios):
      try:
        return (sum(r**p for r in ratios) - 1.0)**2
      except OverflowError:
        return float('inf')
        
    def _getRallPow(parentR, daughterRs):
      ratios = tuple(d / parentR for d in daughterRs)
      checkPows = [-log(len(ratios)) / log(r) for r in ratios if r != 1.0]
        
      try:
        return brentq(_rallLaw, min(checkPows), max(checkPows), args=ratios)
      except ValueError:
        #print('Rall-incompatible branch ratios: %s'
        #      % ' '.join('%.2f' % r for r in ratios))
        return fmin(_rallLawTrouble, 0.0, args=ratios, disp=False)[0]
    
    def _overallRall(p, ratiosList):
      return sum(_rallLawTrouble(p, *ratios) for ratios in ratiosList)
    def _getOverallRallPow(ratiosList):
      return fmin(_overallRall, 0.0, args=(ratiosList,), disp=False)[0]
    
    # check connectivity
    self.checkConnectivity(removeDisconnected=True, removeLoops=True)
    self.findBranches()
    
    if display:
      print("number of connected nodes: %d" % len(self.nodes))
      print("number of connected compartments: %d" % len(self.compartments))
      print("number of connected segments: %d" % len(self.segments))
      print("number of branches: %d" % len(self.branches))
      print('Surface area = %g mm^2' % self.surfaceArea)
      print('Volume = %g mm^3' % self.volume)
      print('Surface to volume ratio = %g mm^-1'
            % (self.surfaceArea/self.volume))

    # make a path distance finder centered at the soma
    pDF = PathDistanceFinder(self, self.soma)
    # find all the neuron tips
    tips, tipPositions = self.getTips()
    # measure path lengths from Soma to tips
    pathLengths = [pDF.distanceTo(tip, pos)
                   for tip, pos in zip(tips, tipPositions)]
    _dispListStats(pathLengths, display=display,
                   printName='Path length from Soma to tips')
    # measure tortuosities from Soma to tips
    tortuosities = [pDF.tortuosityTo(tip, pos)
                    for tip, pos in zip(tips, tipPositions)]
    _dispListStats(tortuosities, display=display,
                   printName='Tortuosity of path from Soma to tips')
  
    # measure branch tortuosities
    bTortuosities = [branch.tortuosity for branch in self.branches
                     if branch.tortuosity < float('inf')]
    _dispListStats(bTortuosities, display=display,
                   printName='Tortuosity of neuron branches')
  
    if self.soma.branchOrder is None:
      self.calcBranchOrder(doPlot=False)

    self.mergeBranchesByDistanceToEdge(makePlots=makePlots)

    branchAngles = [getBranchAngle(branch, neighbor, segLoc, nLoc, node)
                    for branch in self.branches
                      for neighbor, (segLoc, nLoc, node)
                        in zip(branch.neighbors, branch.neighborLocations)
                        if neighbor.branchOrder > branch.branchOrder]
    _dispListStats(branchAngles, display=display,
                   printName='For neuron branches, branch angle')
    
    radList = [] # added 11.12.2014
    rallRatios = []
    daughterRatios = []
    rallPowers = []
    ratiosList = []
    DPratios = [] # added 11.12.2014
    daughterdaughter = [] # daughter-daughter ratio
    
    for segment in self.branches:
      #if segment.branchOrder < 4:
      #  continue
      daughters = [n for n in segment.neighbors
                   if n.branchOrder > segment.branchOrder]
      if daughters:
        rallRatio = \
          sum(n.avgRadius**1.5 for n in daughters) / segment.avgRadius**1.5
        rallRatios.append(rallRatio)
        daughterRatios.extend(n.avgRadius / segment.avgRadius
                              for n in daughters)
        # added 11.12.2014
        daughtRads, daughtCount = [], 0
        DDR = []
        for n in daughters:
          daughtRads.append(n.avgRadius)
          daughtCount = daughtCount + 1
          DDR.append(n.avgRadius)
        DPratios.append( (sum(daughtRads)/daughtCount)/segment.avgRadius )
        daughterdaughter.append(np.mean([DDR[d]/DDR[d+1] for d in 
                                         range(len(DDR)-1)]))
        
        
        rallPowers.append(_getRallPow(segment.avgRadius,
                                      [n.avgRadius for n in daughters]))
        ratiosList.append([n.avgRadius / segment.avgRadius for n in daughters])
        radList.append(segment.avgRadius)
        
    _dispListStats(rallRatios, display=display,
                   printName='For neuron branches, Rall ratio')
    _dispListStats(daughterRatios, display=display,
                   printName='For neuron branches, daughter branch ratio')
    _dispListStats(rallPowers, display=display,
                   printName='For neuron branches, Rall power')
  
    properties = {
      'Num Nodes' : len(self.nodes),
      'Num Compartments' : len(self.compartments),
      'Num Segments' : len(self.segments),
      'Num Branches' : len(self.branches),
      'Surface Area' : self.surfaceArea,
      'Volume' : self.volume,
      'Area-To-Volume Ratio' : self.surfaceArea / self.volume,
      'Path Length' : pathLengths,
      'Tortuosity' : tortuosities,
      'Branch Tortuosity' : bTortuosities,
      'Branch Angles' : branchAngles,
      'Rall Ratio' : rallRatios,
      'Daughter/Parent Radius' : daughterRatios,
      'Overall Rall Power' : _getOverallRallPow(ratiosList),
      'Radius List': radList,
      'DP Ratio': DPratios,
      'DD Ratio': daughterdaughter
    }
    units = {
      'Num Nodes' : '',
      'Num Compartments' : '',
      'Num Segments' : '',
      'Num Branches' : '',
      'Surface Area' : 'mm^2',
      'Volume' : 'mm^3',
      'Area-To-Volume Ratio' : 'mm^-1',
      'Path Length' : 'um',
      'Tortuosity' : '',
      'Branch Tortuosity' : '',
      'Branch Angles' : 'degrees',
      'Rall Ratio' : '',
      'Daughter/Parent Radius' : '',
      'Overall Rall Power' : '',
      'Radius List': 'um',
      'DP Ratio': '',
      'DD Ratio': ''
    }
  
    if passiveFile:
      from neuron_simulateGeometry import makeModel, simulateModel
      import peelLength
      import json
      # get the properties
      with open(passiveFile, 'r') as fIn:
        passiveProperties = json.load(fIn)
      # make a demo model
      model = makeModel(self, passiveProperties)
      # simulation model on specified geometry
      timeTrace, vTraces, textOutput = simulateModel(self, model)
      if makePlots:
        _plotTraces(timeTrace, vTraces)
       
      somaV = max(vTraces[self.soma.name])
      rIn = somaV / model['stimulus']['amplitude']
      properties['Input resistance'] = rIn
      units['Input resistance'] = 'MOhm'
      if display:
        print('Input resistance = %g MOhm' % rIn)
      tipsV = [max(vTraces[segment.name]) for segment in self.segments
              if 'Soma' not in segment.tags and segment.isTerminal]
      tipsTransfer = [tipV / somaV for tipV in tipsV]
      _dispListStats(tipsTransfer, display=display,
                     printName='Coupling coefficient from soma to tips')
      properties['Coupling Coefficient'] = tipsTransfer
      units['Coupling Coefficient'] = ''
      model, vErr, vResid = \
      peelLength.modelResponse(timeTrace, vTraces[self.soma.name],
                               verbose=False, findStepWindow=True,
                               plotFit=False, debugPlots=False,
                               displayModel=display)
      tauM = model[0][0]
      if display:
        print('membrane tau = %6.2f ms' % tauM)
      properties['Membrane Time Constant'] = tauM
      units['Membrane Time Constant'] = 'ms'
  
    if makePlots:
      self.shollAnalysis()
      
    return properties, units
    
  
  def findBranches(self):
    """
    Break up geometry into segments defined by branch points, starting at the
    soma
    """
    if self.branches:
      return
    
    self.checkConnectivity(removeDisconnected=True)
    if not self._somaBranch:
      self._findSoma()
    somaBranch, somaNeighbors0, somaNeighbors1 = self._somaBranch
    # This can cause weird errors if not fixed:
    if somaBranch.neighbors:
      warn('Some routine cleared self.branches without removing '
           + 'somaBranch.neighbors')
      somaBranch.neighbors = []
    
    self.branches = [somaBranch]
    
    openBranches = [(somaBranch, 0, somaNeighbors0), \
                     (somaBranch, 1, somaNeighbors1)]
    openCompartments = set(self.compartments).difference(
      somaBranch.compartments)
    
    while openBranches:
      # check an open branch to see if it has any neighbors branching off
      checkBranch, side, neighbors = openBranches.pop()
      
      # find neighbors on specified side of checkBranch
      commonNeighbors = { (checkBranch, side), }
      
      checkNode = checkBranch.nodes[-side]
      for segment, pos, compartment in neighbors:
        if compartment not in openCompartments:
          continue

        # for each neighboring segment, find the branch it's in, based on
        # compartment (ignore segment and pos)
        branch, neighbors0, neighbors1 = self._getBranch(compartment)
        
        # add that branch to geometry
        self.branches.append(branch)
        
        # remove the compartments in branch from openCompartments
        openCompartments.difference_update(branch.compartments)
        
        # add branch to dict of common neighbors, at appropriate side
        # and add to openBranches with appropriate neighbors
        if branch.nodes[-1] == checkNode:
          if branch.nodes[0] == checkNode:
            # branch is a loop with both ends connected to _check
            commonNeighbors.add( (branch, 0) )
            commonNeighbors.add( (branch, 1) )
          else:
            # side 1 of branch connects to checkBranch at checkNode
            if branch.nodes[-1] != checkNode:
              checkTags = ' '.join(checkBranch.tags)
              checkNodeInd = self.nodes.index(checkNode)
              branchTags = ' '.join(branch.tags)
              branchNodes =str(tuple(self.nodes.index(n) for n in branch.nodes))
              raise AssertionError(('Node mismatch. %s (with tags %s) should'
                + ' connects to %s (with tags %s) at node %d, but %s has '
                + 'nodes %s') % (checkBranch.name, checkTags, branch.name,
                                 branchTags, checkNodeInd, branch.name,
                                 branchNodes))

            commonNeighbors.add( (branch, 1) )
            # side 0 is still open
            openBranches.append((branch, 0, neighbors0))
        else:
          # side 0 of branch connects to checkBranch at checkNode
          if branch.nodes[0] != checkNode:
            checkTags = ' '.join(checkBranch.tags)
            checkNodeInd = self.nodes.index(checkNode)
            branchTags = ' '.join(branch.tags)
            branchNodes = str(tuple(self.nodes.index(n) for n in branch.nodes))
            raise AssertionError(('Node mismatch. %s (with tags %s) should'
              + ' connects to %s (with tags %s) at node %d, but %s has '
              + 'nodes %s') % (checkBranch.name, checkTags, branch.name,
                               branchTags, checkNodeInd, branch.name,
                               branchNodes))
            
          commonNeighbors.add( (branch, 0) )
          # side 1 is still open
          openBranches.append((branch, 1, neighbors1))
      
      # update the neighbors of _check and all the new branches
      while commonNeighbors:
        n1, n1Side = commonNeighbors.pop()
        for n2, n2Side in commonNeighbors:
          _makeNeighbors(n1, n2, n1Side, n2Side, checkNode)
  
  
  def _plotBranchStat(self, branchStat, yLabel, title, \
                      fontSize=22, barWidth=0.25):
    ### plot collected statistic along with number of branches
    ### branchStat should be a dictionary with:
    ###  -each key is a branch order (an integer)
    ###  -each item is a list of y-values (the list for all branches with
    ###   that branch order)
    order = list(branchStat.keys())
    order.sort()
    y = [branchStat[o] for o in order]
    x = list(range(len(order)))
    orderStr = [str(o) for o in order]
    numBranches = [len(y_n) for y_n in y]
    
    # make new figure
    fig = pyplot.figure()
    # plot number of branches as bar plot
    ax1 = pyplot.gca()
    pyplot.bar(x, numBranches, width=barWidth, color='g')
    pyplot.ylabel('# branches', fontsize=fontSize)
    pyplot.xlabel('Branch Order', fontsize=fontSize)
    pyplot.xticks(x, orderStr)
    
    # plot y statistics as a box and whisker plot
    positions = [x_n - barWidth/2.0 for x_n in x]
    ax2 = pyplot.twinx()
    pyplot.boxplot(y, positions=positions, widths=barWidth)
    pyplot.title(title, fontsize=fontSize)
    pyplot.ylabel(yLabel, fontsize=fontSize)
    pyplot.xlabel('Branch Order', fontsize=fontSize)
    pyplot.xticks(x, orderStr)
    
    # set the numBranches y-axis and labels on right, main on left
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position('right')
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position('left')
    pyplot.tight_layout()
    return fig


  def _plotBranchOrderStatistics(self):
    ### Visualize various statistics dependent upon branch order

    # collect data
    branchRadius = {}
    branchLength = {}
    for branch in self.mergedBranches:
      order = branch.centripetalOrder
      if order not in branchRadius:
        branchRadius[order] = [branch.maxRadius]
        branchLength[order] = [branch.length]
      else:
        branchRadius[order].append(branch.maxRadius)
        branchLength[order].append(branch.length)

    self._plotBranchStat(branchRadius, \
                         'Radius (um)', 'Branch Order vs. Radius')
    self._plotBranchStat(branchLength, \
                         'Length (um)', 'Branch Order vs Length')
    
  
  def calcBranchOrder(self, doPlot=True):
    self.calcForewardBranchOrder(doPlot=False, printAxonInfo=False)
    self.calcCentripetalOrder(doPlot=doPlot, network=self.branches)
    self.calcCentripetalOrder(doPlot=doPlot, network=self.segments)
  
  
  def calcForewardBranchOrder(self, doPlot=True, printAxonInfo=False):
    somaPos = self.soma.centroidPosition(mandateTag='Soma')
    pDF = PathDistanceFinder(self, self.soma, somaPos)
    for segment in self.segments:
      segment.branchOrder = pDF.branchOrder(segment)
    
    self.findBranches()
    somaPos = self.somaBranch.centroidPosition(mandateTag='Soma')
    pDF = PathDistanceFinder(self, self.somaBranch, somaPos)
    for branch in self.branches:
      branch.branchOrder = pDF.branchOrder(branch)
  
  
  def calcCentripetalOrder(self, doPlot=True, network=None):
    """
    Define neurite "ends" to be segments that have locally maximal branchOrder
    (ends are terminal segments unless there are loops). 
    Label each end with centripetal order 0.
    Every other segment is labeled with centripetal order equal to the length
    of the longest path from an end to that segment, provided that allowable
    paths ALWAYS move towards the soma.
    
    geometry.calcCentripetalOrder() sets segment.centripetalOrder set to this
    value for all segments in network
    """

    if network is None or not network:
      self.findBranches()
      network = self.branches

    def _isEnd(segment):
      return all(segment.branchOrder >= n.branchOrder
                 for n in segment.neighbors)
    ends = [segment for segment in network if _isEnd(segment)]
    for segment in ends:
      if 'Soma' in segment.tags:
        print('FUCK')
        print(segment.branchOrder)
        print([n.branchOrder for n in segment.neighbors])
        print(all(segment.branchOrder >= n.branchOrder
                 for n in segment.neighbors))

    for segment in network:
      segment.centripetalOrder = -1
    for segment in ends:
      # mark each end as having centripetal order 0
      segment.centripetalOrder = 0
      # now find paths from each end towards the soma. Each segment's
      # centripetal order is the length of the longest path from an end to the
      # segment, PROVIDED that the path ALWAYS moves towards the soma
      openSegs = [segment]
      while openSegs:
        currentSeg = openSegs.pop()
        neighborCentripOrder = currentSeg.centripetalOrder + 1
        for neighbor in currentSeg.neighbors:
          # look for new paths to soma, insisting that:
          #   1. The path ALWAYS moves closer to soma
          #   2. There isn't already a longer path through this area
          if neighbor.branchOrder < currentSeg.branchOrder and \
             neighbor.centripetalOrder < neighborCentripOrder:
            # this is a new (good) path, so mark the centripetal order and
            # continue it
            neighbor.centripetalOrder = neighborCentripOrder
            openSegs.append(neighbor)
  
  
  def checkConnectivity(self, removeDisconnected=False, checkObjects=None,
                        debugInfo=True, removeLoops=False):
    """
    Compute the connectivity of the network:
      -The number/members of connected subgraphs
      -The presence of any loops
    if removeDisconnected is True, remove all but largest subgraph from network
    
    Return the list of subgraphs
    """
    
    if checkObjects is None:
      checkObjects = self.segments
    checkHash = hash(str(checkObjects)+str(removeDisconnected))
    if checkHash in self._connectivityChecked:
      # don't need to check again
      return

    # check to be sure that neighborhood at a location/node is transitive
    for segment in checkObjects:
        for neighbor, (pos, nPos, node) in zip(segment.neighbors,
                                               segment.neighborLocations):
          assert segment in node.segments, \
            "%s should be in node %d's list of segments, but is not" \
            % (segment.name, self.nodes.index(node))
          assert neighbor in node.segments, \
            "%s should be in node %d's list of segments, but is not" \
            % (neighbor.name, self.nodes.index(node))
          for n2, (pos2, nPos2, node2) in zip(neighbor.neighbors,
                                              neighbor.neighborLocations):
            if node2 != node:
              continue
            assert n2 == segment or n2 in segment.neighbors, \
                   "%s and %s are neighbors at node %d, and so are %s and %s,"\
                   " but %s and %s are not" % (segment.name, neighbor.name,
                   self.nodes.index(node), neighbor.name, n2.name,
                   segment.name, n2.name)

    subGraphs = []
    checkObjs = {obj for obj in checkObjects}
    while checkObjs:
      # start checking new subgraph
      start = checkObjs.pop()
      connected = { (start, None) }
      subGraph = { start }
      pathSegNames = { start : [start.name] }
      paths = {start : []}
      
      while connected:
        # find all the elements connected to this subgraph
        segment, startNode = connected.pop()
        for neighbor, (pos, nPos, node) in zip(segment.neighbors,
                                               segment.neighborLocations):
          if node == startNode and neighbor != segment:
            # this is just backtracking
            continue
          if neighbor in checkObjs:
            connected.add((neighbor, node))
            subGraph.add(neighbor)
            checkObjs.remove(neighbor)
            pathSegNames[neighbor] = pathSegNames[segment] + [neighbor.name]
            paths[neighbor] = paths[segment] + [(segment, node, neighbor)]
          else:
            # there is a loop!
            names1 = pathSegNames[segment]
            try:
              names2 = pathSegNames[neighbor]
            except KeyError as err:
              print(segment.name, neighbor.name, neighbor in self.segments,
                    neighbor in checkObjs, neighbor in subGraph,
                    neighbor in pathSegNames)
              print([neighbor in s for s in subGraphs])
              raise err
                
            ind = 0
            for name1, name2 in zip(names1, names2):
              if name1 != name2:
                ind -= 1
                break
              ind += 1
            loopSegNames = names1[ind:] + names2[:ind-1:-1]
            if removeLoops:
              #loops.append(loopSegNames)
              warn('Have not implement loop removal.\nLoop detected',
                    '->'.join(loopSegNames))
            else:
              pass
              # warn('Loop detected', '->'.join(loopSegNames))
      
      subGraphs.append(subGraph)
    
    # sort the subgraphs so that the largest is first
    if isinstance(checkObjects[0], Segment):
      checkType = 'Segment'
      subGraphs.sort(key=lambda x: sum([len(y.compartments) for y in x]))
    elif isinstance(checkObjects[0], Compartment):
      checkType = 'Compartment'
      subGraphs.sort(key=lambda x: len(x))
    else:
      raise RuntimeError("Can't sort type: %s" % str(type(checkObjects[0])))

    if debugInfo:
      print('Number of subgraphs = %d / size of graphs: %s'
            % (len(subGraphs), str([len(graph) for graph in subGraphs])))
    
    if removeDisconnected and len(subGraphs) > 1:
      badGraphs, subGraphs = subGraphs[:-1], subGraphs[-1]

      badSegs = set()
      badComps = set()
      badNodes = set()
      if checkType == 'Segment':
        while badGraphs:
          subGraph = badGraphs.pop()
          
          # find unwanted objects in the subgraph
          badSegs.update(subGraph)
          
          for seg in subGraph:
            badComps.update(seg.compartments)
            badNodes.update(seg.nodes)
      else:
        while badGraphs:
          subGraph = badGraphs.pop()
          
          # find unwanted objects in the subgraph
          badComps.update(subGraph)
          
          for comp in subGraph:
            badNodes.update(comp.nodes)
            badSegs.add(comp.segment)

      self.segments[:] = \
        [seg for seg in self.segments if seg not in badSegs]
      self.compartments[:] = \
        [comp for comp in self.compartments if comp not in badComps]
      self.nodes[:] = [node for node in self.nodes if node not in badNodes]
      self.branches = []
      if self._somaBranch is not None:
        self._somaBranch[0].neighbors = []
      
      checkHash = hash(str(checkObjects)+str(removeDisconnected))
      
      print("Removed all but largest subgraphs")
    
    # record that the connectivity is already checked
    self._connectivityChecked.add(checkHash)
    
    return subGraphs


  def _plotShollGraph(self, distances):
    """
    Plot the number of neurites at a given distance
    """
    
    # neuriteDistance starts at zero, and has two data points for each
    # distance: one with the previous (running) number of compartments, and one
    # with the change added in (running +1 or -1)
    runningNum = 0
    neuriteDistance = [0.0]
    numIntersections = [0]
    lastNeuriteDistance = 0.0
    for d in distances:
      if d[0] > lastNeuriteDistance:
        neuriteDistance.append(lastNeuriteDistance)
        numIntersections.append(runningNum)
        neuriteDistance.append(d[0])
        numIntersections.append(runningNum)
        lastNeuriteDistance = d[0]
      runningNum += d[1]
    neuriteDistance.append(lastNeuriteDistance)
    numIntersections.append(runningNum)
    
    fig = pyplot.figure()
    pyplot.plot(neuriteDistance, numIntersections, 'k-')
    pyplot.title('Sholl Analysis', fontsize=22)
    pyplot.xlabel('Distance from soma', fontsize=22)
    pyplot.ylabel('Number of compartments', fontsize=22)
    ax = pyplot.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(16)
    pyplot.tight_layout()

  
  def shollAnalysis(self, straightenNeurites=True):
    """
    Find the number of neurites that intersect a sphere of a given radius
    """
    # get the centroid of the soma, weighting each compartment's contribution
    # by volume
    
    # define how distance from centroid to compartment is measured
    if straightenNeurites:
      # get distance traveled along neurites (e.g. as though neuron was
      # straightened out
      
      centroid = self.soma.centroidPosition(mandateTag='Soma')
      
      # compute distance from soma to each segment
      somaPaths = PathDistanceFinder(self, self.soma, centroid)
      
      # store results in an array
      distances = []
      for s in self.segments:
        d0, d1 = somaPaths.distanceTo(s, 0.0), somaPaths.distanceTo(s, 1.0)
        if d1 < d0:
          d0, d1 = d1, d0
        
        distances.append((d0, 1))
        distances.append((d1, -1))
      
    else:
      # get euclidean distance from soma centroid to each compartment
      # (must be done compartment by compartment, because segments curve)

      centroid = self.soma.centroid(mandateTag='Soma')
      
      # define distance from centroid to compartment
      def _centroidDist(c):
        def _tupleDist(_t1, _t2):
          return sqrt( (_t1[0] - _t2[0])**2 + \
                       (_t1[1] - _t2[1])**2 + \
                       (_t1[2] - _t2[2])**2 )

        d0 = _tupleDist(centroid, (c.x0, c.y0, c.z0))
        d1 = _tupleDist(centroid, (c.x1, c.y1, c.z1))
        if d0 <= d1:
          return d0, d1
        else:
          return d1, d0
      
      # compute the distance from the soma centroid to each compartment that's
      # not in the soma
      distances = []
      for c in self.compartments:
        if 'Soma' in c.tags:
          continue
        d0, d1 = _centroidDist(c)
        distances.append((d0, 1))  
        distances.append((d1, -1))
    
    # sort the distances by increasing distance
    distances.sort(key=lambda x: x[0])
    
    self._plotShollGraph(distances)
  
    
  def _addSegment(self, name, segList=None):
    """
    Add a new segment to the model
    """
    if name in self.tags:
      raise IOError('Tried to create segment with existing name/tag')

    newSeg = Segment(self)
    newSeg.name = name
    if name not in self.tags:
      self.tags[name] = 0
    
    if segList is None:
      segList = self.segments
    segList.append(newSeg)
    return newSeg
  
  
  def _addNode(self, segment, _x, _y, _z, _r1, _r2=None, _r3=None,
               _theta = 0.0, _phi=0.0):
    """
    Define and add a new node to the model within the specified segment
    """
    newNode = Node(_x, _y, _z, _r1, _r2, _r3, _theta, _phi)
    newNode.segments.append(segment)
    newNode.tags.update(segment.tags)
    newNode.tags.add(segment.name)
    self.nodes.append(newNode)
    segment.nodes.append(newNode)
    return newNode
  
  
  def _addCompartment(self, segment, node0, node1=None, append=False):
    """
    Define and add compartment to geometry within specified segment
    """
    if type(node0) is int:
      # passed an index, get the node object
      node0 = self.nodes[node0]

    if node1 is None:
      # define one-node compartment and add it to the segment
      newComp = OneNodeCompartment(node0)
      node0.compartments.append(newComp)
      if node0 == segment.nodes[0]:
        segment.compartments.insert(0, newComp)
      else:
        assert node0 == segment.nodes[-1]
        segment.compartments.append(newComp)
    else:
      # define two-node compartment and add it to the segment
      if type(node1) is int:
        node1 = self.nodes[node1]
      newComp = TwoNodeCompartment(node0, node1)
      node0.compartments.append(newComp)
      node1.compartments.append(newComp)
      if append or _node1 == segment.nodes[-1]:
        segment.compartments.append(newComp)
      else:
        compInd = segment.nodes.index(node0)
        if len(segment.compartments) and \
           isinstance(segment.compartments[0], OneNodeCompartment):
          compInd += 1
        segment.compartments.insert(compInd, newComp)
    
    # add segment information to compartment
    newComp.tags.update(segment.tags)
    newComp.tags.add(segment.name)
    newComp.segment = segment
    
    # add compartment to geometry
    self.compartments.append(newComp)
    
    # update tag counts
    self.tags['*'] += 1
    for tag in newComp.tags:
      self.tags[tag] += 1
    
    # update geometry totals
    self.surfaceArea += newComp.surfaceArea
    self.volume += newComp.volume
    return newComp
  
  
  def _connectSegments(self, segment0, location0, segment1, location1,
                       implicitConnect=True):
    """
    Connect two segments in the geometry, managing the connecting info in the
    segments and their nodes and compartments
    if implicitConnect is set to True:
      if A connects to B at position (x,y,z)
          AND _connectSegments is called to connect A and C and (x,y,z):
        then ensure that A, B, and C are all connected at (x,y,z)
    """
    node0 = segment0.nodeAt(location0)
    node1 = segment1.nodeAt(location1)
    # check to make sure the two nodes are identical
    try:
      assert (node0.x, node0.y, node0.z, node0.r1) == \
             (node1.x, node1.y, node1.z, node1.r1), \
        'Tried to connect %s (at %f) and %s (at %f) but connecting nodes ' \
          % (segment0.name, location0, segment1.name, location1 ) \
          + 'were different'
        
    except AssertionError:
      # buggy .hoc file with non-matching connecting nodes
      # check to see if there is a location where these segments COULD connect
      
      # find location where segment0 and segment1 connect
      numMatches = 0
      oldLoc0, oldLoc1 = location0, location1
      node0, location0, node1, location1 = None, None, None, None
      if not segment0.nodeLocations:
        segment0._setNodeLocations()
      if not segment1.nodeLocations:
        segment1._setNodeLocations()
      for n0, l0 in zip(segment0.nodes, segment0.nodeLocations):
        for n1, l1 in zip(segment1.nodes, segment1.nodeLocations):
          if (n0.x, n0.y, n0.z, n0.r1) == (n1.x, n1.y, n1.z, n1.r1):
            numMatches += 1
            node0, location0, node1, location1 = n0, l0, n1, l1
        
      if numMatches == 0:
        # no valid connection
        print(len(segment0.nodeLocations), len(segment1.nodeLocations))
        print('No way to connect %s and %s' % (segment0.name, segment1.name))
        raise
      elif numMatches > 1:
        # ambiguous as to where to connect
        print('Multiple possible connections between %s and %s, giving up'
              % (segment0.name, segment1.name))
      #else:
        #warn('incorrect connection location for %s (at %.1f) and %s (at %.1f)'
        #     % (segment0.name, oldLoc0, segment1.name, oldLoc1),
        #     'actual connection at %.1f and %.1f' % (location0, location1))
    
    
    def _replaceNode(_oldNode, _replacementNode):
      """
      swap _oldNode with _replacementNode in all objects that use _oldNode
      """
      # swap nodes in all relevant segments
      #_replacementNode.segments.extend(_oldNode.segments)
      for _seg in _oldNode.segments:
        _ind = _seg.nodes.index(_oldNode)
        _seg.nodes[_ind] = _replacementNode
        if _seg not in _replacementNode.segments:
          _replacementNode.segments.append(_seg)
        for _nInd in range(len(_seg.neighbors)):
          _loc, _nLoc, _node = _seg.neighborLocations[_nInd]
          if _node == _oldNode:
            _seg.neighborLocations[_nInd] = (_loc, _nLoc, _replacementNode)
      
      # swap nodes in all relevant compartments
      #_replacementNode.compartments.extend(_oldNode.compartments)
      for _comp in _oldNode.compartments:
        _ind = _comp.nodes.index(_oldNode)
        _comp.nodes[_ind] = _replacementNode
        if _comp not in _replacementNode.compartments:
          _replacementNode.compartments.append(_comp)
      
      # remove node1 from geometry
      self._removeNodes.add(_oldNode)

    
    if implicitConnect:
      # loop through neighbors of segment0 at location0
      for nSeg, (loc, nLoc, node) in zip(segment0.neighbors[:],
                                         segment0.neighborLocations[:]):
        if node == node0:
          assert loc == location0
        if loc == location0:
          if nSeg == segment1:
            assert nLoc != location1, \
              'Tried to connect %s and %s, which are already connected' \
              % (segment0.name, segment1.name)
            
            # segment0 (at location0) connects to segment1 in two locations!
            #  location1 (the new connection), and
            #  nLoc (an older connection)
            assert node == node0, 'Weird node mismatch near loop segment'
            # this is unneccessary if the assertion is true (it should be):
            # _replaceNode(node, node0)
            
          # make nSeg neighbors with segment1
          _makeNeighbors(segment1, nSeg, location1, nLoc, node0,
                         checkDuplicate=True)
          # make nSeg neighbors with all of segment1's neighbors at loc1
          for nSegN, (locN, nLocN, nodeN) in zip(segment1.neighbors[:],
                                                segment1.neighborLocations[:]):
            if nodeN != node1:
              continue
            assert locN == location1
            _makeNeighbors(nSeg, nSegN, loc, nLocN, node0)
      
      for nSeg, (loc, nLoc, node) in zip(segment1.neighbors[:],
                                         segment1.neighborLocations[:]):
        if node == node1:
          assert loc == location1
        if loc == location1:
          if nSeg == segment0:
            assert nLoc != location0, \
              'Tried to connect %s and %s, which are already connected' \
              % (segment0.name, segment1.name)
            # segment1 (at location1) connects to segment0 in two locations!
            #  location0 (the new connection), and
            #  nLoc (an older connection)
            assert node == node1, 'Weird node mismatch near loop segment'
            # this is unneccessary if the assertion is true (it should be):
            #_replaceNode(node, node0)

          # make nSeg neighbors with segment0
          _makeNeighbors(segment0, nSeg, location0, nLoc, node0,
                         checkDuplicate=True)
    
    # make segment0 and segment1 neighbors, with node0 as the connecting node
    _makeNeighbors(segment0, segment1, location0, location1, node0)

    # remove node1 from the geometry and replace it everywhere with node0
    _replaceNode(node1, node0)

  
  def _mergeSegments(self, segmentA, segmentB, _segList):
    # Merge segmentA and segmentB into one segment, preserving their neighbor
    # information
    #
    assert segmentA in _segList
    assert segmentB in _segList
    assert segmentB in segmentA.neighbors and \
           segmentA in segmentB.neighbors, \
          "Merged segments must be neighbors"
    
    # find the connection location
    ind = segmentA.neighbors.index(segmentB)
    (locationA, locationB, nodeAB) = segmentA.neighborLocations[ind]
    assert locationA in [0.0, 1.0] and locationB in [0.0, 1.0], \
           "Tried to merge with locationA=%g, locationB=%g, but both locations must be an end point" \
           % (locationA, locationB)
    
    lenA = segmentA.length
    lenB = segmentB.length
    newLen = lenA + lenB
    
    if locationA == 1.0:
      if locationB == 0.0:
        # define function to get location on the new segment
        # add segmentB compartments + nodes to segmentA
        segmentA.compartments += segmentB.compartments
        segmentA.nodes += segmentB.nodes[1:]
      else: #locationB == 1.0:
        # add segmentB compartments + nodes to segmentA
        segmentA.compartments += reversed(segmentB.compartments)
        segmentA.nodes += reversed(segmentB.nodes[:-1])
    else:
      if locationB == 0.0:
        # add segmentB compartments + nodes to segmentA
        segmentA.compartments = list(reversed(segmentB.compartments)) \
                              + segmentA.compartments
        segmentA.nodes = list(reversed(segmentB.nodes[1:])) \
                       + segmentA.nodes
      else: #locationB == 1.0:
        # add segmentB compartments + nodes to segmentA
        segmentA.compartments = segmentB.compartments + segmentA.compartments
        segmentA.nodes = segmentB.nodes[:-1] + segmentA.nodes

    # update node locations
    segmentA._setNodeLocations()
    
    # remove segmentB from segmentA's neighbors
    _removeNeighbor(segmentA, segmentB)
    # remove segmentA from segmentB's neighbors
    _removeNeighbor(segmentB, segmentA)
    
    # recompute the location of segmentA's old neighbors
    for ind in range(len(segmentA.neighbors)):
      neighbor = segmentA.neighbors[ind]
      loc, nLoc, node = segmentA.neighborLocations[ind]
      nodeInd = segmentA.nodes.index(node)
      newLoc = segmentA.nodeLocations[nodeInd]

      # update the connection location in segmentA
      segmentA.neighborLocations[ind] = (newLoc, nLoc, node)
      # update the connection location in the neighbor
      nInd = zip(neighbor.neighbors,
                 neighbor.neighborLocations).index((segmentA,
                                                    (nLoc, loc, node)))
      neighbor.neighborLocations[nInd] = ((nLoc, newLoc, node))
    
    # replace segmentB with segmentA in other segments' neighbors
    #  and add those neighbors as new neighbors to segmentA
    for neighbor, (loc, nLoc, node) in zip(segmentB.neighbors,
                                           segmentB.neighborLocations):
      # update the location
      nInd = zip(neighbor.neighbors,
                 neighbor.neighborLocations).index((segmentB,
                                                    (nLoc, loc, node)))
      nodeInd = segmentA.nodes.index(node)
      newLoc = segmentA.nodeLocations[nodeInd]

      # make segmentB no longer neighbor's neighbor
      #_removeNeighbor(neighbor, segmentB)
      # make segmentA and neighbor neighbors
      #_makeNeighbors(segmentA, neighbor, newLoc, nLoc, node)
      
      # replace segmentB by segmentA in neighbor.neighbors, at new location
      neighbor.neighbors[nInd] = segmentA
      neighbor.neighborLocations[nInd] = (nLoc, newLoc, node)
      
      # add neighbor to segmentA's list of neighbors
      segmentA.neighbors.append(neighbor)
      segmentA.neighborLocations.append((newLoc, nLoc, node))
    
    # update segmentB's tags into segmentA
    segmentA.tags.update(segmentB.tags)
    segmentA.tags.add(segmentB.name)
    
    # remove segmentB from the list of segments
    #[SLOW]
    #_segList.remove(segmentB)
    self._removeSegments.add(segmentB)
  
  
  def mergeBranchesByDistanceToEdge(self, makePlots=True):
    if not hasattr(self.soma, 'centripetalOrder') or \
        self.soma.centripetalOrder is None:
      self.calcBranchOrder(doPlot=False)
    
    def _getMergePath(startB, considerBranches):
      startLen = startB.length
      if startB.branchOrder == 0:
        startLen *= 0.5
      finishedPaths = []
      # start at startB
      openPaths = [([startB], startLen)]
      # find all the paths from startB to a neurite tip
      while openPaths:
        path, length = openPaths.pop()
        current = path[-1] ; nextOrder = current.branchOrder + 1
        
        nextBranches = [n for n in current.neighbors
                        if n.branchOrder == nextOrder]
        if nextBranches:
          for n in nextBranches:
            openPaths.append((path + [n], length + n.length))
        else:
          finishedPaths.append((path, length))
      
      # return the longest path
      return max(finishedPaths, key=lambda x:x[1])[0]
    
    self.mergedBranches = []
    considerBranches = {b for b in self.branches}
    
    while considerBranches:
      minOrder = min(b.branchOrder for b in considerBranches)
      startB = [b for b in considerBranches
                if b.branchOrder == minOrder].pop()
      mergePath = _getMergePath(startB, considerBranches)
      
      merged = Segment(self)
      merged.name = 'MergedBranch%d' % len(self.mergedBranches)
      # add the compartments, nodes, tags, and branch order
      merged.compartments = startB.compartments[:]
      merged.nodes = startB.nodes[:]
      merged.tags = {t for t in startB.tags}
      merged.branchOrder = startB.branchOrder
      merged.centripetalOrder = startB.centripetalOrder
      # add .merged element to startB
      startB.merged = merged
      
      considerBranches.remove(startB)
      previous = startB
      
      for b in mergePath[1:]:
        # add the compartments, nodes, and tags
        merged.compartments += [TwoNodeCompartment(previous.nodes[-1],
                                                   b.nodes[0])]
        merged.compartments += b.compartments
        merged.nodes += b.nodes
        merged.tags.update(b.tags)
        # add .merged element to b
        b.merged = merged
        
        considerBranches.remove(b)
        previous = b
      
      # set nodeLocations in merged
      merged._setNodeLocations()
      
      # add merged to mergedBranches
      self.mergedBranches.append(merged)
      
    # find the neighbors and neighborLocations of the merged branches
    for branch in self.branches:
      merged = branch.merged
      for n, (loc, nLoc, node) in zip(branch.neighbors,
                                      branch.neighborLocations):
        if n == branch:
          # this is a loop
          assert (loc == 0 and nLoc == 1) or (loc == 1 and nLoc == 0)
          merged.neighbors.append(merged)
          merged.neighborLocations.append((loc, nLoc, node))
        else:
          nMerged = n.merged
          if nMerged == merged:
            # both branches are part of merged, so no neighbor need be added
            continue
          # find the location of the neighbors on each merged branch
          mergedLoc = merged.nodeLocations[merged.nodes.index(node)]
          nMergedLoc = nMerged.nodeLocations[nMerged.nodes.index(node)]
          # make them neighbors
          _makeNeighbors(merged, nMerged, mergedLoc, nMergedLoc, node)
    
    if makePlots:
      self._plotBranchOrderStatistics()

  
  def mergeBranchesByOrder(self, makePlots=True):
    # merge branches together if one of a branch's neighbors is clearly a
    # continuation of it (using centripetalOrder as method of determining
    # continuation)
    
    if not hasattr(self.soma, 'centripetalOrder') or \
        self.soma.centripetalOrder is None:
      self.calcBranchOrder(doPlot=False)
    
    def _getMergeBranch(currentBranch):
      # find the branch that should be merged with current branch
      # and what the order of the merged branch should be
      # if no branch should be merged, return None
      orders0 = {}
      orders1 = {}
      
      mergeBranch = None
      mergeCentripetalOrder = currentBranch.centripetalOrder + 1
      mergeBranchOrder = currentBranch.branchOrder - 1
      for n, (locC, locN, nodeC) in zip(currentBranch.neighbors,
                                        currentBranch.neighborLocations):
        assert locC in [0, 1] and locN in [0, 1]
        #if locC not in [0, 1] or locN not in [0, 1]:
        #  # only merge at branch end points
        #  continue
        
        
        if n.centripetalOrder == mergeCentripetalOrder and \
           n.branchOrder == mergeBranchOrder:
          if mergeBranch is not None:
            # too many potential merges => this branch can't merge
            return None
          else:
            # potential merge
            mergeBranch = n
      
      if mergeBranch is not None:
        # check to make sure no other branch can merge onto mergeBranch:
        for n in mergeBranch.neighbors:
          if n != currentBranch and n.branchOrder == currentBranch.branchOrder\
             and n.centripetalOrder == currentBranch.centripetalOrder:
            # another branch could equally merge onto mergeBranch, which means
            # that no merge can be preferred
            return None
      
      return mergeBranch
    
    self.mergedBranches = []
    visited = {b : False for b in self.branches}
    openBranches = {b for b in self.branches if b.centripetalOrder == 0}
    
    # find the new merged branches
    while openBranches:
      # get the next branch that could start a merged branch
      current = openBranches.pop()
      visited[current] = True

      # mark all the neighbors as open branches
      #openBranches.update(n for n in current.neighbors if not visited[n])

      # create the merged branch that will hold current and anything that
      # continues it
      merged = Segment(self)
      merged.name = 'MergedBranch%d' % len(self.mergedBranches)
      # add the compartments, nodes, tags, and branch order
      merged.compartments = current.compartments[:]
      merged.nodes = current.nodes[:]
      merged.tags = {t for t in current.tags}
      merged.branchOrder = current.branchOrder
      merged.centripetalOrder = current.centripetalOrder
      # add .merged element to current
      current.merged = merged
      barf = (current.name, current.branchOrder, current.centripetalOrder)
      # find the continuation, if any
      current = _getMergeBranch(current)
      while current is not None:
        # don't allow any other branches to use current branch
        openBranches.discard(current)
        if visited[current]:
          print('%s with tags: %s' % (current.name,
                                      ' '.join(t for t in current.tags)))
          print(current.branchOrder, current.centripetalOrder)
          for n in current.neighbors:
            print('%s with tags: %s' % (n.name,
                                      ' '.join(t for t in n.tags)))
            print(n.branchOrder, n.centripetalOrder)
          print(barf)
          if current.merged is not None:
            print(current.merged.name)
            for b in self.branches:
              if hasattr(b, 'merged') and b.merged == current.merged:
                print(b.name)
          else:
            print('None')
          
        visited[current] = True
        
        # add the compartments, nodes, and tags
        merged.compartments[0:0] = [TwoNodeCompartment(current.nodes[-1],
                                    merged.nodes[0])]
        merged.compartments[0:0] = current.compartments
        merged.nodes[0:0] = current.nodes
        
        merged.tags.update(current.tags)
        merged.branchOrder = current.branchOrder
        merged.centripetalOrder = current.centripetalOrder
        
        # add .merged element to current
        current.merged = merged
        
        # find the continuation, if any
        barf = (current.name, current.branchOrder, current.centripetalOrder)
        current = _getMergeBranch(current)
      
      # set nodeLocations in merged
      merged._setNodeLocations()
      # add merged to mergedBranches
      self.mergedBranches.append(merged)
      # add any missed branches
      if not openBranches:
        try:
          minUnvisitedOrder = min(b.centripetalOrder for b in self.branches
                                  if not visited[b])
          openBranches = {b for b in self.branches if not visited[b] and
                          b.centripetalOrder == minUnvisitedOrder}
        except ValueError:
          openBranches = {}

    # find the neighbors and neighborLocations of the merged branches
    for branch in self.branches:
      merged = branch.merged
      for n, (loc, nLoc, node) in zip(branch.neighbors,
                                      branch.neighborLocations):
        if n == branch:
          # this is a loop
          assert (loc == 0 and nLoc == 1) or (loc == 1 and nLoc == 0)
          merged.neighbors.append(merged)
          merged.neighborLocations.append((loc, nLoc, node))
        else:
          nMerged = n.merged
          if nMerged == merged:
            # both branches are part of merged, so no neighbor need be added
            continue
          # find the location of the neighbors on each merged branch
          mergedLoc = merged.nodeLocations[merged.nodes.index(node)]
          nMergedLoc = nMerged.nodeLocations[nMerged.nodes.index(node)]
          # make them neighbors
          _makeNeighbors(merged, nMerged, mergedLoc, nMergedLoc, node)
    
    if makePlots:
      self._plotBranchOrderStatistics()

  
  def _addOneNodeCompartments(self):
    """
    Add extra compartments to account for regions covered by nodes with only
    one compartment attached (the end of the soma is probably the most
    important)
    """
    for segment in self.segments:
      if not segment.nodes:
        # empty segment
        warn("Warning, empty segment: %s" % segment.name)
      
      node0 = segment.nodes[0]
      if len(node0.compartments) == 1:
        self._addCompartment(segment, node0)
      
      node1 = segment.nodes[-1]
      if len(node1.compartments) == 1:
        self._addCompartment(segment, node1)
  
  
  def _getBranch(self, branchStart):
    """
    self._getBranch(branchStart)
      branchStart: starting compartment
    returns branch, neighbors0, neighbors1
      branch: a Segment with 0 or >= 2 neighbors at each endpoint, and no
        neighbors in the middle
      neighbors0, neighbors1: list of neighbors at 0 and 1 end respectively
        each neighbor is a tuple (segment, position, compartment)
          segment: the neighbor segment
          position: the position (in the neighbor segment) of the connection
          compartment: the neighbor compartment
    Find the branch (bounded by compartments with more than two neighbors)
    that starts with startCompartment. The end of each branch is a node with
    a number of attached compartments not equal to 2.
    NOTE: must add compartment-by-compartment (instead of by segments) because
      segments can in principle have neighbors midway instead of at endpoints
    NOTE: don't add _neighbors0 or _neighbors1 to _branch.neighbors, because
      they are *Compartment* neighbors, and _branch will want *branch*
      neighbors
    
    """
    
    assert isinstance(branchStart, Compartment), \
      'Must start a branch with a compartment'
    # determine starting segment, position of starting compartment in segment
    startSeg = branchStart.segment
    startPos = branchStart.length / 2.0
    for c in startSeg.compartments:
      if c == branchStart:
        break
      startPos += c.length
      
    # create the branch
    branch = Segment(self)
    branch.name = 'branch%d' % len(self.branches)
    branch.tags = {t for t in startSeg.tags}
    branch.tags.add(startSeg.name)
    
    def _getBranchNeighbors(segment, startPos):
      neighbors0, neighbors1 = [], []
      pos0, nPos0, pos1, nPos1 = 0, None, 1, None
      for neighbor, (location, nLocation, node) in zip(segment.neighbors,
                                                    segment.neighborLocations):
        nComps = [c for c in node.compartments if c in neighbor.compartments]
        # this assertion fails when a segment is its own neighbor:
        #assert len(nComps) == (2 - (nLocation in [0, 1])), \
        #  'Found %d neighbor compartments at %f' % (len(nComps), nLocation)
        if location == 0:
          neighbors0.append((neighbor, nLocation, nComps[0]))
          if nLocation not in [0, 1]:
            neighbors0.append((neighbor, nLocation, nComps[-1]))
        elif location == 1:
          neighbors1.append((neighbor, nLocation, nComps[0]))
          if nLocation not in [0, 1]:
            neighbors1.append((neighbor, nLocation, nComps[-1]))
        elif location < startPos:
          assert nLocation in [0, 1]
          if location > pos0:
            pos0 = location
            farComp = segment.compartmentAt(location, choice=0)
            neighbors0 = [(neighbor, nLocation, nComps[0]),
                          (segment, location, farComp)]
          elif location == pos0:
            farComp = segment.compartmentAt(location, choice=0)
            neighbors0 += [(neighbor, nLocation, nComps[0]),
                          (segment, location, farComp)]
        else:
          assert nLocation in [0, 1]
          if location < pos1:
            pos1 = location
            farComp = segment.compartmentAt(location, choice=1)
            neighbors1 = [(neighbor, nLocation, nComps[0]),
                          (segment, location, farComp)]
          else:
            farComp = segment.compartmentAt(location, choice=1)
            neighbors1 += [(neighbor, nLocation, nComps[0]),
                          (segment, location, farComp)]
      
      return neighbors0, pos0, neighbors1, pos1
      
    
    def _getBranchPart(segment, pos0, pos1):
      # return the compartments and nodes in segment between pos0 and pos1
      if pos0 == 0:
        if pos1 == 1:
          return segment.compartments[:], segment.nodes[:]
        else:
          pos0L = 0.0
      else:
        pos0L = pos0 * segment.length
      
      if not segment.nodeLocations:
        locs = list(cumsum(comp.length for comp in segment.compartments))
        segment.nodeLocations = [loc / locs[-1] for loc in locs]
      n0 = bisect_left(segment.nodeLocations, pos0)
      assert segment.nodeLocations[n0] == pos0
      n1 = bisect_left(segment.nodeLocations, pos1)
      assert segment.nodeLocations[n1] == pos1
      n1 -= 1
      
      return segment.compartments[n1:n2], segment.nodes[n1:n2+1]
    
    # find neighbors of startSeg
    neighbors0, pos0, neighbors1, pos1 = \
      _getBranchNeighbors(startSeg, startPos)
    # get the part of startSeg between pos0 and pos1
    branch.compartments, branch.nodes = _getBranchPart(startSeg, pos0, pos1)
    
    while len(neighbors0) == 1:
      # extend branch in 0 direction
      
      # get the neighbor segment
      segment, nPos0, nComp = neighbors0.pop()
      # add the segment's name to branch tags
      branch.tags.add(segment.name)
      if nPos0 == 1:
        # normal case, the neighbor is oriented in the same direction
        
        # get the new neighbors0
        neighbors0, pos0, nDummy, posDummy = _getBranchNeighbors(segment, 1)
        # update branch with neighbor's segment, removing duplicate node
        branch.compartments[0:0], branch.nodes[0:1] = \
          _getBranchPart(segment, pos0, 1)
      else:
        # the neighbor is oriented backwards relative to startSeg
        
        # get the new neighbors0
        nDummy, posDummy, neighbors0, pos0 = _getBranchNeighbors(segment, 0)
        # update branch with neighbor's segment, removing duplicate node
        branch.compartments[0:0], branch.nodes[0:1] = (
          reversed(L) for L in _getBranchPart(segment, 0, pos0))
    
    while len(neighbors1) == 1:
      # extend branch in 1 direction

      # get the neighbor segment
      segment, nPos1, nComp = neighbors1.pop()
      # add the segment's name to branch tags
      branch.tags.add(segment.name)
      if nPos1 == 0:
        # normal case, the neighbor is oriented in the same direction
        
        # get the new neighbors1
        nDummy, posDummy, neighbors1, pos1 = _getBranchNeighbors(segment, 0)
        # update branch with neighbor's segment, removing duplicate node
        numC, numN = len(branch.compartments), len(branch.nodes)
        branch.compartments[numC:], branch.nodes[numN-1:] = \
          _getBranchPart(segment, 0, pos1)
      else:
        # the neighbor is oriented backwards relative to startSeg
        
        # get the new neighbors1
        neighbors1, pos1, nDummy, posDummy = _getBranchNeighbors(segment, 1)
        # update branch with neighbor's segment, removing duplicate node
        numC, numN = len(branch.compartments), len(branch.nodes)
        branch.compartments[numC:], branch.nodes[numN-1:] = (
          reversed(L) for L in _getBranchPart(segment, pos1, 1))
  
    return branch, neighbors0, neighbors1

  
  @property
  def soma(self):
    if self._soma is None:
      self._findSoma()
    return self._soma
  
  @property
  def somaBranch(self):
    if self._somaBranch is None:
      self._findSoma()
    return self._somaBranch[0]
  
  
  def _findSoma(self):
    """
    Find the compartment with the largest diameter (which is presumably in the
    soma). Return the segment or branch that contains this compartment.
    
    """
    # locate the largest radius compartment -it should be the soma
    somaCenter = max(self.compartments, key=lambda c: c.maxRadius)
    
    # get the segment that contains the soma (plus extra stuff)
    soma = somaCenter.segment
    self._soma = soma
    
    # add the Soma tag to the segment, but not its compartments
    soma.tags.add('Soma')
    
      # find a cutoff on radius that defines the Soma proper
    rCutoff = (soma.minRadius * soma.maxRadius**3)**0.25
    
    # find the index to the center of the Soma
    centerInd = soma.compartments.index(somaCenter)
    
    # apply the Soma tag to contiguous stretch of large compartments near
    # the center
    self.tags['Soma'] = 0
    for ind in range(centerInd, len(soma.compartments)):
      c = soma.compartments[ind]
      if c.avgRadius >= rCutoff:
        c.tags.add('Soma')
        self.tags['Soma'] += 1
      else:
        break
    for ind in range(centerInd - 1, -1, -1):
      c = soma.compartments[ind]
      if c.avgRadius >= rCutoff:
        c.tags.add('Soma')
        self.tags['Soma'] += 1
      else:
        break

    self._somaBranch = self._getBranch(somaCenter)
  
  
  def findAxons(self, findBranch=False, debugInfo=False, minLength=10,
                edgeSafety=2.0):
    """
    Locate any axon branches
    To be an axon branch, it must:
      1. be terminal, but not contain the Soma
      2. be longer then minLength
      3. have its terminal node less than edgeSafety * node radius from the
         furthest extent of the neuron in x, y, or z
      4. AND, for pyramidal, be the longest axon
    """
    
    for segment in self.segments:
      if len(segment.nodes) != len(segment.compartments) + 1:
        raise AssertionError('%s has %d nodes and %d compartments'
          % (segment.name, len(segment.nodes), len(segment.compartments)))
    
    if self._axons is not None:
      # already found the axons, return the results
      if findBranch:
        return self._axonBranches
      else:
        return self._axons

    # need the Soma to be tagged
    somaBranch = self.somaBranch
    # need to know the range of the whole geometry
    self._findGeometryRange()
    # need the branches, because one criterion for axons is length, and
    # one axon may be made up of many consecutive segments
    self.findBranches()

    def _isEdgeNode(_node, _safety):
      _radius = _node.avgRadius
      def _edgeCoord(_coord, _dim):
        return _coord - _safety * _radius < self.minRange[_dim] or \
               _coord + _safety * _radius > self.maxRange[_dim]
      return _edgeCoord(_node.x, 0) or \
             _edgeCoord(_node.y, 1) or \
             _edgeCoord(_node.z, 2)

    if debugInfo:
      print('Edges: %6.1f %6.1f %6.1f' % tuple(self.minRange))
      print('       %6.1f %6.1f %6.1f' % tuple(self.maxRange))

    self._axons = list()
    self._axonBranches = list()
    for branch in self.branches:
      if branch == somaBranch:
        # it's the Soma
        continue
      if branch.length < minLength:
        # too short to be sure it's an axon
        continue
      
      n0 = branch.neighborsAt(0)
      n1 = branch.neighborsAt(1)
      if n0 and n1:
        # not terminal
        continue

      if debugInfo:
        print('Possible axon: %s (%s)' % (branch.name,
                                        ' '.join(t for t in branch.tags)))
        if not n0:
          print('\t%d neighbors at 0.0' % len(n0))
          node = branch.nodeAt(0)
          print('\t node at 0.0: %.1f, %.1f, %.1f, %f' %
            (node.x, node.y, node.z, node.r1))
        if not n1:
          print('\t%d neighbors at 1.0' % len(n1))
          node = branch.nodeAt(1)
          print('\t node at 1.0: %.1f, %.1f, %.1f, %f' %
            (node.x, node.y, node.z, node.r1))
      
      if ((not n0 and _isEdgeNode(branch.nodes[0], edgeSafety)) or 
          (not n1 and _isEdgeNode(branch.nodes[-1], edgeSafety))):
        
        ## This way only the longest axon is kept
        if len(self._axonBranches) > 0:
          if branch.length > self._axonBranches[0].length:
            self._axonBranches.pop(0)
          #  branch.addTag('Axon')
            self._axonBranches.append(branch)
        else: # if there aren't any axons yet, just keep it
          self._axonBranches.append(branch)
        
        # after that, find the last segment in the only axon
        for b in self._axonBranches:
          b.addTag('Axon')
        print('self._axonBranches is length: %i' %len(self._axonBranches))
        segments = {c.segment for c in self._axonBranches[0].compartments}
        if debugInfo:
          print('Found axon with segments %s' %
                ' '.join(s.name for s in segments))
        for s in segments:
          s.addTag('Axon')
          if s.isTerminal:
            if len(self._axons) < 1:
              self._axons.append(s)
    
    if findBranch:
      return self._axonBranches
    else:
      return self._axons
  
  
  def _findGeometryRange(self):
    """
    Find the physical extent of the geometry in x,y,z
    """
    if (not any(isnan(x) for x in self.minRange) and
        not any(isnan(x) for x in self.maxRange)):
      # the range has already been found/specified
      return
    
    self.minRange = [float('inf'), float('inf'), float('inf')]
    self.maxRange = [float('-inf'), float('-inf'), float('-inf')]
    def _updateRange(_nodeCoord, _rangeInd):
      if _nodeCoord < self.minRange[_rangeInd]:
        self.minRange[_rangeInd] = _nodeCoord
      if _nodeCoord > self.maxRange[_rangeInd]:
        self.maxRange[_rangeInd] = _nodeCoord
    for _node in self.nodes:
      _updateRange(_node.x, 0)
      _updateRange(_node.y, 1)
      _updateRange(_node.z, 2)


def _makeNeighbors(segment1, segment2, location1, location2, node,
                   checkDuplicate=False):
  # make segment1 and segment2 neighbors
  if checkDuplicate:
    preExist = any(neighbor == segment2
                   and (loc, nLoc, n) == (location1, location2, node)
                   for neighbor, (loc, nLoc, n) in zip(segment1.neighbors,
                                                   segment1.neighborLocations))
    if not preExist:
      segment1.neighbors.append(segment2)
      segment1.neighborLocations.append((location1, location2, node))

    preExist = any(neighbor == segment1
                   and (loc, nLoc, n) == (location2, location1, node)
                   for neighbor, (loc, nLoc, n) in zip(segment2.neighbors,
                                                   segment2.neighborLocations))
    if not preExist:
      segment2.neighbors.append(segment1)
      segment2.neighborLocations.append((location2, location1, node))

  else:
    segment1.neighbors.append(segment2)
    segment1.neighborLocations.append((location1, location2, node))
    segment2.neighbors.append(segment1)
    segment2.neighborLocations.append((location2, location1, node))


def _removeNeighbor(segment, neighbor):
  # remove neighbor from list of segment's neighbors
  ind = segment.neighbors.index(neighbor)
  segment.neighbors.pop(ind)
  segment.neighborLocations.pop(ind)


def getBranchAngle(segment, neighbor, segLoc, nLoc, node):
  # calculate angle between segment and its neighbor
  
  # 2-node segment progression (produced from xml knossos skeletons)
  if len(segment.nodes) == 2 and len(neighbor.nodes) == 2:
    if segLoc == 0:
      segNode = segment.nodes[0]
    elif segLoc == 1:
      segNode = segment.nodes[-1]
    else:
      ind = segment.nodes.index(node)
      segNode = segment.nodes[ind] # why -1?
    if nLoc == 0:
      nNode = neighbor.nodes[0]
    elif nLoc == 1:
      nNode = neighbor.nodes[1]
    else:
      ind = neighbor.nodes.index(node)
      nNode = neighbor.nodes[ind]

  
  # normal progression
  else:
    if segLoc == 0:
      segNode = segment.nodes[1]
    elif segLoc == 1:
      segNode = segment.nodes[-2]
    else:
      ind = segment.nodes.index(node)
      segNode = segment.nodes[ind-1]
    
    if nLoc == 0:
      nNode = neighbor.nodes[1]
    elif nLoc == 1:
      nNode = neighbor.nodes[-2]
    else:
      ind = neighbor.nodes.index(node)
      nNode = neighbor.nodes[ind+1]
  
  segVec = (node.x - segNode.x, node.y - segNode.y, node.z - segNode.z)
  nVec = (nNode.x - node.x, nNode.y - node.y, nNode.z - node.z)
  dot = sum(sC * nC for sC, nC in zip(segVec, nVec))
  segMag = sum(sC * sC for sC in segVec)
  nMag = sum(nC * nC for nC in nVec)
  if len(segment.nodes) > 2:
    #cosAngle = max(-1.0, min(1.0, dot / sqrt(segment.length)))
    cosAngle = max(-1.0, min(1.0, dot / sqrt(segMag * nMag))) ########## MAJOR EDIT
  else:
    print('2-node segment')
    cosAngle = max(-1.0, min(1.0, dot / sqrt(segment.length)))
  try:
    angle = (180/pi) * acos(cosAngle)
  except ValueError:
    print(dot, segMag, nMag)
    print(dot / sqrt(segMag * nMag))
    raise
  return angle
  

class Segment:
  def __init__(self, geometry):
    self.geometry = geometry
    
    self.name = None
    self.tags = set()
    
    self.compartments = []
    self.nodes = []
    self.neighbors = []
    self.nodeLocations = []
    # neighborLocations are tuples:
    #   (location in this segment [0.0 - 1.0],
    #    location in neighbor segment [0.0 - 1.0],
    #    connecting Node)
    self.neighborLocations = []
    #self.volume = None
    self.branchOrder = None


  
  def neighborsAt(self, location):
    return [n for n, (loc, nLoc, node)
            in zip(self.neighbors, self.neighborLocations) if loc == location]
  
  
  def _setNodeLocations(self):
    try:
      locs = list(cumsum(comp.length for comp in self.compartments))
      self.nodeLocations = [loc / locs[-1] for loc in locs]
    except:
      self.nodeLocations = [0.0, 1.0]


  def compartmentAt(self, location, choice=None):
    if not self.nodeLocations:
      self._setNodeLocations()
    ind = bisect_left(self.nodeLocations, location)
    if self.nodeLocations[ind] == location:
      if choice == 0:
        return self.compartments[ind - 1]
      elif choice == 1:
        return self.compartments[ind]
      else:
        assert self.nodeLocations[ind] != location, \
          'Two compartments touch requested location'
    else:
      return self.compartments[ind - 1]


  def nodeAt(self, location):
    if not self.nodeLocations:
      self._setNodeLocations()
    ind = bisect_left(self.nodeLocations, location)
    assert self.nodeLocations[ind] == location, 'No node at requested location'
    return self.nodes[ind]


  def coordAt(self, location):
    if not self.nodeLocations:
      self._setNodeLocations()
    ind = bisect_left(self.nodeLocations, location)
    if self.nodeLocations[ind] == location:
      # exactly at a node:
      n = self.nodes[ind]
      return (n.x, n.y, n.z)
    else:
      # need to interpolate
      i0 = ind - 1
      n0 = self.nodes[i0]
      n1 = self.nodes[ind]
      cN0 = ( (self.nodeLocations[ind] - location) /
              (self.nodeLocations[ind] - self.nodeLocations[i0]) )
      cN1 = 1.0 - cN0
      return (cN0 * n0.x + cN1 * n1.x, cN0 * n0.y + cN1 * n1.y,
              cN0 * n0.z + cN1 * n1.z)

  def clear(self):
    if self.compartments:
      for c in self.compartments:
        for tag in c.tags:
          self.geometry.tags[tag] -= 1
      self.geometry.compartments = [c for c in self.geometry.compartments
                                    if c not in self.compartments]
      self.compartments = []
    if self.nodes:
      delNodes = []
      for n in self.nodes:
        n.segments.remove(self)
        if not n.segments:
          delNodes.append(n)
      self.geometry.nodes = [n for n in self.geometry.nodes
                             if n not in delNodes]
      self.nodes = []

  def addTag(self, newTag):
    """
    Add a new tag to the segment and all its compartments
    """
    if newTag not in self.geometry.tags:
      self.geometry.tags[newTag] = 0
    
    self.tags.add(newTag)
    for c in self.compartments:
      c.tags.add(newTag)
      self.geometry.tags[newTag] += 1
  
  @property
  def length(self):
    return sum([c.length for c in self.compartments])
  
  @property
  def surfaceArea(self):
    return sum([c.surfaceArea for c in self.compartments])
  
  @property
  def maxRadius(self):
    # compute maximum radius
    return max(c.maxRadius for c in self.compartments)
  
  @property
  def minRadius(self):
    # compute minimum radius
    return min(c.minRadius for c in self.compartments)
  
  @property
  def avgRadius(self):
    # compute average radius, weighted by volume
    return sum(c.avgRadius * c.volume for c in self.compartments) / \
           sum(c.volume for c in self.compartments)
  
  @property
  def volume(self):
    return sum(c.volume for c in self.compartments)
  
  @property
  def tortuosity(self):
    n0 = self.nodes[0] ; n1 = self.nodes[-1]
    if n0 == n1:
      return float('inf')
    euclideanD = sqrt((n0.x - n1.x)**2 + (n0.y - n1.y)**2 + (n0.z - n1.z)**2)
    if euclideanD == 0:
      return 1
    else:
      return self.length / euclideanD
    
  @property
  def isTerminal(self):
    n0, n1 = False, False
    for loc, nLoc, node in self.neighborLocations:
      if loc == 0.0:
        if n1:
          return False
        n0 = True
      elif loc == 1.0:
        if n0:
          return False
        n1 = True
    return True
    

  
  def lengthPerArea(self, _x1, _x2 = 0.5):
    """
    Compute length per cross sectional area for segment
    """
    if _x2 < _x1:
      _temp = _x1
      _x1 = _x2
      _x2 = _temp
    
    _lengths = [c.length for c in self.compartments]
    _cumLengths = list(cumsum(_lengths))
    def _findComp(_l):
      for n in range(len(_lengths)):
        if _cumLengths[n] >= _l:
          _c = self.compartments[n]
          _x = 1.0 - (_cumLengths[n] - _l) / _c.length
          return (n, _c, _x)
   
    # find 1st compartment, and proportion of distance across c1
    (_compInd1, _c1, _cX1) = _findComp(_x1 * _cumLengths[-1])
    # find 2nd compartment, and proportion of distance across c2
    (_compInd2, _c2, _cX2) = _findComp(_x2 * _cumLengths[-1])
      
    # compute length per area across distance
    if _compInd2 == _compInd1:
      return _c1.lengthPerArea(_cX1, _cX2)
    else:
      lPA = _c1.lengthPerArea(_cX1, 1.0) + _c2.lengthPerArea(0.0, _cX2)
      _compInd2 -= 1
      while _compInd2 > _compInd1:
        lPA += self.compartments[_compInd2].lengthPerArea(0.0, 1.0)
        _compInd2 -= 1
      return lPA
  
  def centroid(self, mandateTag=None):
    # return the centroid of a Segment, weighted by volume
    def _weightedC(_c):
      # return the weighted centroid of a compartment
      return tuple(_t * _c.volume / v for _t in _c.centroid)

    if mandateTag is None:
      v = self.volume
      centroids = [_weightedC(c) for c in self.compartments]
    else:
      v = sum(c.volume for c in self.compartments if mandateTag in c.tags)
      centroids = [_weightedC(c) for c in self.compartments if \
                   mandateTag in c.tags]
    
    # sum up the centroids
    centroid = (0.0, 0.0, 0.0)
    for pos in centroids:
      centroid = tuple(a + b for a,b in zip(centroid, pos))
    
    return centroid
  
  def centroidPosition(self, mandateTag=None):
    # return the position of the centroid of a Segment, weighted by volume
    #   (in this case, position is a number from 0 to 1, denoting proportion of
    #    distance along the segment)

    if mandateTag is None:
      halfV = self.volume / 2
      segLen = self.length
      centroidLen = 0.0
      for c in self.compartments:
        if c.volume < halfV:
          halfV -= c.volume
          centroidLen += c.length
        else:
          cFrac = halfV / c.volume
          centroidLen += c.length * cFrac
          break
    else:
      v = sum(c.volume for c in self.compartments if mandateTag in c.tags)
      halfV = v / 2
      segLen = self.length
      centroidLen = 0.0
      for c in self.compartments:
        if mandateTag in c.tags:
          if c.volume < halfV:
            halfV -= c.volume
            centroidLen += c.length
          else:
            cFrac = halfV / c.volume
            centroidLen += c.length * cFrac
            break
        else:
          centroidLen += c.length
    
    # return the position of the centroid
    return centroidLen / segLen
   

class Node:
  def __init__(self, _x, _y, _z, _r1, \
               _r2=None, _r3=None, _theta=0.0, _phi=0.0):
    self.x = _x
    self.y = _y
    self.z = _z
    self.r1 = _r1
    if self.r1 <= 0.0:
      if self.r1 < 0:
        raise ValueError('Encountered negative radius')
      else:
        raise ValueError('Encountered radius=0')
    if _r2 is None:
      if _r3 is not None:
        raise ValueError(\
          'Specify 1 radius for spherical nodes, 3 for ellipsoidal nodes')
      self.r2 = _r1
      self.r3 = _r1
    else:
      self.r2 = _r2
      self.r3 = _r3
    
    self.theta = _theta # angle from z axis
    self.phi = _phi # angle of azimuth (from x axis to semi-major axis)
    self.surface_area = None
    self.volume = None
    
    
    self.compartments = []
    self.segments = []
    self.tags = set()
  
  @property
  def maxRadius(self):
    return min(self.r1, self.r2, self.r3)
  
  @property
  def minRadius(self):
    return min(self.r1, self.r2, self.r3)
  
  @property
  def avgRadius(self):
    return (self.r1 * self.r2 * self.r3)**(1.0/3.0)
  
  def getElipse(self, node1):
    """
    Return details of elipse from intersection with compartment defined by
    connecting this node and node1
    """
    if self.r2 != self.r1 or self.r3 != self.r1:
      raise IOError('Currently cannot handle ellipsoidal nodes')
    
    # unnecessary when dealing with spheres:
    #(xAxis, yAxis, zAxis) = (node1.x - x, node1.y - y, node1.z - z)
    
    return (self.r1, self.r1, 0.0, self.x, self.y, self.z)



class Compartment:  
  def __init__(self):
    # do nothing, this is a pure virtual class
    self._surfaceArea = None
    self._volume = None
    self._length = None
    self.x0 = None
    self.y0 = None
    self.z0 = None
    self.x1 = None
    self.y1 = None
    self.z1 = None
    self.tags = set()
    self.name = None
    self.segment = None
    self.nodes = None
  
  @property
  def neighbors(self):
    """
    Find the compartment's neighbors
    """
    raise RuntimeError('Compartment is a pure virtual class')
  
  @property
  def length(self):
    if self._length is None:
      self._length = self._calcLength()
    return self._length
  
  def _calcLength(self):
    raise RuntimeError('Compartment is a pure virtual class')
  
  @property
  def surfaceArea(self):
    if self._surfaceArea is None:
      # compute surface area and convert from um^2 to mm^2
      self._surfaceArea = 1.0e-6 * self._calcSurfaceArea()
    return self._surfaceArea
  
  def _calcSurfaceArea(self):
    raise RuntimeError('Compartment is a pure virtual class')

  @property
  def volume(self):
    if self._volume is None:
      # compute volume and convert from um^3 to mm^3
      self._volume = 1.0e-9 * self._calcVolume()
    return self._volume
  
  def _calcVolume(self):
    raise RuntimeError('Compartment is a pure virtual class')
  
  def lengthPerArea(self, _x1, _x2 = 0.5):
    raise RuntimeError('Compartment is a pure virtual class')
  
  @property
  def maxRadius(self):
    return max(n.maxRadius for n in self.nodes)
  
  @property
  def minRadius(self):
    return min(n.minRadius for n in self.nodes)
  
  @property
  def avgRadius(self):
    """
    Return an estimate of average radius, pretending compartment was a cylinder
    or sphere (as appropriate to subclass)
    """
    raise RuntimeError('Compartment is a pure virtual class')
  
  @property
  def centroid(self):
    """
    return centroid of compartment as a tuple
    """
    raise RuntimeError('Compartment is a pure virtual class')


class OneNodeCompartment(Compartment):
  def __init__(self, node):
    
    # init Compartment object
    Compartment.__init__(self)
    
    self.nodes = [node]
    
    connectCompartment = node.compartments[0]
    if connectCompartment.node0 == node:
      node1 = connectCompartment.node1
    else:
      node1 = connectCompartment.node0
    (self.semiMajor, self.semiMinor, self.theta, self.x0, self.y0, self.z0) = \
      node.getElipse(node1)
    
    _direction = (node.x - node1.x, node.y - node1.y, node.z - node1.z)
    _norm = sqrt(sum([x*x for x in _direction]))
    _scale = self.length / _norm
    
    self.x1 = self.x0 + _scale * _direction[0]
    self.y1 = self.y0 + _scale * _direction[1]
    self.z1 = self.z0 + _scale * _direction[2]
    
    _scale *= 3.0/8.0
    self._centroid = (self.x0 + _scale * _direction[0], \
                      self.y0 + _scale * _direction[1], \
                      self.z0 + _scale * _direction[2])

  @property
  def node(self):
    return self.nodes[0]
  
  @property
  def neighbors(self):
    """
    Return a list of this compartments neighbors
    """
    return [comp for comp in self.node.compartments if comp != self]
  
  def _calcLength(self):
    """
    compute and set length
    """
    if self.node.r2 != self.node.r1 or self.node.r3 != self.node.r1:
      raise IOError('Currently cannot handle ellipsoidal nodes')
    
    return self.semiMajor
  
  def _calcSurfaceArea(self):
    """
    compute and set surface area
    """
    if self.semiMajor == self.semiMinor:
      return 4 * pi * self.semiMajor * self.semiMinor
    else:
      a = node.r1
      b = node.r2
      c = node.r3
      p = 1.6075
      # approximate formula accurate to relative error of <= 1.061%
      return 4 * pi * ( ((a*b)**p + (a*c)**p + (b*c)**p)/3 )**(1/p)
  
  def _calcVolume(self):
    """
    compute and set volume
    """
    return 4.0 * pi / 3.0 * self.node.r1 * self.node.r2 * self.node.r3

  def lengthPerArea(self, _x1, _x2 = 0.5):
    """
    compute and return length per area between relative positions _x1 and _x2
    0 <= _x <= 1, _x2 defaults to 0.5
    """
    if self.node.r2 != self.node.r1 or self.node.r3 != self.node.r1:
      raise IOError('Currently cannot handle ellipsoidal nodes')

    if _x2 < _x1:
      _temp = _x2
      _x2 = _x1
      _x1 = _temp
    
    return ( log((1.0 - _x1*_x1) / (1.0 - _x2*_x2)) /
             (2 * pi * self.node.r1) )
  
  @property
  def avgRadius(self):
    """
    Return an estimate of average radius, pretending compartment was a cylinder
    or sphere (as appropriate to subclass)
    """
    return (1.0e9 * self.volume / pi)**(1.0/3.0)
  
  @property
  def centroid(self):
    """
    return centroid of compartment as a tuple
    """
    return self._centroid



class TwoNodeCompartment(Compartment):
  def __init__(self, node0, node1):
    
    # init Compartment object
    Compartment.__init__(self)
    
    self.nodes = [node0, node1]
    
    (self.semiMajor0, self.semiMinor0, self.theta0, self.x0, self.y0, self.z0)\
       = node0.getElipse(node1)
    (self.semiMajor1, self.semiMinor1, self.theta1, self.x1, self.y1, self.z1)\
       = node1.getElipse(node0)
      
    self._centroid = None
  
  @property
  def node0(self):
    return self.nodes[0]
  
  @property
  def node1(self):
    return self.nodes[-1]
  
  @property
  def neighbors(self):
    """
    Return a list of this compartments neighbors
    """
    neighbors = [comp for comp in self.node0.compartments if comp != self]
    neighbors.extend(comp for comp in self.node1.compartments if \
                     comp != self and comp not in neighbors)
    return neighbors
  
  def _calcLength(self):
    """
    compute and set length
    """
    return sqrt((self.x1 - self.x0)**2 +
                (self.y1 - self.y0)**2 +
                (self.z1 - self.z0)**2)
  
  def _calcSurfaceArea(self):
    """
    compute and set surface area
    """
    ratio0 = self.semiMinor0 / self.semiMajor0
    ratio1 = self.semiMajor1 / self.semiMajor1
    if isnan(ratio0):
      if isnan(ratio1):
        raise IOError('Degenerate (zero radius) compartment')
      ratio0 = ratio1
    elif isnan(ratio1):
      ratio1 = ratio0
    elif ratio0 != ratio1:
      raise IOError('Don''t have formula for arbitrary eliptical frustrum')
    
    if self.semiMinor0 == self.semiMinor1:
      # cylinder
      coneFactor = 1.0
    elif self.length == 0:
      # disk with hole
      warn('Compartment with zero length')
      return pi * abs(self.semiMajor0 * self.semiMinor0
                      - self.semiMajor1 * self.semiMinor1)
    else:
      # cone
      coneFactor = \
        sqrt(1 + ((self.semiMinor0 - self.semiMinor1)/self.length)**2)
    
    if ratio0 == 1.0:
      # circular cross section
      angleFactor = pi
    else:
      # eliptical cross section
      angleFactor = 2.0*special.ellipe(sqrt((1.0 - ratio0*ratio0)) /coneFactor)
    
    return coneFactor * angleFactor * \
      self.length * (self.semiMajor0 + self.semiMajor1)
  
  def _calcVolume(self):
    """
    compute and set volume
    """
    return (pi / 3.0) * self.length * \
      (self.semiMajor0 * self.semiMinor0 + self.semiMajor1 * self.semiMinor1 +\
       0.5 * (self.semiMajor0 * self.semiMinor1 + \
              self.semiMajor1 * self.semiMinor0))

  def lengthPerArea(self, _x1, _x2 = 0.5):
    """
    compute and return length per cross sectional area between relative
    positions _x1 and _x2 (0 <= _x <= 1, _x2 defaults to 0.5)
    """
    if _x2 < _x1:
      _temp = _x2
      _x2 = _x1
      _x1 = _temp
    
    if self.semiMajor0 == self.semiMajor1:
      semiMajor = self.semiMajor0
      if self.semiMinor0 == self.semiMinor1:
        # cylinder
        semiMinor = self.semiMinor0
        coneFact = _x2 - _x1
      else:
        # semi-minor axis is changing, semi-major is constant
        semiMinor = 0.5 * (self.semiMinor0 + self.semiMinor0)
        minorRatio = semiMinor / (self.semiMinor1 - self.semiMinor0)
        coneFact = minorRatio * \
          log( (minorRatio + _x2 - 0.5) / (minorRatio + _x1 - 0.5) )
        if coneFact < 0:
          coneFact = -coneFact
    else:
      semiMajor = 0.5 * (self.semiMajor0 + self.semiMajor0)
      majorRatio = semiMajor / (self.semiMajor1 - self.semiMajor0)
      if self.semiMinor0 == self.semiMinor1:
        # semi-major axis is changing, semi-minor is constant
        semiMinor = self.semiMinor0
        coneFact = majorRatio * \
          log( (majorRatio + _x2 - 0.5) / (majorRatio + _x1 - 0.5) )
        if coneFact < 0:
          coneFact = -coneFact
      else:
        semiMinor = 0.5 * (self.semiMinor0 + self.semiMinor0)
        minorRatio = semiMinor / (self.semiMinor1 - self.semiMinor0)
        ratioProd = majorRatio * minorRatio
        avgRatio = 0.5 * (minorRatio + majorRatio)
        scale1 = (_x1 - 0.5) / ratioProd + avgRatio
        scale2 = (_x2 - 0.5) / ratioProd + avgRatio
        
        if avgRatio > 1.0:
          # answer in terms of logs
          ratioRoot = sqrt(avgRatio*avgRatio - 1.0)
          scale1 /= ratioRoot
          scale2 /= ratioRoot
          coneFact = 0.5 * ratioProd / ratioRoot * \
            log((scale2 - 1) * (scale1 + 1) / ((scale2 + 1) * scale1 - 1))
        elif avgRatio < 1.0:
          # answer in terms of atan
          ratioRoot = sqrt(1.0 - avgRatio*avgRatio)
          scale1 /= ratioRoot
          scale2 /= ratioRoot
          coneFact = ratioProd / ratioRoot * (atan(scale2) - atan(scale1))
        else:
          # answer in terms of 1/x
          coneFact = ratioProd / scale1 - ratioProd / scale2
    
    return coneFact * self.length / (pi * semiMajor * semiMinor)

  @property
  def avgRadius(self):
    """
    Return an estimate of average radius, pretending compartment was a cylinder
    or sphere (as appropriate to subclass)
    """
    try:
      return sqrt(1.0e9 * self.volume / self.length / pi)
    except ZeroDivisionError as err:
      if self.length == 0:
        return 0.5 * (self.semiMajor0 + self.semiMajor1)
      raise err
  
  
  @property
  def centroid(self):
    """
    Return centroid of the compartment as a tuple
    """
    if self._centroid is None:
      # need to calculate centroid location
      weightedLength = 0.5 * pi * self.length*self.length * \
        ((self.semiMajor0 * self.semiMinor0 + \
          self.semiMajor0 * self.semiMinor1 + \
          self.semiMajor1 * self.semiMinor0) / 6.0 + \
         0.5 * self.semiMajor1 * self.semiMinor1)
      ratio1 = weightedLength / (1.0e9 * self.volume) / self.length
      ratio0 = 1.0 - ratio1
      self._centroid = (self.x0 * ratio0 + self.x1 * ratio1, \
                        self.y0 * ratio0 + self.y1 * ratio1,
                        self.z0 * ratio0 + self.z1 * ratio1)
    
    return self._centroid
