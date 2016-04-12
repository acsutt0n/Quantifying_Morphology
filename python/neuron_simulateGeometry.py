#!/usr/bin/python
"""
This script provides resources for simulating a passive neuron model on a
NeuronGeometry
"""

from NeuronGeometry import *
from neuron_readExportedGeometry import HocGeometry
import scipy
import peelLength
from matplotlib import pyplot
import sys
import json


###############################################################################
def _runningOnSubmitNode():
  """
  return True if running on grid engine submit node, false otherwise
  """
  if sys.version_info[0] == 3:
    import subprocess
  else:
    import commands as subprocess
  retVal, outputStr = subprocess.getstatusoutput('qstat')
  return (retVal == 0)




###############################################################################
def getPassiveProperties(passiveFile=None, parameters=None):
  if passiveFile is not None and os.access(passiveFile, os.R_OK):
    with open(passiveFile, 'r') as fIn:
      properties = json.load(fIn)
    return properties
  elif parameters is None:
    parameters = (100.0, 1.0, 1.29e-3,
                  100.0, 1.0, 2.7e-6,
                  100.0, 1.0, 2.7e-6);
  
  properties = [
    {
      'matchProp' : 'branchOrder',  # soma has branchOrder == 0
      'matchVal' : 0,
      'name' : 'Soma',
      'values' : {
        'Ra' : parameters[0],
        'cm' : parameters[1]
      },
      'channels' : {
        'pas' : {
          'g' : parameters[2],
          'e' : 0.0
        }
      }
    },
    {
      'matchProp' : 'isTerminal',  # tip is not soma, but isTerminal
      'matchVal' : True,
      'name' : 'Tip',
      'values' : {
        'Ra' : parameters[3],
        'cm' : parameters[4]
      },
      'channels' : {
        'pas' : {
          'g' : parameters[5],
          'e' : 0.0
        }
      }
    },
    {
      'matchProp' : None,  # match the remainder
      'matchVal' : None,
      'name' : 'Middle',
      'values' : {
        'Ra' : parameters[6],
        'cm' : parameters[7]
      },
      'channels' : {
        'pas' : {
          'g' : parameters[8],
          'e' : 0.0
        }
      }
    }
  ]
  return properties


###############################################################################
def makeModel(geometry, properties):
  # given geometry and passive properties, make an object holding model info
  
  def _setProperties(segment):
    for neuronSection in properties:
      if 'matchTag' in neuronSection:
        if neuronSection['matchTag'] in segment.tags:
          # this is the default section, match anything:
          segment.tags.add(neuronSection['name'])
          return neuronSection['values'], neuronSection['channels']
        else:
          continue
      elif neuronSection['matchProp'] is None:
        # this is the default section, match anything:
        segment.tags.add(neuronSection['name'])
        return neuronSection['values'], neuronSection['channels']
      segmentVal = getattr(segment, neuronSection['matchProp'])
      if (hasattr(neuronSection['matchVal'], "__len__") and
          segmentVal in neuronSection['matchVal']) or\
         segmentVal == neuronSection['matchVal']:
        # this section matches segment's properties, so return it
        segment.tags.add(neuronSection['name'])
        return neuronSection['values'], neuronSection['channels']
    raise RuntimeError('No match for segment')
    
  model = {
    'stimulus' : { 'amplitude' : 1.0, # nA
                   'duration'  : 1000, # ms
                   'delay'     : 100, # ms
                   'segment'  : geometry.soma,
                   'location' : geometry.soma.centroidPosition(mandateTag=
                                                               'Soma')
                   
                 },
    'tFinal'   : 1800, # ms
    'dT'       : 0.2, # ms
    'v0'       : 0.0, # mV
    'properties' : _setProperties
  }
  return model


    

###############################################################################
def _simulateModel(geometry, model, child_conn=None):
  """
  Do the nuts and bolts of neuron simulation, typically called in a separate
  Process by simulateModel
  """
  ##-------------------------------------------------------------------------##
  def _createSegment(segment, geometry):
    segName = segment.name
    if '[' in segName and ']' in segName:
      ind1 = segName.index('[')
      ind2 = segName.index(']')
      baseName = segName[:ind1]
      index = int(segName[ind1+1:ind2])
      if index == 0:
        # this is the first time this baseName is used, create all at once
        maxInd = 0
        for seg in geometry.segments:
          if seg.name.startswith(baseName + '['):
            segIndex = int(seg.name[ind1+1:-1])
            maxInd = max(maxInd, segIndex)
        neuron.h('create %s[%d]' % (baseName, maxInd + 1))
      # get the hoc segment object, and add it to NeuronGeometry segment object
      segment.hSeg = getattr(neuron.h, baseName)[index]
    else:
      neuron.h('create %s' % segName)
      segment.hSeg = getattr(neuron.h, segName)
    return segment.hSeg  
  ##-------------------------------------------------------------------------##
  def _getHSeg(segName):
    if '[' in segName and ']' in segName:
      ind1 = segName.index('[')
      ind2 = segName.index(']')
      baseName = segName[:ind1]
      index = int(segName[ind1+1:ind2])
      return getattr(neuron.h, baseName)[index]
    else:
      return getattr(neuron.h, segName)
  ##-------------------------------------------------------------------------##
  def _addGeometryToHoc(geometry):
    # define all the segments within hoc
    for segment in geometry.segments:
      # create the segment, get hoc segment object
      hSeg = _createSegment(segment, geometry)
      # set it as the currently active segment
      hSeg.push()
      # add all the nodes to define the segment geometry
      for node in segment.nodes:
        neuron.h.pt3dadd(node.x, node.y, node.z, 2 * node.r1)
      # no longer editing this segment
      neuron.h.pop_section()
      
    # connect the segments
    for index, segment in enumerate(geometry.segments):
      hSeg = segment.hSeg
      # to avoid repeating the same connections:
      #  find all the nodes with neighbors, and get the neighbor segment indexes
      #  only specify a connetion if index > all neighbor segment indexes
      neighborNodes = { node for location, nLocation, node
                        in segment.neighborLocations }
      connectNodes = []
      for node in neighborNodes:
        neighborInds = [geometry.segments.index(neighbor) for neighbor in
                        node.segments]
        if min(neighborInds) == index:
          # make all the connections
          connectNodes.append(node)
      
      if not connectNodes:
        continue
      
      for neighbor, (location, nLocation, node) in zip(segment.neighbors,
                                                      segment.neighborLocations):
        if node not in connectNodes:
          # don't add this connection now
          continue
        nSeg = neighbor.hSeg
        nSeg.connect(hSeg, location, nLocation)
  ##-------------------------------------------------------------------------##
  def _setProperties(geometry, model):
    
    # first find branch orders, because they can be used to target properties
    if geometry.soma.branchOrder is None:
      geometry.calcBranchOrder(doPlot=False)
    
    for segment in geometry.segments:
      properties, channels = model['properties'](segment)
      hSeg = segment.hSeg
      for prop, val in properties.items():
        setattr(hSeg, prop, val)
      for channel, chanPropDict in channels.items():
        hSeg.insert(channel)
        for prop, val in chanPropDict.items():
          setattr(hSeg, prop + '_' + channel, val)
  ##-------------------------------------------------------------------------##
  def _initStimulusAndRecording(geometry, model):
    # set up the stimulus
    stimInfo = model['stimulus']
    iClamp = neuron.h.IClamp(stimInfo['segment'].hSeg(stimInfo['location']))
    iClamp.amp = stimInfo['amplitude']
    iClamp.dur = stimInfo['duration']
    iClamp.delay = stimInfo['delay']
    
    # record voltage in each segment
    vTraces = {segment.name : neuron.h.Vector()
               for segment in geometry.segments}
    for segment in geometry.segments:
      trace = vTraces[segment.name]
      hSeg = segment.hSeg
      hPos = 0.5 # for now, just record in the middle of all segments
      trace.record(hSeg(hPos)._ref_v, model['dT'])
    
    return iClamp, vTraces
  ##-------------------------------------------------------------------------##
  def _runSimulation(model):
    neuron.h.dt = model['dT']
    neuron.h.finitialize(model['v0'])
    neuron.h.fcurrent()
    tFinal = model['tFinal']
    while neuron.h.t < tFinal:
      neuron.h.fadvance()


  # redirect stdout and stderr
  import os
  import tempfile
  if sys.version_info[0] == 3:
    from io import StringIO
  else:
    from cStringIO import StringIO
  import traceback
  
  def _redirect_stdout():
    tempStdFile = tempfile.NamedTemporaryFile(delete=False).name
    tempStdOut = os.open(tempStdFile, os.O_WRONLY)
    sys.stdout.flush()
    sys.stderr.flush()
    temp1 = os.dup(1)
    temp2 = os.dup(2)
    os.dup2(tempStdOut, 1)
    os.dup2(tempStdOut, 2)
    os.close(tempStdOut)
    return tempStdFile, temp1, temp2
  
  textOutput = StringIO()
  stdout = sys.stdout ; stderr = sys.stderr
  sys.stdout = textOutput ; sys.stderr = textOutput
  tempStdOutFile, temp1, temp2 = _redirect_stdout()
  #temp1 = None ; temp2 = None ; tempStdOutFile = ""
  err = None ; tb = ""
  
  try:
    import neuron
    
    geometry.checkConnectivity(removeDisconnected=True, removeLoops=True)
    firstSeg = geometry.segments[0]
    if hasattr(firstSeg, 'hSeg') and firstSeg.hSeg is not None:
      for segment in geometry.segments:
        segment.hSeg = None
    _addGeometryToHoc(geometry)
    _setProperties(geometry, model)
    iClamp, vTraces = _initStimulusAndRecording(geometry, model)
    _runSimulation(model)
    # convert traces to python arrays
    for segment in geometry.segments:
      vTraces[segment.name] = scipy.array(vTraces[segment.name])
    # create time trace
    numT = len(vTraces[firstSeg.name])
    timeTrace = scipy.array([n * model['dT'] for n in range(numT)])
          
  except BaseException as err:
    timeTrace = [] ; vTraces = {}
    tb = traceback.format_exc()

  if temp1 is not None:
    os.dup2(temp1, 1)
    os.close(temp1)
  if temp2 is not None:
    os.dup2(temp2, 2)
    os.close(temp2)

  textOutput = textOutput.getvalue()
  # Get the output of stderr, stdout and clean up temporary files
  if tempStdOutFile:
    with open(tempStdOutFile, 'r') as fOut:
      textOutput += fOut.read()
    os.remove(tempStdOutFile)
  sys.stdout = stdout ; sys.stderr = stderr

  # report results
  if child_conn is None:
    # if called directly, return results  
    return timeTrace, vTraces, textOutput, err, tb
  else:
    # otherwise send results back via pipe
    child_conn.send((timeTrace, vTraces, textOutput, err, tb))
    child_conn.close()


###############################################################################
def simulateModel(geometry, model):
  from multiprocessing import Pipe, Process
  from time import sleep
  parent_conn, child_conn = Pipe()
  p = Process(target=_simulateModel, args=(geometry, model, child_conn))
  try:
    p.start()
    while not parent_conn.poll():
      sleep(0.1)
    timeTrace, vTraces, textOutput, err, tb = parent_conn.recv()
    p.join()
  except BaseException:
    if p.is_alive():
      p.terminate()
    timeTrace = [] ; vTraces = [] ; textOutput = "" ; err = None ; tb = ""
    raise
  
  if err is not None:
    print(tb)
    raise err
  return timeTrace, vTraces, textOutput


###############################################################################
def plotTraces(timeTrace, vTraces):
  pyplot.figure()
  for segment, trace in vTraces.items():
    pyplot.plot(timeTrace, trace)
  pyplot.show()


###############################################################################
def _parseArguments():
  import argparse
  parser = argparse.ArgumentParser(description=
    "Simulate a neuron with geometry exported in a .hoc file, and passive "
    + "properties specified in a separate json .txt file",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("geoFile", help="file specifying neuron geometry",
                      type=str)
  parser.add_argument("passiveFile", nargs="?",
                      default="passive_properties.txt",
                      help="file specifying passive properties", type=str)
  return parser.parse_args()


###############################################################################
if __name__ == "__main__":
  # get the geometry file
  options = _parseArguments()
  # create geometry from the file
  geometry = HocGeometry(options.geoFile)
  # get passive properties
  properties = getPassiveProperties(options.passiveFile)
  # make neuron model
  model = makeModel(geometry, properties)
  # simulation model on specified geometry
  timeTrace, vTraces, textOutput = simulateModel(geometry, model)
  if textOutput:
    print(textOutput.rstrip())
  # peel length
  model, vErr, vResid = \
    peelLength.modelResponse(timeTrace, vTraces[geometry.soma.name],
                             verbose=False, findStepWindow=True,
                             plotFit=False, debugPlots=False,
                             displayModel=False)
  print('Fitting exponentials:')
  peelLength.printModel(model, vErr=vErr, vResid=vResid)
  
  # plot results
  plotTraces(timeTrace, vTraces)
  
  #exit
  sys.exit(0)
