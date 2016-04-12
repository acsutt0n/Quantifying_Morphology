# amiraGeometry.py -- creates a class that can import all information
#                 from an Amira spatial graph file
#                 Can then pass it to another function for whatever
#                 (generally for import, but can return a features 
#                 or a hoc if needed)
# usage: python amiraGeometry.py amiraFile options
#        options: save as hoc (-hoc)
#                 return simple features; num of segment, etc (-summary)

import sys, os
import math
import numpy as np



class SptGr:
  
  def __init__(self, amiraFileName):
    self.inFile = amiraFileName
    self.endpoints = []
    self.getEndpoints = False
    self.nodes = []
    self.getNodes = False
    self.lengths = []
    self.getLengths = False
    self.rads = []
    self.getRads = False
    self.connections = []
    self.getConnections = False
    self.sections = {}
    self.secRads = {}
    self.loops = []
    self.getLoops = False
    self.properties = {}
    self.loopSecs = []
    self.neighbors = {}
    
    # decide what to do here
    self.readAmiraFile()
    #self.sanitizeLists()
    self.createSections()
    self.findNeighbors()
    self.checkEndpoints()
    self.sectionRads()
    self.printSummary()
  
  
  
  def readAmiraFile(self):
    # read lines from spatial graph file (may need sub-functions)
    
    lineNum = 0
    print('Loading file %s ...' % self.inFile)
    
    with open(self.inFile, 'r') as fIn:
      # try
      
      for line in fIn:
        lineNum = lineNum + 1
        
        splitLine = line.split(None)
        if splitLine:

          if splitLine[0] == 'Parameters':
            self.loopsOn = True
          elif splitLine[0] == 'define':
            self.properties[splitLine[1]] = int(splitLine[2])
          elif splitLine[0] == '@1':
            self.getEndpoints = True
          elif splitLine[0] == '@2':
            self.getEndpoints = False
            self.getConnections = True
          elif splitLine[0] == '@3':
            self.getConnections = False
            self.getLengths = True
            #print('Value of self.getLengths: %s' % self.getLengths)
          elif splitLine[0] == '@4' and self.loopsOn == True:
            self.getLoops = True
            self.getLengths = False
          # @4 and @5 depend on whether there are loops present
          elif (splitLine[0] == '@4' and self.loopsOn == False) or \
               (splitLine[0] == '@5' and self.loopsOn == True):
            print('Lengths retrieved: %i' % len(self.lengths))
            self.getLoops = False
            self.getNodes = True
          elif (splitLine[0] == '@5' and self.loopsOn == False) or \
               (splitLine[0] == '@6' and self.loopsOn == True):
            print('Number of nodes retrieved: %i' % len(self.nodes))
            self.getNodes = False
            self.getRads = True
          
          elif self.getEndpoints == True:
            self.endpoints.append([float(i) for i in splitLine])
          elif self.getConnections == True:
            self.connections.append([int(i) for i in splitLine])
          elif self.getLengths == True:
            self.lengths.append(int(splitLine[0]))
          elif self.getLoops == True:
            self.loopSecs.append(int(splitLine[0]))
          elif self.getNodes == True:
            self.nodes.append([float(i) for i in splitLine])
          elif self.getRads == True:
            self.rads.append([float(splitLine[0])])
          
        
    print(' done. Lengths is %s items' % len(self.lengths))
    
    return self
  
  
  
  def sanitizeLists(self):
    """
    remove empty [] elements from the lists; these arise from
    newline characters that are included due to the \n between
    the last element in the sptgr and the next section ('@...')
    """
    print('Sanitizing lists ...')
    # find and dispose of all empty indices in self.lengths
    mt_inds = [i for i,x in enumerate(self.lengths) if x == [] ]
    for i in xrange(len(mt_inds)):
      self.lengths.pop(mt_inds[i - i])
    # turn eliminate "second" dimension of self.lengths [[x],[y]] -> [x,y]
    self.lengths = [self.lengths[i][0] 
                    for i in xrange(len(self.lengths))]
      
    # find and dispose of all empty indices in self.nodes
    mt_inds = [i for i,x in enumerate(self.nodes) if x == [] ]
    for i in xrange(len(mt_inds)):
      self.nodes.pop(mt_inds[i - i])
      
    # find and dispose of all empty indices in self.endpoints
    mt_inds = [i for i,x in enumerate(self.endpoints) if x == [] ]
    for i in xrange(len(mt_inds)):
      self.endpoints.pop(mt_inds[i - i])
      
    # find and dispose of all empty indices in self.connections
    mt_inds = [i for i,x in enumerate(self.connections) if x == [] ]
    for i in xrange(len(mt_inds)):
      self.connections.pop(mt_inds[i - i])
      
    # find and dispose of all empty indices in self.rads
    mt_inds = [i for i,x in enumerate(self.rads) if x == [] ]
    for i in xrange(len(mt_inds)):
      self.rads.pop(mt_inds[i - i])
    self.rads = [self.rads[i][0] for i in xrange(len(self.rads))]
    print(' done.')
    
    return self



  def createSections(self):
    # use lengths and nodes to create node-detailed sections
    if len(self.sections) > 0:
      print('Warning: sections already exist in self.sections!')
   # if [] in self.lengths:
   #   self.sanitizeLists()
    # print(self.lengths)
    print('Creating sections with specific nodes ... ')
    for sec in xrange(len(self.lengths)):
      startnode = int( np.sum(self.lengths[:sec]) )
      endnode = int( startnode + self.lengths[sec] )
      self.sections[sec] = self.nodes[startnode:endnode]
    print(' done. ')
    
    return self
        
    

  def whereConnect(self, seg0, seg1):
    # find where the two segments are connected (arbitrary based on
    # location of node in the list
    node00, node01 = self.sections[seg0][0], \
                   self.sections[seg0][-1]
    node10, node11 = self.sections[seg1][0], \
                   self.sections[seg1][-1]
    print(node00, node01, node10, node11)
    if node00 in [node10,node11]:
      if node00 == node10:
        return [0,0]
      else:
        return [0,1]
    if node01 in [node10,node11]:
      if node01 == node10:
        return [1,0]
      else:
        return [1,1]
    else:
      return False



  def findNeighbors(self):
    # find neighbors from connections list and overlappying node tuples
    # structure is self.neighbors[seg] = [neighborseg[end], ...]
    
    if len(self.neighbors) > 0:
      print('Warning: self.neighbors already populated with %i entries' 
            % len(self.neighbors) )
    
    # print(self.connections)
    
    # first entry in list is probably real and not an island
    for cxn in xrange(len(self.connections)):
      seg0, seg1 = self.connections[cxn][0], self.connections[cxn][1]
      connects = self.whereConnect(seg0, seg1)
      if connects is not False:
        # if an entry does not exist for a section, create it
        if seg0 not in self.neighbors.keys():
          self.neighbors[seg0] = [seg1,[connects[0]]]
        else:
          self.neighbors[seg0].append([seg1,connects[0]])
        if seg1 not in self.neighbors.keys():
          self.neighbors[seg1] = [seg0, connects[1]]
        else:
          self.neighbors[seg1].append([seg0,connects[1]])
      else:
        print('No valid connection betweel segs %i and %i.'
             % (seg0, seg1) )
    
    # make sure all neighbors are represented for both segments
    #for sec in self.neighbors.keys():
    #  for neb in self.neighbors[sec]:
    #    if self.neighbors[sec][neb][0] not in self.neighbors.keys():
    #      self.neighbors[self.neighbors[sec][neb][0]] = \
    #        [sec, se
    
    return self



  """
  def returnNeighbors(self, segnum):
    # call with self a chosen segment number, will return all neighbors
    inds0 = [i for i, x in enumerate(self.connections) if x[0]==segnum]
    inds1 = [i for i, x in enumerate(self.connections) if x[1]==segnum]
    neighborsegs = [self.connections[i][0] for i in inds0] + \
                   [self.connections[i][1] for i in inds1] 
    
    return neighborsegs
  """
  
  
  
  def checkEndpoints(self):
    """
    Make sure that the class of [0] and [-1] elements of skelpoints
    is the same as the endpoints list
    """
    print('Checking endpoints...')
    allendpoints = []
    # print('keys of self.sections are %s' %self.sections.keys())
    for seg in self.sections.keys():
      # print('Populating section %i' % seg)
      allendpoints.append([self.sections[seg][0]])
      allendpoints.append([self.sections[seg][-1]])
    # keep only unique endpoints
    allendpoints = np.array(allendpoints)
    
    
    if len(np.shape(allendpoints))>2:
      allendpoints = np.array([allendpoints[i][0] for \
                               i in xrange(len(allendpoints))])
    print(allendpoints[:10])
    # print(len(allendpoints))
    print(np.shape(allendpoints))
    
    ends = np.ascontiguousarray(allendpoints).view(np.dtype((np.void, \
           allendpoints.dtype.itemsize * allendpoints.shape[1])))
    _, idx = np.unique(ends, return_index=True)
    unique_endpoints = np.copy(allendpoints[idx])
    
    if len(self.endpoints) == len(unique_endpoints):
      selfendpoints = self.endpoints[self.endpoints[:,0].argsort()]
      uniqueendpoints = unique_endpoints[unique_endpoints[:,0].argsort()]
      boolary = selfendpoints == uniqueendpoints
      if False in boolary:
        print('Not all endpoint tuples are equal. Running elementwise.')
        for pt in self.endpoints:
          if pt not in unique_endpoints:
            print('unique_endpoints missing pt (%.5f, %.5f, %.5f)' \
                  % (pt[0],pt[1],pt[2]))
      else:
        print('All endpoints present in both arrays.')
    else:
      print('Unequal numbers of endpoints. Running elementwise.')
      for pt in self.endpoints:
        if pt not in unique_endpoints:
          print('unique_endpoints missing pt (%.5f, %.5f, %.5f)' \
                % (pt[0],pt[1],pt[2]))
    # no errors means seg[0] and [-1] can be used for endpoints
    print('Done checking endpoints.')
    
    return self
  
  
  
  def ptmode(self, ptlist):
    # find the most common point in a list of points (N x 3)
    countpt = [ptlist.count(i) for i in ptlist]
    winneris = max(countpt)
    loc = countpt.index(winneris)
    val = ptlist[loc]
    
    return self, val
  
  
  
  def radiusOf(self, secnum):
    # returns the radius at the ends of a section
    neighb0 = [self.neighbors[secnum][i][0] for i in 
               xrange(len(self.neighbors[secnum])) if 
               self.neighbors[secnum][i][1] == 0]
    neighb1 = [self.neighbors[secnum][i][0] for i in 
               xrange(len(self.neighbors[secnum])) if 
               self.neighbors[secnum][i][1] == 1]
    rads0 = [self.secRads[i][0] for i in neighb0]
    for i in xrange(len(neighb0)):
      rad0.append(self.secRads[i][-1])
    rads1 = [self.secRads[i][0] for i in neighb1]
    for i in xrange(len(neighb1)):
      rad1.append(self.secRads[i][-1])
    
    rad0 = [rad0[i] for i in xrange(len(rad0)) if rad0[i] > 0]
    rad1 = [rad1[i] for i in xrange(len(rad1)) if rad1[i] > 0]
      
    consensus0, consesnsus1 = self.ptmode(rad0), self.ptmode(rad1)
    
    return rad0, rad1
  
  
  
  def sectionRads(self):
    # assign radiuses to the sections
    if len(self.secRads) > 0:
      print('Radiuses already assigned')
    print('Creating radiuses ...')
    for sec in xrange(len(self.lengths)):
      startnode = int( np.sum(self.lengths[:sec]) )
      endnode = int( startnode + self.lengths[sec] )
      self.secRads[sec] = self.rads[startnode:endnode]
    print('SecRads draft created with %i rads.' % len(self.secRads))
    
    for sec in self.secRads.keys():
      self.secRads[sec] = [self.secRads[sec][i][0] \
                           for i in xrange(len(self.secRads[sec]))]
    
    meanrad = np.mean(self.rads)
    for sec in self.secRads.keys():
      for node in xrange(len(self.secRads[sec])):
        # print(self.secRads[sec])
        if self.secRads[sec][node] > 10000 * meanrad or \
          self.secRads[sec][node] <= 0:
            
            rad0, rad1 = self.radiusOf(sec)
            self.secRads[sec] = [rad0, rad1]
    
    print(' done.')
    return self
  
  
  
  def printSummary(self, options=False):
    """
    Print some basic results, if options = False a bunch of normal shit
    is printed.
    """
    # some properties
    meanrad = np.mean(self.rads)
    meanlength = np.mean(self.lengths)
    numsegs = len(self.sections.keys())
    numnodes = len(self.nodes)
    
    if options == False:
      print('Properties of %s:' % self.inFile)
      print('%i nodes' % numnodes)
      print('%i segments' % numsegs)
      print('Mean segment length: %.5f' % meanlength)
      print('Mean segment width: %.5f' % meanrad)
  
  
    
    
################################
if __name__ == '__main__':
  arguments = sys.argv
  amiraFileName = arguments[1]
  if len(arguments) > 2:
    options = arguments[2]
  else:
    options = False
  geometry = SptGr(amiraFileName)


