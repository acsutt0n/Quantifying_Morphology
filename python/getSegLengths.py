# get the lengths of each segment and its centroid
# usage: python getSegLengths.py hocFile

import numpy as np


class HocThing:

    
  
  def __init__(self, peh):
    self.Connections = None
    self.secList = []
    self.secLengths = []
    self.secNodes = None
    self.openSection = 0
    self.currentSection = None
    self.currentNodes = None
    self.connectNum = 0
    
    if not self.secList:
#      print('received file %s' % hocFile)
      with open(hocFile, 'r') as fIn:
      # read the geometry file
        secCount = 0
        lineNum = 0
        for line in fIn:
#          print('reading %i' % lineNum)
          # loop through each line in the file
          splitLine = line.strip().split(None)
          try:
      
            splitLine = splitLine[0].split('_')
            
            if splitLine[0] == 'section':
              self.currentSection = int(splitLine[1])
              secCount = secCount + 1
              
            """
            Create secList section list
            """
            if self.secList and self.currentSection:
              self.secList = self.secList.append(currentSec)
            elif not self.secList and self.currentSection:
              self.secList = [currentSection]
            elif not self.currentSection:
              pass
          except:
            pass
              
          lineNum = lineNum + 1
 #       except:
 #         print('could not read line %i' % lineNum)
        
    self.readHoc()
#    self.connectSection()
    
    
    # end of def __init__
  
  
  def getNumSecs(hocFile):
    with open(hocFile, 'r') as fIn:
      # read the geometry file
      secCount = 0
      try:
        for line in fIn:
          # loop through each line in the file
          splitLine = line.strip().split(None)
          splitLine = splitLine[0].split('_')
          
          if splitLine[0] == 'section':
            currentSec = int(splitLine[1])
            secCount = secCount + 1
          
          if self.secList:
            self.secList = secList.append(currentSec)
          else:
            self.secList = [currentSec]
      except:
        return
      
      return self.secList
  
 
  
  def readHoc(self):
    """
    load hoc line to be sorted
    """
    lineNum = 0
    self.openSection = 0
    with open(hocFile, 'r') as fIn:
      # read the geometry file
#      try:
      for line in fIn:
#        print('reading line %i' % lineNum)
        # loop through each line in the file
        
        lineNum = lineNum + 1
        
        splitLine = line.split(None)
        self.parseHocLine(splitLine)
        
        if not splitLine:
          pass
#      except:
#        return
  
  
  def parseHocLine(self, splitLine):        
    """
    Find what each line is doing and send the line to be parsed
    """
#    print('reading line')
    try:
      if splitLine[0] == 'connect':
        # create a new connection
        self.connectSection(splitLine)
      elif splitLine[0].split('_')[0] == 'section' :
        # get current section, open a new section
        self.currentSection = self.createSection(splitLine)
      elif splitLine[0].split('(')[0] == 'pt3dadd':
        # log each node into the section's contents
        self.readSections(splitLine)
      elif splitLine[0] == '}' and self.openSection == 1:
 #       print('parseHocLine is calling closeSection')
        # close the section
        self.closeSection()
        self.openSection = 0
        self.currentSection = None
    except:
      pass

  
  def connectSection(self, splitLine):
#    print('connectSection called')
    if self.openSection == 1:
      self.openSection = 0
    seccols0, seccols1 = splitLine[1].split('_'), \
                         splitLine[2].split('_')
    
    # get section # and which end they connect to
    sec0, sec1 = int(seccols0[1].split('(')[0]), \
                 int(seccols1[1].split('(')[0])
    sec0end = int(seccols0[1].split('(')[1].split(')')[0])
    sec1end = int(seccols1[1].split('(')[1].split(')')[0])
    if sec0end == 0:
      connectEntry = [sec0, sec1]
    elif sec1end == 0:
      connectEntry = [sec1, sec0]
    else:
      print('bad sections found in connection matrix')
#    print(connectEntry)
    
    # if this connection doesn't already exist, create it
    if self.Connections is None:
#      print('starting Connections')
      self.Connections = {0: connectEntry}
      self.connectNum = 1
    elif self.Connections is not None and connectEntry in self.Connections.values():
      print('connection listed twice')
    else:
#      print('adding connection')
      self.Connections[self.connectNum] = [connectEntry]
      self.connectNum = self.connectNum + 1
  
  
  def createSection(self, splitLine):
    """
    set openSection = 1 and get the current section
    """
    self.currentSection = int(splitLine[0].split('_')[1])
    self.openSection = 1
    return self.currentSection
    # open the section for creation
    
    
  def readSections(self, splitLine):
    """
    Now create the actual sections as numpy arrays:
    Nodes={sec1: [(x1,y1,z1,r1), (x2,y2,z2,r2),..., sec2: [...]}
    """
    # NP array syntax: nodes = np.array( [[ x, y, z, r ]] )
    #                  nodes = np.append( nodes, [[ x2,y2,z2,r2 ]], 0 )
    # log each point
#    print(splitLine)
    pt3dcols = splitLine[0].split('(')[1]

    pt3dcols = pt3dcols.split(',')

    coords = [float(pt3dcols[c]) for c in [0,1,2]]
    rad = float(pt3dcols[3].split(')')[0])
    
    if self.currentNodes is None:
      # section doesn't yet exist
      self.currentNodes = np.array([[coords[0],coords[1],coords[2],rad]])
    elif self.currentNodes is not None:
      self.currentNodes = np.append( self.currentNodes,
                          [[coords[0],coords[1],coords[2],rad]], 0)
    

    
    # create dictionary of nodes for each section while 
    # openSection == 1
    # then close and move on to next currentSection

  
  def closeSection(self):
    """
    Create the dictionary of node lists
    """
#    print('closeSection called')
    if len(self.currentNodes) < 1:
      print('Bad section found')
    else:
      pass
    
    if self.secNodes is not None and self.currentSection in self.secNodes.keys():
      print('Tried to add node list for %i which already exists' 
            % self.currentSection)
    else:
      pass
    
#    print(self.currentSection, len(self.currentNodes))
          
    if self.secNodes is None and self.currentSection is not None:
#      print('Creating secNodes at %i' % self.currentSection)
      # add the current section to newly created dict object
      self.secNodes = {self.currentSection: self.currentNodes}
    elif self.secNodes is not None and self.currentSection is not None:
#      print('Added section %i to secNodes' % self.currentSection)
      self.secNodes[self.currentSection] = self.currentNodes
    else:
      print('Error in assigning nodes to %i' % self.currentSection)
    
#    print('section closed')
    # close section
    self.openSection = 0
    self.currentSection = None
    self.currentNodes = None
#    if self.currentNodes is None:
#      print('currentNodes erased')
#    else:
#      print('currentNodes has %i nodes, should have zero' 
#            % len(self.currentNodes))
     
  
def printIfExist(obj):
  if obj is not None:
    if type(obj) is int:
      return obj
    else:
      return len(obj)
  else:
    try:
      return len(obj)
    except:
      return 0



########################################################################
#                          control
def getLengths(hocFile):
  neuron = HocThing(hocFile)
  attrs = vars(neuron)
  print(', '.join("%s: %s" % item for item in attrs.items()))      ###
#  print(type(neuron))
  lengths = [(len(neuron.secNodes[lgs]), neuron.secNodes.keys()[lgs]) 
                   for lgs in range(len(neuron.secList))]
#  print('Section numbers and their lengths:')
  for c in range(len(neuron.secList)):
    print(secs[c], lengths[c])
  if neuron.currentNodes is not None:
    print('There are still %i nodes in currentNodes' % len(neuron.currentNodes))
  # then do more stuff after, create a HocThing object
  
  ## print summary
  if neuron.secNodes is not None:
    print(len(neuron.secNodes))







###############################################################################
if __name__ == "__main__":
  import sys
  hocFile = sys.argv[1]
  getLengths(hocFile)
  # exit
  sys.exit(0)
