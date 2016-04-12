# createHocFromMatlab.py - this takes the newSurfPoints output from
#   basicNeuronGeometryRadius.m and creates a hoc for Ted's code to read
#
# usage: python createHocFromMatlab.py surfPointsFileName 
#               skeletonPointsFile newHocFileName
#

import os, sys, re, math
import numpy as np



class MatNeuron:
  
  def __init__(self, surfFileName, skelFileName, hocFileName):
    self.inFile = surfFileName
    self.outFile = hocFileName
    self.skelFile = skelFileName
    self.skelSections = {}
    self.connections = [] # deal with this later
    #self.surfSections = {} # list of sections and their nodes
    self.rawAreas = {} # raw areas that will be averaged
    self.sumAreas = {} # summed areas
    self.radii = {} # dict of radii for each segment
    self.lengths = {}
    
    # load matlab file
    self.readsurfFile()
    self.readskelFile()
    self.getLengths()
    self.findRadius()
    self.writeFile()
    self.status()
  
  
  
  def readsurfFile(self):
    # read and parse input file from matlab
    # lines are of form:   x   y   z   area   secNum
    lineNum = 0
    
    with open(self.inFile, 'r') as fIn:
      # try
      
      for line in fIn:
        lineNum = lineNum + 1
        # print(lineNum)  # can comment out once it works
        splitLine = line.split(None)
        lineInt = [float(i) for i in splitLine]
        
        # if the current section doesn't already exist, create it
        if lineInt[4] not in self.rawAreas.keys():
         # self.surfSections[lineInt[4]] = [[[lineInt[i] for i in xrange(3)]]]
          self.rawAreas[lineInt[4]] = [[lineInt[3]]]
          
        else:
         # self.surfSections[lineInt[4]].append([[lineInt[i] for i in xrange(3)]])
          self.rawAreas[lineInt[4]].append([lineInt[3]])
    
    # after reading in the matlab file, create the output hoc file
    print('Read %i lines of surface points file' %lineNum)
    return self
  
  
  
  def readskelFile(self):
    # read the skelfile points
    lineNum = 0
    
    with open(self.skelFile, 'r') as sIn:
      
      for line in sIn:
        lineNum = lineNum + 1
        splitLine = line.split(None)
        lineFl = [float(i) for i in splitLine]
        
        # assign skelpoints to each section in dict 
        if lineFl[0] not in self.skelSections.keys():
          self.skelSections[lineFl[0]] = [[[lineFl[i] for i in xrange(1,4)]]]
        else:
          self.skelSections[lineFl[0]].append([[lineFl[i] for i in xrange(1,4)]])
    
    print('Read %i lines of skeleton points file' % lineNum)
    return self
  
  
  
  def getLengths(self):
    # get lengths based on euclidean distance
    print('Getting lengths of segments...')
    for sec in self.skelSections.keys():
      # loops through all the sections 
      currentDist = 0
      for h in xrange(len(self.skelSections[sec])-1):
        # loop through all the nodes for each section
        currentDist = currentDist + \
                      np.linalg.norm( np.array(self.skelSections[sec][h]) -
                                      np.array(self.skelSections[sec][h+1]))
        # all node distances summed
      
      if sec not in self.lengths.keys():
        self.lengths[sec] = currentDist
      else:
        print('Warning: sec num %i already in self.lengths!' %int(sec))
    # all sections summed
    
    return self
    
  
  
  def findRadius(self):
    # sum areas and calculate radius based on right cylinders
    print('Finding radii...')
    for sec in self.rawAreas.keys():
      # loop through all sections
      if sec not in self.sumAreas.keys():
        #print(self.rawAreas[sec])
        self.sumAreas[sec] = np.sum(self.rawAreas[sec])
      else:
        print('Mutliple areas added to section %i .' % int(sec))
        self.sumAreas[sec] = self.sumAreas[sec] + \
                             np.sum(self.rawAreas[key])
    
    # apply the formula to each section
    for sec in self.sumAreas.keys():
      #print(sec)
      #print(self.skelSections.keys())
      # SA = 2*pi*r*l , r = SA/l*pi*2
      self.radii[sec] = self.sumAreas[sec] / \
                        (self.lengths[sec] * np.pi * 2)
    
    return self
  
  
  
  def startSection(self, secNum):
    self.outfile.write('create section_%i\n' %secNum)
    self.outfile.write('section_%i {\n' %secNum)
    self.outfile.write('pt3dclear()\n')
    return
  
  def pointSection(self, x,y,z,rad):
    stringPt = [str(x), str(y), str(z), str(rad)]
    stringPt = ','.join(stringPt)
    self.outfile.write('pt3dadd(%s)\n' % stringPt)
    return
  
  def endSection(self):
    self.outfile.write('}\n\n')
    return
  
  
  
  def writeFile(self):
    # write the hocFile
    print('Writing the output file %s ...' % self.outFile)
    self.outfile = open(self.outFile, 'w')
    
    # write each section to the hoc file
    for sec in self.skelSections.keys():
      self.startSection(sec)
      try:
        rad = self.radii[sec]
      except:
        print('Missing radius')
        rad = 0
      
      # write each point to hoc file
      for pt in xrange(len(self.skelSections[sec])):
        pts = [i for i in self.skelSections[sec][pt]]
        # print(pts)
        self.pointSection(pts[0][0],pts[0][1],pts[0][2],rad)
      # all points for current section written
      
      self.endSection()
    # all sections written
  
  
  def status(self):
    # print some shit to show everything is done
    
    if len(self.skelSections.keys()) == len(self.radii.keys()) == \
                                 len(self.lengths.keys()):
      print('Everything looks good.\n')

        






def hocControl(surfFileName, skelFileName, hocFileName):
  # control the flow through the program
  geometry = MatNeuron(surfFileName, skelFileName, hocFileName)





############################################################
if __name__ == '__main__':
  arguments = sys.argv
  surfFileName = arguments[1]
  skelFileName = arguments[2]
  if len(arguments) > 3:
    hocFileName = arguments[3]
  else:
    hocFileName = 'newHocFile.txt'
  hocControl(surfFileName, skelFileName, hocFileName)
