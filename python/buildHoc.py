# Builds a hoc file (in Amira format) from an SWC file, such as
#    that created by the swc_tree.m function in TREES toolbox

# Usage: python buildHoc.py swcFileName

import os, sys, re, math


class HocThing:
  
  def setFileName(self, fileName):
    self.fileName = fileName
    self.name = fileName.split('.')[0]
  
  def __init__(self, _fileName=None):
    self.fileName = None
    self._connections = []
    self._name = None
    self._nodes = []
    self._radius = []
    self._center = []
    
    if _fileName is not None:
      self.setFileName(_fileName)
      self.readSWC()
    
    
  def _parseSWC(self, line):
    """
    Read line from SWC file and create a new hoc entry
    """
    splitLine = line.split(None) # same as split(' ')
    if not splitLine:
      return
    
    else:
      
  
  
  
    
  def readSWC(self):
    
    lineNum = 0
    with open(self.fileName, 'r') as fIn:
      # read swc file
      try:
        for line in fIn:
          # loop through each line in the file
          # inc line number
          lineNum = lineNum + 1
          # parse the line in geometry file, add info to info
          self._parseSWC(line)






















def runHoc(swcFile):
  hocThing = HocThing(swcFile)
  
  











#################################################################
if __name__ == "__main__":
  import sys
  swcFile = sys.argv[1]
  runHoc(swcFile)
