# hocNodeSections.py - turn each node into its own section using 
#                      float rather than integer section numbers

# usage: python hocNodesections.py hocFile outFile (optional)

import numpy as np
import os, sys, re, math
from NeuronGeometry import *













def hocNodeControl(hocFile, outFile):
  # Control the flow of the program
  



##############################
if __name__ == "__main__":
  if len(sys.argv) > 2:
    outFile = sys.argv[2]
  else:
    outFile = 'hocNodeOutFile.txt'
  hocFile = sys.argv[1]
  hocNodeControl(hocFile, outFile)
  
  





"""

class HocNeuron:
  def __init__(self, hocFileName=None):
    self.connections = []
    self.intSections = {} # dictionary of lists of nodes for each int
                          # section
    if hocFileName is not None:
      self.setFileName(
    readhocFile(hocFile)



"""
