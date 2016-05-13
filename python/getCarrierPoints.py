# find 'carrier points' - places where dendrites terminate
# usage: python getCarrierPoints.py neuronFile.hoc


def getCarrierPoints(hocFile):
  return  


class Section:
  def __init__(self, secName):
    self.secName = secName
    self.nodes = dict()
    self.neighbors = []
  def addPoint(self, point):
    self.nodes[len(self.nodes)] = point
  def addNeighbor(self, neighborSec):
    append.self.neighbors(neighborSec)
  def returnNodes(self):
    return self.nodes[0] self.nodes[-1]
    
  
# ------------- parse connections ---------
# this is called first (?? maybe later ??)
def parseConnections(hocFile):
  with open(hocFile, 'r') as infile:
    zerothEnd = []
    onethEnd = []
    for line in infile:
      cols = line.strip().split('_')
      
      if cols[0] == 'connect section':
        current_zeroth = cols[1].split('(')
        zerothEnd.append(float(current_zeroth[0]))
        current_oneth = cols[2].split('(')
        onethEnd.append(float(current_oneth[0]))
      else:
        continue
  return zerothEnd, onethEnd
  
  
# ------------- parse sections -----------
# this is called second
def parseSections(hocFile):
  with open(hocFile, 'r') as infile:
    
    for line in infile:
      cols = line.strip().split('_')
      if cols[0] == 'connect section':
        continue
      
      elif cols[0] == 'create section':
        current_section = cols[1]
        S = Section(current_section)
        
        # start a new section
      
      elif cols[0] == 'section':
        continue
      
      elif cols[0] == '}':
        # end section
        current_section = -1
      
      else:
        pt3dcols = line.strip().split('(')
        if pt3dcols[0] == 'pt3dclear':
          continue
        elif pt3dcols[0] == 'pt3dadd' && current_section >= 0:
          # add new points
          
          S.addPoint
  












































#####################################################
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]

  getCarrierPoints(hocFile)


