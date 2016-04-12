# hocTortuosity.py -- this program calculates the tortuosity
#    of a given branch/segment
#    tortuosity = euclidean distance / path length



import os
import scipy
from scipy import special
import matplotlib.pyplot as pyplot



def getTortuosity(hocFile):
  

  ### open file for writing
  outfile = open(tortFile, 'w')
    

  class Section:
    def __init__(self):
      self.name = None
      self.nodes = {}
      self.length = 0
      self.pathdist = 0
      #self.surfaceArea = 0
      #self.volume = 0
      self.tort = 0
      self.centroid = 0
    
    def getLength(self):
      zeroNode = self.nodes[0]
      lastnode = len(self.nodes) - 1
      oneNode = self.nodes[lastnode]
      self.length = scipy.sqrt( (zeroNode[0] - oneNode[0])**2 + \
                                (zeroNode[1] - oneNode[1])**2 + \
                                (zeroNode[2] - oneNode[2])**2 )
      return self.length
    
    def getPathDist(self):
      total_dist = 0
      numNodes = len(self.nodes)
      for n in range(numNodes - 1):
        total_dist = total_dist + pointDist(self.nodes[n], self.nodes[n+1])
      self.pathdist = total_dist
      return self.pathdist
    
    def pointDist(pt0, pt1):
      dist = scipy.sqrt( (pt0[0] - pt1[0])**2 + (pt0[1] - pt1[1])**2 + \
                         (pt0[2] - pt1[2])**2 )
      return dist
  
    def findTort(self):
      self.tort = self.length / self.pathdist
      return self.tort
    
    def getCentroid(self):
      zeroNode = self.nodes[0]
      oneNode = self.nodes[-1]
      self.centroid = [ (zeroNode[0] + oneNode[0])/2, \
                        (zeroNode[1] + oneNode[1])/2, \
                        (zeroNode[2] + oneNode[2])/2 ]
      return self.centroid
    
    
    ### write outFile
    def writeTort(self):
      tortline = [self.name, self.length, self.pathdist, self.tort, \
                  self.centroid[0], self.centroid[1], self.centroid[2]]
      newline = ['0','0','0','0','0','0','0']
      for c in range(len(newline)):
        newline[c] = str(tortline[c])
      tortString = ' '.join(newline)
      outfile.write(newline)
      outfile.write('\n')


  
  ######## parseHocSegments
  def parseHocSegments(hocFile):
    
    openSection = 0
    with open(hocFile, 'r') as infile:
      
      
      for line in infile:
        cols = line.strip().split('_') 
        if cols[0] == 'connect section':
          continue # skip this part now
      
        elif cols[0] == 'section' and openSection == 0:
          newcols = cols[1].split('{')
          current_section = int(newcols[0])
          openSection = 1
          newSec = Section()
          newSec.name = current_section
          # print('found section')
        
        elif cols[0] == 'section' and openSection == 1:
          newcols = cols[1].split('{')
          current_section = int(newcols[0])
          print('error trying to create section while another section is open')
          print('section %i' % current_section)
          
        # if end of a section is found, write section data
        elif cols[0] == '}' and openSection:
          openSection = 0
          newSec.length = newSec.getLength()
          # omit for now
          #newSec.surfaceArea = getArea(newSec.nodes)
          #newSec.volume = getVolume(newSec.nodes)
          newSec.pathdist = newSec.getPathDist()
          
          newSec.writeTort()
          print('new section %i written' % int(newSec.name))
          
          
        
        # parse lines to get points
        else:
          pt3dcols = cols[0].split('(')
          if pt3dcols[0] == 'pt3dclear': # skip it
            # print('pt3dclear skipped')
            continue
            
          # gotta add these points and rads if openSection==1
          elif pt3dcols[0] == 'pt3dadd' and openSection: 
            pointcols = pt3dcols[1].split(',')
            # into ['309.060028', '174.105026', '52.000000', '0.090000)']
            # coordinates
            
            X = float(pointcols[0])
            Y = float(pointcols[1])
            Z = float(pointcols[2])
            radcols = pointcols[3].split(')')
            rad = float(radcols[0])
            
            new_node = [X,Y,Z,rad]
            
            if newSec.nodes == {}:
              newSec.nodes = {0: new_node}
            else:
              numNodes = len(newSec.nodes)
              newSec.nodes[numNodes] = new_node

            # print('point added')
            
          elif pt3dcols[0] == '\n':
            continue
          else:
            continue # print('line not recognized')


  parseHocSegments(hocFile)
























#####################################################
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]
  tortFile = sys.argv[2]

  getTortuosity(hocFile)
