# loadSTL.py -- import STL (ascii) file and save the mesh points in a 
#               in N x 4 array (x, y, z, area; N points total)
# usage: python loadSTL.py -options

import sys, os
import numpy as np
import zipfile
import math



class STL:
  def __init__(self, stlFileName, newFileName=False):
    self.inFile = stlFileName
    self.outFile = newFileName
    self.tempvertices = []
    self.vertices = []
    self.areas = []
    self.iszip = None
    self.archiveName = None
    self.properties = {}
    
    
    # run the program
    self.archive()
    self.readFile()
    self.writeData()
    self.checkHeader()
    
  
  
    
  def archive(self):
    # create name of inner unzipped file
    name = str(self.inFile)
    splitname = name.split('.')
    if splitname[-1] == 'zip':
      self.archiveName = '.'.join(splitname[:2])
      self.iszip = True
    else:
      self.archiveName = self.inFile
      self.iszip = False
    
    return self
  
  
  
  def computeVertex(self, dline):
    # create the vertex from the dataline when 3 points are received
    
    def sidelength(A,B):
      length = math.sqrt( (A[0]-B[0])**2 +
                          (A[1]-B[1])**2 +
                          (A[2]-B[2])**2 )
      return length
    
    def triarea(A,B,C):
      a, b, c = sidelength(B,C), sidelength(A,C), sidelength(A,B)
      p = 0.5 * (a+b+c)
      area = math.sqrt( p*(p-a)*(p-b)*(p-c) )
      return area
    
    self.tempvertices.append(dline)
    if len(self.tempvertices) == 3:
      meanX = np.mean([self.tempvertices[i][0] for i in xrange(3)])
      meanY = np.mean([self.tempvertices[i][1] for i in xrange(3)])
      meanZ = np.mean([self.tempvertices[i][2] for i in xrange(3)])
      self.vertices.append([meanX, meanY, meanZ])
      curarea = triarea(self.tempvertices[0], self.tempvertices[1],
                        self.tempvertices[2])
      self.areas.append(curarea)
      self.tempvertices = []
    
    return self
  
  
  
  def parseLine(self, splitLine):
    # parse each line of the unzipped file
    # calls computeVertex
    if not splitLine:
      pass
    
    elif splitLine[0] == 'solid':
      self.properties['triangles'] = int(splitLine[2].split('(')[1])
      self.properties['nodes'] = int(splitLine[4])
    
    elif splitLine[0] == 'vertex':
      dline = [float(i) for i in splitLine[1:]]
      self.computeVertex(dline)
    
    return self
    
  
  
  def readFile(self):
    # read the lines of the zipped file
    # calls parseLine, which calls computeVertex
    
    lineNum = 0
    
    if self.iszip == True:
      print('Unzipping to read file %s ...' % self.inFile)
      with zipfile.ZipFile(self.inFile) as z:
        print(' done. ')
        print('Reading file %s ...' % self.archiveName)
        with z.open(self.archiveName) as fIn:
          for line in fIn:
            
            splitLine = line.split(None)
            self.parseLine(splitLine)
            lineNum = lineNum + 1
            
            # print a progress every 10 million lines
            if lineNum % 10000000 == 0:
              print('%i lines read' % lineNum)
              
    elif self.iszip == False:
        print('Reading file %s ...' % self.archiveName)
        with open(self.archiveName) as fIn:
          for line in fIn:
            
            splitLine = line.split(None)
            self.parseLine(splitLine)
            lineNum = lineNum + 1
            
            # print a progress every 10 million lines
            if lineNum % 10000000 == 0:
              print('%i lines read' % lineNum)
      
    
    print(' done.')
    return self



  def writeData(self):
    # write the N x 4 array to a .txt doc
    if self.outFile == False:
      self.outFile = 'surface.txt'
    fOut = open(self.outFile, 'w')
    
    print('Writing text file of vertices and areas...')
    for pt in xrange(len(self.vertices)):
      strline = self.vertices[pt]
      strline.append(self.areas[pt])
      strline = [str(i) for i in strline]
      strline = ' '.join(strline)
      fOut.write(strline)
      fOut.write('\n')
    print(' done.')
    
    
  
  def checkHeader(self):
    # check number of triangles/nodes
    print('Number of triangles according to header: %i'
          % self.properties['triangles'])
    print('Number of triangles with areas computed: %i'
          % len(self.areas))
    print('Number of nodes according to header: %i' 
          % self.properties['nodes'])
    print('3 x num. of areas computed: %i' % (3*len(self.areas)))
    
    return self




############################
if __name__ == '__main__':
  arguments = sys.argv
  stlFile = arguments[1]
  if len(arguments) > 2:
    newFile = arguments[2]
  else:
    newFile = False
  centers = STL(stlFile, newFile)

















