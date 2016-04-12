# loadHoc.py -- loads in hoc file 
#

import sys, os
import numpy as np



class HocGeometry:
  def __init__(self, fileName, newFile=False):
    self.fileName = fileName
    self.sections = {}
    self.connections = []
    self.openSection = None
    self.currentSection = None
    self.secRads = {}
    self.outFile = newFile
    self.medrad = False
    self.nodes = {}
    
    # run the code
    self.readFile()
    # self.matchRadius()
    self.showInfo()

    
  
  
  def addNode(self, splitLine):
    # add a node to the current section
    #ptsplit = splitLine.split(',')
    #print(ptsplit)
    pt1 = float(splitLine[0].split('(')[1].split(',')[0])
    pt2 = float(splitLine[1].split(',')[0])
    pt3 = float(splitLine[2].split(',')[0])
    rad = float(splitLine[3].split(')')[0])
    current_node = [pt1, pt2, pt3]
    if self.currentSection not in self.sections.keys():
      self.sections[self.currentSection] = [current_node]
    else:
      self.sections[self.currentSection].append(current_node)
    if self.currentSection not in self.secRads.keys():
      self.secRads[self.currentSection] = [rad]
    else:
      self.secRads[self.currentSection].append(rad)
    # add node to log
    self.nodes[len(self.nodes)] = [pt1, pt2, pt3, rad]
    
    return self
  
  
  
  def addConnection(self, splitLine):
    # add a connection item
    prevsec = int(splitLine[1].split('[')[1].split(']')[0])
    nextsec = int(splitLine[2].split('[')[1].split(']')[0])
    connect = [prevsec, nextsec]
    self.connections.append(connect)
    
    return self
  
  
  
  def readFile(self):
    # read in file
    
    print('Reading %s ...' %self.fileName)
    with open(self.fileName, 'r') as fIn:
      
      lineNum = 0
      for line in fIn:
        lineNum = lineNum + 1
        splitLine = line.split(None)
        
        if splitLine:

          if splitLine[0].split('_')[0] == 'filament':
            print('Found new segment')
            self.openSection = True
          elif splitLine[0] == 'connect':
            print('Found new connection entry')
            self.addConnection(splitLine)
          elif splitLine[0] == '}':
            # print('Added section %i' % self.currentSection)
            self.openSection = False
          elif self.openSection == True:
            if len(splitLine) > 1 and splitLine[1] == '{':
              self.currentSection = int(splitLine[0].split('[')[1].split(']')[0])
            elif splitLine[0].split('(')[0] == 'pt3dadd':
              print('Found new point')
              self.addNode(splitLine)
    print('Finished reading %i lines.' % lineNum)
    
    return self
  

  
  def showInfo(self):
    # display some basic info
    print('Number of sections: %i' % len(self.sections))
    rads, nodes = [], 0
    for sec in self.sections.keys():
      for n in self.sections[sec]:
        nodes = nodes+1
    print('Number of nodes: %i' % nodes)
    for sec in self.secRads.keys():
      for r in self.secRads[sec]:
        rads.append(r)
    if not self.medrad:
      self.medrad = rads[int(len(rads)/2)]
    print('Number of radii: %i' % len(rads))
    print('Median radius: %.5f , min radius: %.5f, max radius: % .5f' \
          % (self.medrad, min(rads), max(rads)) )
    # print(np.shape(self.connections))

    
    return self
  
  
  
  def whereConnect(self, refsec, fixsec):
    """
    For each entry in the connection matrix, examine the points that
    are supposed to be connected and if they don't match change the
    connection matrix to reflect the actual connection order.
    """
    ref10 = self.sections[refsec][0]
    ref11 = self.sections[refsec][-1]
    fix00 = self.sections[fixsec][0]
    fix01 = self.sections[fixsec][-1]
    if ref10 == fix00:
      return 0 # use the 0th node's rad of ref sec for new rad
    elif ref10 == fix01:
      return 0
    elif ref11 == fix00:
      return 1
    elif ref11 == fix01:
      return 1
    else:
      print('No valid connection found between sections %i and %i' 
            % (refsec, fixsec))
  

    
  def matchRadius(self):
    """
    If a radius of < 0 is found, its neighbors are used to get the
    correct radius
    """
    print('Making sure all points and radii match...')
    rads = []
    for sec in self.secRads.keys():
      for r in self.secRads[sec]:
        rads.append(r)
    if not self.medrad:
      self.medrad = rads[int(len(rads)/2)]
    print('Median radius is: %.5f' % self.medrad)
    
    self.uniqueNodes = []
    self.uniqueRads = []
    
    for sec in self.sections.keys():
      for n in xrange(len(self.sections[sec])):
        
        # check to see if that point already exists in uniquesecs
        if self.sections[sec][n] in self.uniqueNodes:
          radInd = self.uniqueNodes.index(self.sections[sec][n])
          self.secRads[sec][n] = self.uniqueRads[radInd]
        else:
          self.uniqueNodes.append(self.sections[sec][n])
          if self.secRads[sec][n] <= 0 or self.secRads[sec][n] > 10000:
            print('Bad radius found, section %i node %i' %(sec, n))
            self.uniqueRads.append(self.medrad)
            self.secRads[sec][n] = self.medrad
            print('Replaced with: %.5f' % self.uniqueRads[-1])
          else:
            self.uniqueRads.append(self.secRads[sec][n])

    print('Radii fixed.')
    
    return self
    
    
  
  
  def scanFixRadius(self):
    # fix radii that don't match up
    for sec in self.sections.keys():
      curr_rad0, curr_rad1 = self.secRads[sec][0], self.secRads[sec][-1]
      
      if curr_rad0 > 1000 or curr_rad0 <= 0:
        print('Section %i has a (0) radius of %.5f' % (sec, curr_rad0))
        # find neighbors
        nebs0, rads0 = [], []
        nebs0 = [self.connections[i][0] \
                 for i in xrange(len(self.connections)) \
                 if self.connections[i][1] == sec]
        rads0 = [self.secRads[i][-1] for i in nebs0]
        # print(rads0, nebs0)
        rad_here = self.whereConnect(nebs0[0],sec)
        self.secRads[sec][0] = self.secRads[nebs0[0]][rad_here] #[nebs0.index(neb_min)]
        print('New rad is %.5f' %self.secRads[sec][0])
      
      if curr_rad1 > 1000 or curr_rad1 <= 0:
        print('Section %i has a (-1) radius of %.5f' % (sec, curr_rad1))
        # find neighbors
        nebs1, rads1 = [], []
        nebs1 = [self.connections[i][1] \
                 for i in xrange(len(self.connections)) \
                 if self.connections[i][0] == sec]
        rads1 = [self.secRads[i][0] for i in nebs1]
        # print(rads1, nebs1)
        rad_here = self.whereConnect(nebs1[0], sec)
        self.secRads[sec][-1] = self.secRads[nebs1[0]][rad_here] #rads1[nebs1.index(neb_min)]
        print('New rad is %.5f' %self.secRads[sec][-1])
      
    return self
      
      
    """
    nebs1 = [self.connections[i][1] \
             for i in xrange(len(self.connections)) \
             if self.connections[i][0] == sec]
    # and the associated radii

    rads1 = [self.secRads[i][0] for i in nebs1]
    
    # check radii
    
    chk0, chk1 = 0, 0
    for r in xrange(len(rads0)):
      if rads0[r] != curr_rad0:
        print('sections %i, %i, radii %.6f, %.6f' \
              % (sec, nebs0[r], curr_rad0, rads0[r]))
    for r in xrange(len(rads1)):
      if rads1[r] != curr_rad1:
        print('sections %i, %i, radii %.6f, %.6f' \
              % (sec, nebs1[r], curr_rad1, rads1[r]))
  """
  
  

  def fixRadius(self):
    # correct the wrong radii using it's neighbors or the median radius
    # for each segment, find its neighbors
    rads = []
    for sec in self.secRads.keys():
      for r in self.secRads[sec]:
        rads.append(r)
    rads = np.array(rads)
    #print(np.shape(rads))
    medrad = rads[int(len(rads)/2)]
    
    print('Fixing radii ...')
    for sec in self.secRads.keys():
      if np.mean(self.secRads[sec]) > 1000*medrad:
        # find neighbors at each end
        nebs0, nebs1 = [], []
        print('Bad radius of section %i found, mean: %.5f' \
              % (sec, np.mean(self.secRads[sec])) )
        nebs0 = [self.connections[i][0] \
                 for i in xrange(len(self.connections)) \
                 if self.connections[i][1] == sec]
        nebs1 = [self.connections[i][1] \
                 for i in xrange(len(self.connections)) \
                 if self.connections[i][0] == sec]
        # want the -1th rad for 1-connects and the 0th rad for 0-connects
        rads0 = [self.secRads[i][-1] for i in nebs0]
        rads1 = [self.secRads[i][0] for i in nebs1]
        # print('Sec %i has %i neighbors' % (sec, (len(nebs0)+len(nebs1))))
        
        # if the mean endpoint radius is acceptable, use it
        if np.mean(rads0) > 1000*medrad:
          # otherwise use the min radius (if it's acceptable)
          if min(rads0) > 1000*medrad or min(rads1) < 0:
            # otherwise is the median
            self.secRads[sec][0] = medrad
          else:
            self.secRads[sec][0] = min(rads0)
        else:
          self.secRads[sec][0] = np.mean(rads0)
        # repeat for rads1
        if np.mean(rads1) > 1000*medrad:
          if min(rads1) > 1000*medrad or min(rads1) < 0:
            self.secRads[sec][-1] = medrad
          else:
            self.secRads[sec][-1] = min(rads1)
        else:
          self.secRads[sec][-1] = np.mean(rads1)
        print('New mean: %.5f' % np.mean(self.secRads[sec]))
    
      for node in xrange(len(self.secRads[sec])):
        if self.secRads[sec][node] > 1000*medrad or \
           self.secRads[sec][node] <= 0:
             print('Bad radius of section %i found, radius: %.5f' \
                   % (sec, self.secRads[sec][node]))
             if np.mean(self.secRads[sec]) < 1000*medrad and \
                np.mean(self.secRads[sec]) > 0:
                  self.secRads[sec][node] = np.mean(self.secRads[sec])
             else:
               self.secRads[sec][node] = medrad
             print('New radius: %.5f' % self.secRads[sec][node])
    
    print('Median radius is: %.5f' % medrad)
    return self
  
  
  
  def medianDist(self, secNum):
    # return the median distance between points of a given section
    distlist = []
    count = 0
    # print('Length of section: %i' % len(self.sections[secNum]))
    for node in xrange(len(self.sections[secNum])-1):
      d = [a - b for a,b in zip(self.sections[secNum][count], \
                                self.sections[secNum][count+1])]
      distlist.append(np.linalg.norm(d))
      count = count + 1
    distlist.sort()
    meddist = distlist[int(len(distlist)/2)]
    
    return meddist
  
  
  
  def findLongSections(self, version=2):
    """
    Determines whether points in certain segments are spaced too far
    apart (relative to median point spacing). 
    """
    median_distances = [self.medianDist(sec) for sec in self.sections.keys()]
    median_distances.sort()
    median_dist = median_distances[int(len(median_distances)/2)]
    long_sections, long_distances = [], []
    
    if version == 1:
      for s in xrange(len(median_distances)):
        if median_distances[s] > 2*median_dist:
          long_sections.append(median_distances[s])
          long_distances.append(s)
      print('Found %i sections with points spaced far apart' \
            % len(long_sections))
    elif version == 2:
      for sec in self.sections.keys():
        if len(self.sections[sec]) < 3:
          long_sections.append(sec)
          d = [a - b for a,b in zip(self.sections[sec][0], \
                                    self.sections[sec][-1])]
          long_distances.append(np.linalg.norm(d))
    
    print(np.shape(long_distances))
    # sections_distances = zip(long_sections, long_distances)
    return long_sections, long_distances, median_dist
  
  
  
  def interpPoints(self, interpRad=False):
    """
    Interpolate the points and radii between sections that have too
    few points.
    """
    # print(np.shape(long_distances))
    long_sections, long_distances, meddist = self.findLongSections()
   
    print('Long inter-point distances found: %i' % len(long_sections))
    count = 0
    for sec in long_sections:
      print('Supposed long section %i has %i nodes' \
            % (sec, len(self.sections[sec])))
      # set first and last points for interpolation
      pt0, pt1 = self.sections[sec][0], self.sections[sec][-1]
      # find number of points
      numpts = int(long_distances[long_sections.index(sec)]/meddist)
      Xs = np.linspace(pt0[0], pt1[0], numpts)
      Ys = np.linspace(pt0[1], pt1[1], numpts)
      Zs = np.linspace(pt0[2], pt1[2], numpts)
      
      newpts = np.dstack((Xs, Ys, Zs))
      newpts = [newpts[0][i] for i in xrange(len(newpts[0]))]
      self.sections[sec] = newpts
      count = count + 1
      rad0, rad1 = self.secRads[sec][0], self.secRads[sec][-1]
      # print(rad0, rad1)
      rads = np.linspace(rad0, rad1, numpts)
      # print(rads)
      self.secRads[sec] = rads 
      
    long_sections, long_distances, meddist = self.findLongSections()
    print('Long sections still remaining: %i' % len(long_sections))
    if len(long_sections) > 0:
      print(long_distances, meddist)
    
    return self
    
    
    
  def writeHoc(self):
    """
    Write a hoc file.
    """
    print('Writing output file %s ...' % self.outFile)
    with open(self.outFile, 'w') as fOut:
      
      def createSection(secNum):
        fOut.write('create section_%i\n' %secNum)
        fOut.write('section_%i {\n' %secNum)
        fOut.write('pt3dclear()\n')
        for node in xrange(len(self.sections[secNum])):
          fOut.write('pt3dadd(%.6f, %.6f, %.6f, %.6f)\n' \
                     % (self.sections[secNum][node][0],
                        self.sections[secNum][node][1],
                        self.sections[secNum][node][2],
                        self.secRads[secNum][node]))
        fOut.write('}\n')
      
      def createConnection():
        for c in xrange(len(self.connections)):
          fOut.write('connect section_%i(1), section_%i(0)\n' \
                     % (self.connections[c][0],self.connections[c][1]))
      
      
      for sec in self.sections.keys():
        createSection(sec)
      createConnection()
    
    
    return 
    



##################################
if __name__ == '__main__':
  arguments = sys.argv
  hocFile = arguments[1]
  if len(arguments) > 2:
    newFile = arguments[2]
    geometry = HocGeometry(hocFile, newFile)
  else:
    geometry = HocGeometry(hocFile)
    
########################################################################
def writeHoc(geo, fName):
  """
  Write a hoc file.
  """
  print('Writing output file %s ...' % fName)
  with open(fName, 'w') as fOut:
    
    def createSection(geo, secNum):
      fOut.write('create section_%i\n' %secNum)
      fOut.write('section_%i {\n' %secNum)
      fOut.write('pt3dclear()\n')
      for node in range(len(geo.segments[secNum].nodes)):
        fOut.write('pt3dadd(%.6f, %.6f, %.6f, %.6f)\n' \
                   % (geo.segments[secNum].nodes[node].x,
                      geo.segments[secNum].nodes[node].y,
                      geo.segments[secNum].nodes[node].z,
                      geo.segments[secNum].nodes[node].r1))
      fOut.write('}\n')
    
    def createConnection(geo):
      for c in range(len(geo.connections)):
        con_list = [geo.connections[c]['filament1'],
                    geo.connections[c]['location1'],
                    geo.connections[c]['filament2'],
                    geo.connections[c]['location2']]
        fOut.write('connect %s(%.1f), %s(%.1f)\n' \
                   % (con_list[0], con_list[1], con_list[2], con_list[3]))
    
    
    for sec in range(len(geo.segments)):
      createSection(geo, sec)
    createConnection(geo)

