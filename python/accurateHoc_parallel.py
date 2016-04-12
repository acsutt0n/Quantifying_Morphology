# accurateHoc.py -- this toolbox creates a very accurate hoc file
# usage: python accurateHoc.py hocFile tiffFolder outFile
# Input: 1. hocFile - hoc-format, really only used for the skeleton
#        2. tiffStack - the directory containing a tiffstack in B&W
#        3. outFile - filename for the completed hoc file
# Output: 1. outFile - a hoc file
# Dependencies: 
#        1. NeuronGeometry.py and neuron_readExportedGeometry.py, and
#           neuron_readExportedGeometry.demoread must return geometry
#        2. imageMatrix.py

import math
import numpy as np
from imageMatrix import *
from neuron_readExportedGeometry import *
from imageMatrix import *
from spiral import *
from timeit import default_timer as timer



class Hoc():
  def __init__(self, geoFile, tiffFolder, voxel=[0.176,0.176,0.38], threads=None):
    if tiffFolder == None:
      self.tiffFolder = '/home/alex/data/morphology/848/848_081/fake_filament/'
    else:
      self.tiffFolder = tiffFolder
    self.threads = threads
    self.geoFile = geoFile
    self.voxel = voxel
    self.geofile = geoFile
    self.varr = None
    self.load_geometry()
    self.skelpoints, self.skel_vectors, self.skel_names = [], [], []
    self.skel_coords() # load skelcoords and skel vectors
    self.segments = {'*':None} # segments and nodes will go here
    self.hoc = [] # tuple of [seg_name,[x,y,z],rad]
    self.get_cross_sections()
    self.get_connections()
    self.writeHoc()
    
    
    
  ##################### Loading functions ##################################
  
  def load_geometry(self):
    # load voxelized image from directory, discard coordinates
    _, self.varr = get_voxel_locations(self.tiffFolder, 'f', self.voxel, False)
    self.geometry = demoRead(self.geoFile)
    return self


  def skel_coords(self):
    # return all skelcoords
    for s in self.geometry.segments:
      for n in range(len(s.nodes)-1):
        pt0 = [s.nodes[n].x, s.nodes[n].y, s.nodes[n].z]
        pt1 = [s.nodes[n+1].x, s.nodes[n+1].y, s.nodes[n+1].z]
        self.skelpoints.append(pt0)
        self.skel_vectors.append(get_vector(pt0,pt1))
        self.skel_names.append(s.name)
        # this method skips the initial connecting plane, which
        # would be really messy anyway
    return self



  ######################## Small functions ################################

  def dist(pt0, pt1):
    if len(pt0) == len(pt1):
      return math.sqrt( (pt0[0]-pt1[0])**2 + 
                        (pt0[1]-pt1[1])**2 +
                        (pt0[2]-pt1[2])**2 )
    else:
      print('Dimension mismatch:')
      print(pt0, pt1)






  ####################### Money functions ##################################
  
  def add_node(self, node):
    # n = str(n)
    print('Writing all nodes to self.segments...')
    if node['segname'] not in self.segments.keys():
      # if it doesn't exist, then the first node doesn't exist either
      self.segments[node['segname']] = {'0': node}
    else:
      new = len(self.segments[node['segname']])
      self.segments[node['segname']][str(new)] = node
    return self

  
  def get_cross_sections(self):
    """
    Do like everything. Paralleize?
    """
    from multiprocessing import Pool
    if self.threads:
      pool = Pool(self.threads)
    else:
      pool = Pool()
    N = range(len(self.skelpoints)-1)
    # nodes = pool.map(self.par_cross_sections, N)
    # or:
    nodes = [self.par_cross_sections(i) for i in N]
    # doesn't seem to be any faster....
    #pool.close()
    #pool.join()
    for node in nodes:
      self.add_node(node)
    return self
    
  def par_cross_sections(self, n):
    plancoords, numpts = \
      gen_plane(self.skelpoints[n], self.skelpoints[n+1],self.voxel) # imageMatrix
    print('Called gen_plane')
    cross_dims = plane_XY(plancoords, numpts) # imageMatrix
    #print('Called plane_XY')
    #vcoords = scale_plane(plancoords, self.skelpoints[n], self.voxel) # imageMatrix
    #print('Called scale_plane')
    sarr, start_pos = return_cross_sec_array(plancoords, self.varr, 
                                             numpts, self.voxel) # imageMatrix
    print('Called return_cross_sec_array')
    if start_pos is not None:
      s = Spiral(sarr, start_pos)#, self.skelpoints[n]) # spiral
    else:
      
    print('Points for node %i done in parallel.' %n)
    node = {'coord': self.skelpoints[n], 
            'area': s.area*cross_dims[0]*cross_dims[1],
            'surface': s.surface['xs']*cross_dims[0] + \
                       s.surface['ys']*cross_dims[1],
            'arr': reproduce_matrix(s.live_pts),
            'segname': self.skel_names[n],
            'start': s.start}
    node['rad'] = np.sqrt(node['area']/np.pi)
    return node

    
    """
    # create cross-sections to populate nodes & segments
    for n in range(len(self.skelpoints)):
      plancoords, numpts = gen_plane(self.skel_vectors[n],self.voxel) # imageMatrix
      cross_dims = plane_XY(plancoords, numpts) # imageMatrix
      vcoords = scale_plane(plancoords, self.skelpoints[n], self.voxel) # imageMatrix
      sarr = return_cross_sec_array(vcoords, self.varr, numpts) # imageMatrix
      s = Spiral(sarr) # spiral
      # populate segments and nodes dicts
      self.add_node(n,s,cross_dims)
      #self.hoc.append([self.skel_names[n], self.skelpoints[n],
      #                 self.segments[self.skel_names[n]][str(len(self.segments[self.skel_names[n]]))]['rad']])
    # done
    return self
    """

  
  ################### Writing the Hoc file ###############################

  def get_connections(self):
    conns = []
    lineNum = 0
    
    with open(self.geofile, 'r') as fIn:
      for line in fIn:
        lineNum = lineNum + 1
        splitLine = line.split(None)
        try:
          if splitLine[0] == 'connect':
            # save everything 
            conns.append(line)
        except:
          pass
    # should have all the connections now as a list in conns
    return conns
        
  
  
  def writeHoc(self, newName='newhoc.hoc'):
    """
    Create a new hoc file, but mantain all of the previous connections
    info since that's a pain to reproduce.
    """
    def parse_segnum(sname):
      # Return the segment number given the segment name: filament_999[0]
      ss = sname.split('[')
      sss = ss[1].split(']')[0]
      return int(sss)
      
    conns = self.get_connections()
    with open(newName, 'w') as fOut: # MAKE SURE this file doesn't exist
      # since it is set to append?? (a instead of w)
      
      def create_filament():
        def get_name():
          for sname in self.segments.keys():
            if sname != '*':
              return sname.split('[')[0]
        s = get_name()+'['+str(len(self.segments.keys())-1)+']'
        s = 'create ' + s
        fOut.write(s)
        fOut.write('\n\n')
        return
      
      def new_filament(name):
        s = name + ' {'
        fOut.write(s)
        fOut.write('\n')
        fOut.write('  pt3dclear()')
        fOut.write('\n')
        return
      
      def pt3dadd(coord, rad):
        s = ','.join([str(c) for c in coord])
        s = '  pt3dadd(' + s + ',' + str(rad) + ')'
        fOut.write(s)
        fOut.write('\n')
        return
      
      def end_filament():
        fOut.write('}')
        fOut.write('\n')
        return
      
      def add_last_node(sname, rad):
        nod = self.geometry.segments[parse_segnum(sname)].nodes[-1]
        pt3dadd([nod.x, nod.y, nod.z], rad)
        return
      
      # start writing
      create_filament()
      for sname in self.segments.keys():  # for each segment
        if sname != '*':
          new_filament(sname)               # start a new segment
          for node in range(len(self.segments[sname].keys())): # for each node
            coord = self.segments[sname][str(node)]['coord']
            rad = self.segments[sname][str(node)]['rad'] 
            pt3dadd(coord, rad)  # write that node + radius
          add_last_node(sname, rad) # add the last node
          end_filament()                    # end that filament
      
      # write the connections
      for conline in conns:
        fOut.write(conline)
    # newhoc now closed
    return
    
    
        
    
    
     
    
    
  

"""
class Seg:
  def __init__(self, 
"""







