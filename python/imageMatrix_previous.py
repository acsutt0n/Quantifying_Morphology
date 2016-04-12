# imageMatrix.py - image processing toolbox for neuron stuff
# usage: python imageMatrix.py tiffDirectory/ outFile

import os, sys, math
import numpy as np
from PIL import Image
from multiprocessing import Pool
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from timeit import default_timer as timer
from spiral import *
#import pickle


# SIMPLE TOOLS
########################################################################
def load_img_array(imFile):
  img = Image.open(imFile)
  arr = np.array(img)
  return arr


def gen_vec(pt0, pt1):
  if len(pt0) == len(pt1):
    return [pt1[i]-pt0[i] for i in range(len(pt0))]
  else:
    print('Dimension mismatch: %i and %i' %(len(pt0),len(pt1)))


def coord_bounds(coords):
  if np.shape(coords)[1] != 3 and np.shape(coords)[1] != 2:
    print('Coords must be N x 3 or N x 2 array')
    return np.shape(coords)
  if np.shape(coords)[1] == 3:
    cx,cy,cz=[c[0] for c in coords],[c[1] for c in coords],[c[2] for c in coords]
    bounds = [min(cx), max(cx), min(cy), max(cy), min(cz), max(cz)]
    print('Xmin: %.2f, Xmax: %.2f, Ymin: %.2f, Ymax: %.2f, Zmin: %.2f, Zmax %.2f'
          %(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]))
  else:
    cx,cy=[c[0] for c in coords],[c[1] for c in coords]
    bounds = [min(cx), max(cx), min(cy), max(cy)]
    print('Xmin: %.2f, Xmax: %.2f, Ymin: %.2f, Ymax: %.2f'
          %(bounds[0],bounds[1],bounds[2],bounds[3]))
  return bounds


def show_matrix(arr):
  show_array(arr)
def show_array(arr):
  img = Image.fromarray(arr)
  img.show()
  return


def show_dist(l, dataPerBin=100):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.hist(l, bins=int(len(l)/dataPerBin))
  plt.show()
  

def make255(arr, T=255):
  zarr = np.zeros(np.shape(arr))
  for i in range(len(arr)):
    for j in range(len(arr[0])):
      if arr[i,j] > 0:
        zarr[i,j] = T
  return zarr


def dist(pt0,pt1):
  return np.sqrt(sum([(pt1[i]-pt0[i])**2 for i in range(len(pt0))]))


# invert and flatten
def invert(arr):
  K = np.shape(arr)
  if len(K)==3:
    zarr = np.zeros(np.shape(arr)[:2])
    for i in range(np.shape(arr)[0]):
      for j in range(np.shape(arr)[1]):
        if arr[i,j,0] == 255:
          zarr[i,j] = 0
        if arr[i,j,0] == 0:
          zarr[i,j] = 255
  elif len(K)==2:
    zarr = np.zeros(np.shape(arr)[:2])
    for i in range(np.shape(arr)[0]):
      for j in range(np.shape(arr)[1]):
        if arr[i,j] == 255:
          zarr[i,j] = 0
        if arr[i,j] == 0:
          zarr[i,j] = 255
  else:
    print('Bad shape of array sent to invert')
  return zarr


def matrix2coords(darr, voxel):
  # A 3-D array is converted to a list of triplets at each voxel
  print('Converting the matrix to coordinates.')
  dims = len(voxel)
  coords = []
  if dims == 3:
    for z in range(len(darr)):
      for x in range(len(darr[z])):
        for y in range(len(darr[z][x])):
          if darr[z][x][y] > 0:
            coords.append([x*voxel[0],y*voxel[1],z*voxel[2]])
      print('z-slice # %i done' %int(z+1))
  elif dims == 2:
    for x in range(len(darr)):
      for y in range(len(darr[x])):
        if darr[x][y] > 0:
          coords.append([x*voxel[0],y*voxel[1]])
  return coords


def coords2matrix(coords, voxel, i=1):
  # x-y-z tuples are converted to a matrix; "i" is the matrix value
  matrix = []
  dims = len(voxel)
  for c in coords:
    matrix.append([int(c[i]/voxel[i]) for i in range(dims)])
  matrix = np.array(matrix)
  maxs = [max(matrix[:,i]) for i in range(dims)]
  print(maxs, dims)
  newmat = np.zeros([i+1 for i in maxs])
  if dims == 2:
    for m in matrix:
      newmat[m[0],m[1]] = i
  elif dims == 3:
    for m in matrix:
      newmat[m[0],m[1],m[2]] = i
  return newmat


def get_vector(pt0, pt1):
  if len(pt0) == len(pt1):
    return [x-y for x,y in zip(pt0,pt1)]



def test_plane(plan, cross_sec, voxel):
  # if plane isn't large enough (all outer 1's adjacent to 0's) it doubles the size
  if len(plan[0]) != 3 or len(voxel) != 3 or len(cross_sec[0]) != 3:
    print('Need 3 values (x,y,z) for plane, cross section and voxel args')
    return False
  
  def distance(pt0,pt1):
    return np.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
  
  def plane_info(plan):
    def farthest_pt(pt,pts):
      long_dist = 0
      for p in pts:
        if distance(p,pt) > long_dist:
          long_dist = distance(p,pt)
      return long_dist
      
    # start here
  
  def check_cs(p_info, cs_info):
    # if the span of plan <= span of cs, the plane needs to be extended
    # so that it is larger than the cross-section
    expand = 1
    for i in range(3):
      if p_info['M'][i] <= cs_info['M'][i]:
        expand = expand*0 
    return [False, True][expand]
  # if even 1 dimension is identical the plane should be extended
  
  #def expand_plan(p_info, voxel):
  #  minmax = [['xmin','xmax'],['ymin','ymax'],['zmin','zmax']]

  p_info = plan_info(plan)
  cs_info = plan_info(cross_sec)
  expand = check_cs(p_info, cs_info)
  return expand
  #if expand:
  #  new_plane = expand_plan(p_info, voxel)
  #  return [expand, new_plane]
  #else:
  #  return [expand]





def switch_binary(arr):
  zarr = np.zeros(np.shape(arr))
  for i in range(len(arr)):
    for j in range(len(arr[0])):
      if arr[i,j] == 1:
        zarr[i,j] = 125
      elif arr[i,j] == 2:
        zarr[i,j] = 255
  return zarr


def clean_filament(darr):
  """ Prints the locations of the skelpoints in the fake filament. """
  #carr = np.zeros(np.shape(darr))
  skelpts = []
  if len(np.shape(darr)) == 3:
    for s in range(len(darr)): # for each slice
      width, start, row  = 0, None, None
      for r in range(len(darr[s])): # for each row
        count = 0
        for c in range(len(darr[s][r])):
          if darr[s][r][c] > 0 and count == 0:
            temp_start = c
          if darr[s][r][c] > 0:
            count = count + 1
        # at the end of each row:
        if count > width:
          width = count
          start = temp_start
          row = r
      # at the end of each slice:
      # skelpoint location: (slice #, row #, start + k)
      print(s, row, start, int(width/2))
      skelpts.append([s, row, start+int(width/2)])
      darr[s][row][start+int(width/2)] = 2
  if len(np.shape(darr)) == 2:
    width, start, row  = 0, None, None
    for r in range(len(darr)): # for each row
      count = 0
      for c in range(len(darr[r])):
        if darr[r][c] > 0 and count == 0:
          temp_start = c
        if darr[r][c] > 0:
          count = count + 1
      # at the end of each row:
      if count > width:
        width = count
        start = temp_start
        row = r
    # at the end of each slice:
    # skelpoint location: (slice #, row #, start + k)
    print(row, start+ int(width/2)) 
    skelpts.append([row, start+int(width/2)])
    darr[row][start+int(width/2)] = 2
  return darr, skelpts



def plot_cross_secs(cs):
  # assumes each element in cs is a cross section with a [x,y,z]
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  cols = ['b','r','g','y','k']*int(len(cs)/5+1)
  for c in range(len(cs)):
    for i in range(len(cs[c])):
      ax.scatter(cs[c][i][0],cs[c][i][1],cs[c][i][2], c=cols[c], 
                 edgecolor = cols[c], alpha = 0.1)
  plt.show()



#########################################################################
######################### fake cross sections ###########################
"""
def fake_cross_section(fake_filament='/home/alex/data/morphology/848/848_081/fake_filament/'):
  def dist(inds0, inds1):
    if len(inds0) != len(inds1): print('Dimension mismatch!')
    else:
      return sum([math.sqrt((x-y)**2) for x,y in zip(inds0, inds1)])
  
  def return_nearest_dist(pt, points):
    mindist = np.inf
    for p in points:
      currdist = dist(pt, p)
      if currdist < mindist:
        mindist = currdist
    return mindist
  
  def get_skelcoords(darr, voxel, skel=2):
    skelcoords = []
    for s in range(len(darr)):
      for i in range(len(darr[s])):
        for j in range(len(darr[s][i])):
          if darr[s][i][j] == 2:
            skelcoords.append([x*y for x,y in zip(voxel,[i,j,s])])
    return skelcoords
  
  def return_cross_section(pt0, pt1, darr, voxel, M=10):
    # so that the plane can be scaled, get two points and instead of vec
    vec = get_vector(pt0, pt1)
    plancoords = gen_plane(vec, voxel, M, False)
    # scale the plancoords to match with pt0
    for p in range(len(plancoords)):
      for i in range(3):
        plancoords[p][i] = plancoords[p][i] + pt0[i]
    voxelcoords = matrix2coords(darr, voxel)
    vdist = 2 * dist([0,0,0],voxel)
    cross_sec = []
    for p in plancoords:
      if return_nearest_dist(p, voxelcoords) <= vdist:
        cross_sec.append(p)
        
    return cross_sec, plancoords
    # plan_mat = coords2matrix(plancoords, voxel)
    # want to find overlap between plane and segment
    
  def get_vector(pt0, pt1):
    if len(pt0) == len(pt1):
      return [x-y for x,y in zip(pt0,pt1)]
    
  voxel = [0.0176,0.0176,0.38]
  darr = gen_segment(fake_filament) # raw image
  carr = clean_filament(darr)
  skelcoords = get_skelcoords(carr, voxel)
  cross_secs = []
  
  for s in range(len(skelcoords)-1):
    pt0, pt1 = skelcoords[s],skelcoords[s+1]
    M=100
    cs, plancoords = return_cross_section(pt0, pt1, carr, voxel, M)
    # notbigenough = test_plane(plancoords, cs, voxel)
    cross_secs.append(cs) # once the plane exceeds the cross section, append
    
  return cross_secs, carr
"""

def display_fake_filament(darr):
  rimgs = [Image.fromarray(darr[i]) for i in range(10)]
  for k in range(10):
    rimgs[k].show()
  return


def gen_segment(fake_filament):
  # default value
  
  def make_1d(arr, invert=True, gray=False): # make 1-D and grayscale, invert
    narr = np.zeros([len(arr), len(arr[0])])
    if gray:
      newscheme = [0, 125, 255] # gray skeleton on black neuron, background white
    else:
      newscheme = [0,2,1] # 'binary', not gray
    for i in range(len(arr)):
      for j in range(len(arr[0])):
        if arr[i,j,0] == 255:
          narr[i,j] = newscheme[0]
        elif arr[i,j,0] == 0:
          narr[i,j] = newscheme[2]
        else:
          narr[i,j] = newscheme[2]
    return narr
    
  def load_fake_segment(imfile):
    img = Image.open(imfile)
    arr = np.array(img)
    return make_1d(arr) # already been done!
    
  fils = os.listdir(fake_filament)
  fils = [fake_filament+i for i in fils]
  fils.sort()
  darr = []
  for f in fils:
    darr.append(load_fake_segment(f))
  #pool = Pool() # these calls were yielding a pickle error related to pool properties
  #darr = pool.map(load_fake_segment, fils) # 
  #pool.close()
  #pool.join()
  return darr # return multi-dimensional array


############### New fake segment creation with numpy ##################

def create_segment(dims=[96,96,96]): # June 2015
  # Create a fake segment with voxels
  stack = [np.zeros([dims[0],dims[1]]) for i in range(dims[2])]
  # process is an 8-by-8 square by default
  for z in range(96):
    for i in range(8):
      for j in range(8):
        try:
          stack[z][z+i][z+j] = 1
        except:
          pass
  # done with fake filament; now make fake skeleton
  skel = []
  for i in range(4,92):
    skel.append([i,i,i])
  return stack, skel


def show_slice(stack, z=int(len(stack)/2)):
  # Show a particular slice of interest
  for i in range(len(stack)):
    for j in range(len(stack[i])):
      for k in range(len(stack[i][j])):
        if stack[i][j][k] > 0:
          stack[i][j][k] = 255
  img = Image.fromarray(stack[z])
  img.show()
  return



def plot_multi_coords(skelcoords, voxelcoords=None, plancoords=None, alfs=None):
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  colors = ['r','b','k']
  if alfs is None:
    alphas =  [0.1,0.01,0.01]
  elif type(alfs) is int:
    alphas=[alfs, alfs, alfs]
  elif type(alfs) is list:
    alphas = alfs
  print('Starting plot...')
  for s in skelcoords:
    ax.scatter(s[0],s[1],s[2],c=colors[0],edgecolor=colors[0],alpha=alphas[0])
  if voxelcoords is not None:
    for v in voxelcoords:
      ax.scatter(v[0],v[1],v[2],c=colors[1],edgecolor=colors[1],alpha=alphas[1])
  if plancoords is not None:
    for p in plancoords:
      ax.scatter(p[0],p[1],p[2],c=colors[2],edgecolor=colors[2],alpha=alphas[2])
  plt.show()
  return
  



# BINARY OPTIONS 
########################################################################

def make_binary(imFile):
  # make a binary image with a moving average
  # this version is currently not being used
  arr = load_img_array(imFile)
  zarr = np.zeros(np.shape(arr))
  on = False
  avg = None
  firstcol, lastcol = 0, len(arr[0])-1
  for row in range(len(arr)):
    brightest = max(arr[row])
    ind = list(arr[row]).index(brightest)
    zarr[row,ind] = 1
    
    # proceed bidirectionally
    if ind != firstcol and ind != lastcol:
      # negative direction
      prev = brightest
      neg_ind = ind
      while neg_ind > firstcol:
        neg_ind = neg_ind - 1
        if abs((arr[row, neg_ind]-arr[row, ind]) \
          / max(arr[row,ind],arr[row,neg_ind])) < 0.1:
            # include in binary array
            zarr[row, neg_ind] = 1
            prev = arr[row, neg_ind]
      # positive direction
      prev = brightest
      pos_ind = ind
      while pos_ind < lastcol:
        pos_ind = pos_ind + 1
        if abs((arr[row, pos_ind]-arr[row,ind]) \
          / max(arr[row,ind], arr[row,pos_ind])) < 0.1:
            # include in binary array
            zarr[row, pos_ind] = 1
            prev = arr[row, pos_ind]
      # else they will be zero as pre-set
    elif ind == firstcol:
      # only need positive direction
      prev = brightest
      pos_ind = ind
      while pos_ind < lastcol:
        pos_ind = pos_ind + 1
        if abs((arr[row, pos_ind]-arr[row,ind]) \
          / max(arr[row,ind], arr[row,pos_ind])) < 0.1:
            # include in binary array
            zarr[row, pos_ind] = 1
            prev = arr[row, pos_ind]
    elif ind == lastcol:
      # only need negative direction
      prev = brightest
      neg_ind = ind
      while neg_ind > firstcol:
        neg_ind = neg_ind - 1
        if abs((arr[row, neg_ind]-arr[row, ind]) \
          / max(arr[row,ind],arr[row,neg_ind])) < 0.1:
            # include in binary array
            zarr[row, neg_ind] = 1
            prev = arr[row, neg_ind]
    else:
      print('Bad column found in row %i, col %i' %(row, ind))
  # marks end of a row
  
  return zarr # binary array



def make_binary_thresh(imfile, thresh=20):
  # uses a hard threshold to determine which pixels to include
  if type(imfile) is str:
    img = Image.open(imfile)
    arr = np.array(img)
  # arr = invert(arr)
  else:
    arr = imfile
  try:
    L = len(arr[0,0])
  except:
    L = 1
  if L > 1:
    print('FLATTENING DID NOT WORK!')
    return False
  # helper function - binarizes the matrix
  def betch(arr):
    bets = np.zeros(np.shape(arr)[:2])
    for i in range(len(arr)):
      for j in range(len(arr[0])):
        if arr[i,j] > thresh:
          bets[i,j] = 1
    return bets
  zarr = betch(arr)
  return zarr



def make_binary_fromhist(imfile,T=1.0, sshow=False):
  img = Image.open(imfile)
  arr = np.array(img)
  # helper function
  def betch(arr):
    bets = []
    for i in range(len(arr)):
      for j in range(len(arr[0])):
        bets.append(arr[i,j])
    return bets
  bins, edges = np.histogram(betch(arr), 255) # one bin for each intensity
  edges = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
  target = bins[-1]*T # last bin should be 255, mult by the provided T val
  curr, ind = np.inf, -1
  while curr > target:
    ind = ind+1
    curr = bins[ind]
  thresh = edges[ind-1]
  # make binary image
  zarr = np.zeros(np.shape(arr))
  for i in range(len(arr)):
    for j in range(len(arr[0])):
      if arr[i,j] > thresh:
        if sshow == True:
          zarr[i,j] = 255
        else:
          zarr[i,j] = 1
  if sshow == True:
    zimg = Image.fromarray(zarr)
    zimg.show()
    
  print('file %s done' %imfile)
  return zarr



# MEAT AND POTATOES
########################################################################

def gen_plane(pt0, pt1, voxel, M=15, sshow=False):
  """
  Assuming that the distance between skelpoints is somewhat related
  to the width of the neurite and therefore the length of the normal 
  vector will reflect this.
  """
  vec = gen_vec(pt0, pt1)
  if len(vec) != 3:
    print('Error: a vector is defined by 3 values: i,j,k')
    return None
  numpts = int(dist([0,0,0],vec)/dist([0,0,0],voxel))*M
  # i also think this is automatically scaled to pt0
  plancoords = []
  xs = np.linspace(vec[0]-voxel[0]*numpts,vec[0]+voxel[0]*numpts, numpts*2) # default is + 10 voxels
  ys = np.linspace(vec[1]-voxel[1]*numpts,vec[1]+voxel[1]*numpts, numpts*2)
  def solve(vec, x, y):
    return -((x*vec[0] + y*vec[1])/vec[2])
  for m in xs: # Xs generated first
    for n in ys: # Ys generated second
      plancoords.append([m+voxel[0],n+voxel[1],solve(vec,m,n)+voxel[2]])
      # this also scales it back to the original location
  means = [np.mean([p[0] for p in plancoords]),
           np.mean([p[1] for p in plancoords]),
           np.mean([p[2] for p in plancoords])]
  # subtract the means, add the first pt
  for p in range(len(plancoords)):
    plancoords[p][0] = plancoords[p][0] - means[0] + pt0[0]
    plancoords[p][1] = plancoords[p][1] - means[1] + pt0[1]
    plancoords[p][2] = plancoords[p][2] - means[2] + pt0[2]
  lin = [[v*i for v in vec] for i in range(-10,10)] 
  if sshow == True:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for p in plancoords:
      ax.scatter(p[0],p[1],p[2], c='b', edgecolor='b', alpha=0.1)
    for p in lin:
      ax.scatter(p[0],p[1],p[2], c='r', edgecolor='r')
    plt.show()
  return plancoords, numpts # ALSO returns numpoints


############################################## The good stuff #########

def gen_vec(pt0, pt1):
  if len(pt0) == len(pt1):
    return [pt1[i]-pt0[i] for i in range(len(pt0))]
  else:
    print('Dimension mismatch: %i and %i' %(len(pt0),len(pt1)))



def gen_plane_index(sk0, sk1, M=15):
  """
  Same as above but uses indices instead of coordinates.
  """
  vec = gen_vec(sk0, sk1)
  numpts = M
  planinds = []
  xs = [int(i) for i in np.linspace(vec[0]-numpts, vec[0]+numpts, numpts*2)]
  ys = [int(i) for i in np.linspace(vec[1]-numpts, vec[1]+numpts, numpts*2)]
  def solve(vec, x, y): # the plane equation
    return -((x*vec[0] + y*vec[1])/vec[2])
  for m in xs:
    for n in ys:
      planinds.append([m, n, int(solve(vec, m, n))])
  means = [np.mean([p[0] for p in planinds]),
           np.mean([p[1] for p in planinds]),
           np.mean([p[2] for p in planinds])]
  for p in range(len(planinds)):
    planinds[p] = [planinds[p][i]+sk0[i] for i in range(3)]
  return planinds, numpts



def get_cross_sec(sk0, sk1, stack):
  """
  Returns a cross section of the intersection; tricky because the stack
  is [z][x][y] but *everything* else is [x][y][z].
  """
  # get plane
  planinds, numpts = gen_plane_index(sk0, sk1, M=15)
  dims = np.shape(stack)
  def pt_ok(pt, dims):
    # returns false if pt is outside shape of stack
    if p[0] < 0 or p[0] > dims[1]-1:
      return False
    if p[1] < 0 or p[1] > dims[2]-1:
      return False
    if p[2] < 0 or p[2] > dims[0]-1:
      return False
    return True
  #
  cs = np.zeros([2*numpts, 2*numpts])
  x_ind, y_ind, startpt = 0, 0, None
  for p in planinds:
    if y_ind >= 2*numpts:
      y_ind = 0
      x_ind = x_ind + 1
    if pt_ok(p, dims): # if point inside stack
      # print(x_ind, y_ind)
      cs[x_ind][y_ind] = stack[p[2]][p[0]][p[1]]
      if p == sk0:
        startpt = [x_ind, y_ind]
    # this should populate the cs with binary 1 for intersections and 0 otherwise
    y_ind = y_ind + 1
  # should now have the whole plane (cross section)
  if startpt is not None:
    print('Got the start point')
  return cs, startpt



def spiral_cross_sec(sk0, sk1, stack):
  """
  Apply spiral algo to the section.
  """
  cs, startpt = get_cross_sec(sk0, sk1, stack)
  s = Spiral(cs, startpt)
  area = s.area
  rad = np.sqrt(area/np.pi)
  return rad



def get_skeleton_rads(skel, stack):
  """
  Get the radius of all cross-sections in the stack.
  """
  rads = [spiral_cross_sec(skel[s],skel[s+1], stack) for s in range(len(skel)-1)]
  if len(rads) != len(skel)-1:
    print('Did not get all the sections')
  else:
    print('got them all!')
  return rads




############################################## end of good stuffs ########


def scale_plane(plancoords, pt0, voxel):
  vcoords = [] # ints for indexing the 3-D array self.varr
  for p in plancoords:
    p = [p[i]+pt0[i] for i in range(3)]
    vcoords.append([int(p[i]/voxel[i]) for i in range(3)])
  return vcoords


def plane_XY(plancoords, numpts):
  xdiff = dist(plancoords[0], plancoords[1])
  ydiff = dist(plancoords[0], plancoords[numpts])
  return [xdiff, ydiff]


def return_cross_sec_array(plane, varr, numpts, switch=False):
  """
  Plane is [i,j,k] triplets of the normal-plane, not [x,y,z] coordinates.
  Varr is the gigantic 3-D array of the voxelized image. This function
  returns a plane-sized 2-D array to be spiraled.
  """
  # plane = [plane[len(plane)-i-1] for i in range(len(plane))]
  sarr = np.zeros((numpts*2, numpts*2)) # get the whole plane, keep skelpoint at center
  means = []
  for i in range(3):
    means.append(int(np.mean([p[i] for p in plane])))
  log = []
  #sarr[means[0]][means[
  # SHOULD switch varr[z,x,y] to match plane[x,y,z], not sure why this works...
  for m in range(len(sarr)): # 'i' values first
    for n in range(len(sarr[m])): # 'j' values of new array
      if switch:
        vals = [int(i) for i in [plane[m*2*numpts+n][1],plane[m*2*numpts+n][0], plane[m*2*numpts+n][1]]]
      else:                                            # vals = [z,x,y]
        vals = [int(i) for i in plane[m*2*numpts + n]] # vals = [x,y,z]
      
      if vals[0] < 0 or vals[0] > np.shape(varr)[0]-1: # varr = [z,x,y]
        sarr[m][n] = 0 # X
      elif vals[1] < 0 or vals[1] > np.shape(varr)[1]-1:
        sarr[m][n] = 0 # Y
      elif vals[2] < 0 or vals[2] > np.shape(varr)[2]-1:
        sarr[m][n] = 0 # Z
      else: ####
        if switch:
          if varr[vals[2]][vals[0]][vals[1]] > 0:
            sarr[m][n] = 1
        elif not switch: ###
          if varr[vals[0]][vals[1]][vals[2]] > 0:
            sarr[m][n] = 1
            log.append([dist(vals, means), m,n])

    # done with cols
  # done with rows
  newlog = sorted(log)
  print(newlog)
  if len(newlog) > 0:
    return sarr, [newlog[0][1],newlog[0][2]]
  else:
    return sarr, None
  """
  sarr = np.zeros([numpts, numpts])
  for m in range(numpts):
    for n in range(numpts):
      # if the array == 1 at the plane voxel value, make a 1, else 0
      shape = np.shape(varr)
      # if the plane is past the reaches of the varr, make it zero
      if plane[(m*numpts)+n][0] > shape[0]-1 or plane[(m*numpts)+n][1] > shape[1]-1 or plane[(m*numpts)+n][2] > shape[2]-1:
        sarr[m][n] = 0
      else:
        try:
          if varr[plane[(m*numpts)+n][0]][plane[(m*numpts)+n][1]][plane[(m*numpts)+n][2]] == 1:
            sarr[m][n] = 1
          else:
            sarr[m][n] = 0
        except:
          sarr[m][n] = 0
  # now 1s should be where the plane normal to the skelpoints vector
  # intersects with supra-threshold image voxel values and 0 otherwise
  return sarr
  """
  


# this mapping function will parallelize image loading and thresholding
def import_images(folder, par=True, ttime=True):
  """
  This function loads images from a folder as PIL Image files and
  thresholds them, creating a list of z-slices to be turned into a matrix
  This version is not currently used.
  """
  fils = [os.listdir(folder)]
  def keep_tifs(rawlist):
    tiflist = []
    for f in rawlist:
      if len(f.split('.'))>1:
        if f.split('.')[1] == 'tif':
          tiflist.append(f)
    return tiflist
  tiflist = keep_tifs(fils)
  newtiflist = [folder+f for f in tiflist].sort() # alphabetize
  tifobjs = [load_img_array(f) for f in tiflist]
  
  # here start parallel stuff
  if par or ttime:
    start_time_par = timer()
    pool = Pool(8)
    results_par = pool.map(show_at_thresh, tifobjs)
    pool.close()
    pool.join()
    total_time_par = timer() - start_time_par
  # or non-parallel stuff
  elif par==False or ttime:
    start_time_nopar = timer()
    results_nopar = [show_at_thresh(f) for f in newtiflist]
    total_time_nopar = timer() - start_time_nopar
  print('Time for parallel: %.2f seconds' % total_time_par)
  print('Time for non-parallel: %.2f seconds' % total_time_nopar)
  
  return results_par, results_nopar
  
  


def save_coords(coords, fname='tempcoords.p'):
  pickle.dump(coords, open(fname, 'wb'))
  print('Coordinates written to %s as pickle' %fname)
  return



def get_voxel_locations(folder, fname, voxel=[0.176,0.176,0.38], ssave=False):
  # uses raw threshold function
  # get images as list
  fils = os.listdir(folder)
  def keep_tifs(rawlist):
    tiflist = []
    for f in rawlist:
      if len(f.split('.'))>1:
        if f.split('.')[1] == 'tif':
          tiflist.append(f)
    return tiflist
  tiflist = keep_tifs(fils) # this indexing may be removed if needed
  newtiflist = [folder+f for f in tiflist]
  newtiflist.sort() # alphabetize
  # commandeer all cores
  stime = timer()
  # pool = Pool()
  darr = [make_binary_thresh(i) for i in newtiflist] # adjust threshold in function
  # pool.close()
  # pool.join()
  print('Time taken for retrieving coordinates: %.2f' %(timer()-stime))
  # send to matrix2coords to get tuples back
  coords = matrix2coords(darr, voxel)
  if ssave:
    # save this
    save_coords(coords, fname)
  
  return coords, darr




### testing shit
"""
start = timer()
zarr = [make255(i) for i in result]
print('time taken: %.2f' %(timer()-start))

start = timer()
pool = Pool()
zzarr = pool.map(make255, result)
print('time taken: %.2f' %(timer()-start))
"""


############################ CONTROL ##############################

# usage: python imageMatrix.py[0] tiffDirectory/[1] outFile[2]
if __name__ == '__main__':
  arguments = sys.argv
  directory = arguments[1]
  if len(arguments) > 2:
    outfile = arguments[2]
  else:
    outfile = 'temp_coords.p'
  get_voxel_locations(directory, outfile, [1,1,1], False)
 







