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
from neuron_readExportedGeometry import *


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



def switch_binary(arr):
  zarr = np.zeros(np.shape(arr))
  for i in range(len(arr)):
    for j in range(len(arr[0])):
      if arr[i,j] == 1:
        zarr[i,j] = 125
      elif arr[i,j] == 2:
        zarr[i,j] = 255
  return zarr



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


def show_slice(stack, z=None):
  # Show a particular slice of interest
  if z is None:
    z = int(len(stack)/2)
  for i in range(len(stack)):
    for j in range(len(stack[i])):
      for k in range(len(stack[i][j])):
        if stack[i][j][k] > 0:
          stack[i][j][k] = 255
  img = Image.fromarray(stack[z])
  img.show()
  return



# BINARY OPTIONS 
########################################################################

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



#######################################################################
############################################## The good stuff #########
#######################################################################


def gen_vec(pt0, pt1):
  if len(pt0) == len(pt1):
    return [pt1[i]-pt0[i] for i in range(len(pt0))]
  else:
    print('Dimension mismatch: %i and %i' %(len(pt0),len(pt1)))



def gen_plane_index(sk0, sk1, M=50, as_int=False):
  """
  Same as above but uses indices instead of coordinates. X,Y are ints
  because they refer to indices; Z can be float for interpolation.
  """
  vec = gen_vec(sk0, sk1)
  numpts = M
  planinds = []
  xs = [int(i) for i in np.linspace(vec[0]-numpts, vec[0]+numpts, numpts*2)]
  ys = [int(i) for i in np.linspace(vec[1]-numpts, vec[1]+numpts, numpts*2)]
  def solve(vec, x, y): # the plane equation
    if vec[2]==0:
      return -((x*vec[0] + y*vec[1])/0.01)
    return -((x*vec[0] + y*vec[1])/vec[2]) # 
  for m in xs:
    for n in ys:
      if as_int is True:
        planinds.append([m, n, int(solve(vec, m, n))])
      elif as_int is False:
        planinds.append([m, n, float(solve(vec, m, n))])
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
  planinds, numpts = gen_plane_index(sk0, sk1, M=1000)
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
      cs[int(x_ind)][int(y_ind)] = stack[int(p[2])][int(p[0])][int(p[1])] # Rearrange for 3D
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
  rad = np.sqrt(area/np.pi) # Assumes cylindrical cross-sec
  return rad



def get_geo_rads(geo, stack, switchxy=True):
  """
  Get the radius of all skeleton points from a geo file.
  """
  rads = np.zeros(len(geo.nodes))
  for b in geo.branches:
    for n in range(len(b.nodes)-1):
      # If this radius has already been retrieved, forget it
      if rads[geo.nodes.index(b.nodes[n])] != 0:
        r = rads[geo.nodes.index(b.nodes[n])]
      else:
        if switchxy:
          sk0 = [int(i) for i in [b.nodes[n].y, b.nodes[n].x, b.nodes[n].z]]
          sk1 = [int(i) for i in [b.nodes[n+1].y, b.nodes[n+1].x, b.nodes[n+1].z]]
        else:
          sk0 = [int(i) for i in [b.nodes[n].x, b.nodes[n].y, b.nodes[n].z]]
          sk1 = [int(i) for i in [b.nodes[n+1].x, b.nodes[n+1].y, b.nodes[n+1].z]]
        print(sk0, sk1)
        r = spiral_cross_sec(sk0, sk1, stack)
        print('Got radius: %.2f' %r)
        rads[geo.nodes.index(b.nodes[n])] = r
      # If it's the last node, give it the same radius as the prev one
      if (n+1) == len(b.nodes)-1:
        rads[geo.nodes.index(b.nodes[n+1])] = r
  return rads



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



def add_radius(hocFile, geo, rads, newHocFile='new_hoc_radius.hoc'):
  """
  The hocFile should have already been loaded, but this is easier. This
  function will match the loaded geo objects (especially nodes) to the
  node list loaded previously. It uses the geo object to map the new radii
  onto the new hocfile.
  """
  geo2 = demoReadsilent(hocFile)
  # Nodelist should be of same length and index as rads...
  nodelist = [[n.x, n.y, n.z] for n in geo.nodes]
  if len(nodelist) != len(geo2.nodes):
    print('Some nodes are repeated!') # Not a big deal
  rcount, norad, yesrad = -1, 0, 0
  # Write the new hocfile; only change the pt3dadd rad values
  with open(newHocFile, 'w') as fOut:
    with open(hocFile, 'r') as fIn:
      for line in fIn:
        if line:
          splitLine = line.split('(')
          # If it's a line with a radius, write the new line
          if splitLine[0] == 'pt3dadd' or splitLine[0] == '  pt3dadd':
            rcount = rcount + 1
            x, y, z, _ = splitLine[1].split(',')
            x, y, z = float(x), float(y), float(z)
            try:
              rad = rads[nodelist.index([x, y, z])]
              yesrad = yesrad + 1
            except:
              norad = norad + 1
              rad = np.mean(rads)
            try:
              fOut.write('pt3dadd(%f, %f, %f, %f)' 
                         %(float(x), float(y), float(z), rad))
            except:
              print(x, y, z, rad)
          # Else, just write the original line
          else:
            fOut.write(line)
  print('Missed %i nodes, got %i nodes' %(norad, yesrad))
  return
  


def get_branch_skelpts(geo, voxel=[1,1,1]):
  """
  Return skelpts organized by branches (list of lists). Voxel is the dims
  of a voxel in the image (default is for indices).
  """
  if type(geo) is str:
    try:
      geo = demoReadsilent(geo)
    except:
      print('A string should point to the hocfile!')
  skelpts = []
  for b in geo.branches:
    t = []
    for n in b.nodes:
      t.append([int(n.x/voxel[0]), int(n.y/voxel[1]), int(n.z/voxel[2])])
    skelpts.append(t)
  return skelpts


########################## FAST ACCESS #################################

def downsample(skelpts, by=3):
  # Downsample the skeleton nodes by an integer, skips 2 nodes for every
  # 2 nodes taken.
  newskel  = [ [] for i in skelpts]
  for s in range(len(skelpts)):
    if len(skelpts[s]) <= by and len(skelpts[s]) >= 2:
      newskel[s] = [skelpts[s][0], skelpts[s][1]]
    elif len(skelpts[s]) < 2:
      pass
    elif len(skelpts[s]) > by and len(skelpts[s]) <= by*2:
      newskel[s] = [skelpts[s][0], skelpts[s][1]]
    else: # More than 2*'by' nodes in this seg
      for u in range(int(len(skelpts[s])/by)):
        newskel[s].append(skelpts[s][u*by])
  return newskel



def cross_sec_locs(sk0, sk1, dims):
  """
  Returns the locations/indices of a cross-section.
  The output is an array (500+500) x (500+500) and each element of that
  array is a 3-tuple that tells the location in the stack to get the 
  value. sk0 = [x,y,z], dims = [z,x,y]
  """
  # get plane
  planinds, numpts = gen_plane_index(sk0, sk1, M=500)
  def pt_ok(pt, dims):
    # returns false if pt is outside shape of stack
    if pt[0] < 0 or pt[0] > dims[1]-1:
      return False
    if pt[1] < 0 or pt[1] > dims[2]-1:
      return False
    if pt[2] < 0 or pt[2] > dims[0]-1:
      return False
    return True
  #
  cs = [ [ [] for i in range(numpts*2)] for j in range(numpts*2) ]
  x_ind, y_ind, startpt = 0, 0, None
  for p in planinds:
    if y_ind >= 2*numpts: # Reset the counter
      y_ind = 0
      x_ind = x_ind + 1
    if pt_ok(p, dims): # if point inside stack
      # print(x_ind, y_ind)
      cs[int(x_ind)][int(y_ind)] = p
      if p == sk0:
        startpt = [x_ind, y_ind]
    # this should populate the cs with 3-tuple location of data values
    y_ind = y_ind + 1
  # should now have the whole plane (cross section)
  if startpt is not None:
    print('Got the start point')
  return cs, startpt




############################################## INTEGER Zs

def fill_cross_sec(cross_secs, root_dir): 
  """
  Pass all the cross-sections and the directory where the image stack
  is contained. cross_secs is a list of 3-tuple arrays from 
  cross_secs_locs.
  """
  def place_values(cs, zslice, slicenum):
    # Populate as much of cross-sec as possible from the given zslice
    for i in range(len(cs)):
      for j in range(len(cs[i])):
        try:
          if cs[i][j][2] == slicenum:
            cs[i][j] = zslice[ cs[i][j][0], cs[i][j][1] ] # ...[2] already used
        except:
          pass
    return cs
  #
  print('Getting which slices are needed for each cross-section ...')
  slices = []
  for cs in cross_secs:
    if cross_secs.index(cs) % 100 == 0:
      print('%i / %i cross-secs examined so far.' %(cross_secs.index(cs), len(cross_secs)))
      slices.append(which_slices(cs))
  # Got all the slice info
  print('Filling in the cross-sections ...')
  fils = os.listdir(root_dir)
  fils = [f for f in fils if f.split('.')[-1] == 'tif']
  fils.sort()
  fils = [root_dir + f for f in fils]
  # For each z slice, find all the cross-secs that need data from it
  for u in range(len(fils)):
    arr = np.asarray(Image.open(fils[u]))
    for which in range(len(slices)):
      if u in slices[which]:
        cross_secs[which] = place_values(cross_secs[which], arr, u)
    if u % 10 == 0:
      print('%i / %i slices processed so far ... ' %(u, len(fils)))
  # Every cross-sec should now be populated
  return cross_secs



def which_slices(cs):
  """
  Determines which z slices are needed to reconstruct the cross-sec
  Pass one cross_section at a time
  """
  wh = []
  for i in cs:
    for j in i:
      try:
        if len(j) == 3: # If the point is within bounds
          wh.append(int(j[2])) # Add the z slice number
      except:
        pass
  return list(set(wh))



############################################## FLOAT Zs

def interp_voxels(cs, slice1num, slice1, slice2num, slice2, voxel):
  """
  Interpolate between two slices to find the correct value.
  Fill as much of the slice as possible.
  """
  num_interp = int(max([voxel[2]/voxel[0], voxel[2]/voxel[1]])) # Num of interp pts
  for i in range(len(cs)):
    for j in range(len(cs[i])):
      if type(cs[i][j]) is list:
        if len(cs[i][j]) == 3:
          if cs[i][j][2] == slice1num: # If this Z is contained in current slice
              column = np.linspace(slice1[cs[i][j][0]][cs[i][j][1]],
                                   slice2[cs[i][j][0]][cs[i][j][1]],
                                   num_interp)
              cs[i][j] = column[int( (cs[i][j][2]%1)*num_interp )]
  return cs
  


def float_cross_secs(cross_secs, root_dir, voxel):
  """
  This populates the cross-secs where x & y values are ints but z is float.
  This is done by interpolating between the z slices.
  """
  w_slices = [] # Find which slices are necessary
  for cs in cross_secs:
    if cross_secs.index(cs) % 100 == 0:
      print('%i / %i cross-secs examined so far.' %(cross_secs.index(cs), len(cross_secs)))
      w_slices.append(which_slices(cs))
  print('Filling and interpolating for cross-sections...')
  fils = os.listdir(root_dir) # Get files from directory
  fils = [f for f in fils if f.split('.')[-1].lower() == 'tif']
  fils.sort()
  fils = [root_dir + f for f in fils]
  #
  # Iterate through the slices to get the Z-values
  curr_slice = np.asarray(Image.open(fils[0]))
  for s in range(1,len(fils)): # For each slice ...
    next_slice = np.asarray(Image.open(fils[s])) # ... read the next slice ...
    for which in range(len(w_slices)): 
      if s-1 in w_slices[which]: # ... if this cross-sec requires the current slice ...
        # Fill in the missing pieces for that slice
        cross_secs[which] = interp_voxels(cross_secs[which], 
                                          slice1num=s-1, slice1=curr_slice,
                                          slice2num=s, slice2=next_slice,
                                          voxel=voxel)
    curr_slice = next_slice # Replace the current slice
    if s % 10 == 0:
      print('%i / %i slices examined ...' %(s, len(fils)))
  #
  # Once all slices have been filled in, replace any missing vals with 0s
  for cs in range(len(cross_secs)):
    for i in range(len(cross_secs[cs])):
      for j in range(len(cross_secs[cs][i])):
        if type(cross_secs[cs][i][j]) is list:
          cross_secs[cs][i][j] = 0
  return cross_secs

############################ end of float & int cross-sects


def integer_cross_sec(cross_sec):
  """
  If a cross-sec tuple pointed to a non-existent value (outside range),
  it was left as an empty list ([]). This replaces those with zeros.
  """
  new_cs = []
  for i in range(len(cross_sec)):
    temp = [] # For each row
    for j in cross_sec[i]:
      try:
        temp.append(int(j))
      except:
        temp.append(0)
    new_cs.append(temp)
  return new_cs

  

# Get number of triplet/tuples in array
def num_of_triplets(arr):
  cnt = 0
  for i in arr:
      for j in i:
          if len(j) == 3:
              cnt = cnt + 1
  return cnt
  





def cs_from_skeleton(skelpts, dims):
  """
  Get the cross-sec indices for all of the relevant skel points
  """
  cross_sections, center_pts = [], [] # Keep track of the center points
  for s in skelpts:
    for k in range(len(s)-1):
      pass #
  #
  return









######################### end of good stuffs #############################



##########################################################################
######################### for beta testing





######################### threshold helpers ##############################


def suggest_thresh(imfile, th_range=None):
  """
  """
  # Helper function to show thresholded images
  def show_t(arr, thresh):
    t_arr = np.zeros(arr.shape)
    for i in range(len(arr)):
      for j in range(len(arr[i])):
        if arr[i,j] > thresh:
          t_arr[i,j] = 1
    t_arr = make255(t_arr)
    return Image.fromarray(t_arr)
  arr = np.array(Image.open(imfile))
  if len(arr.shape) > 2:
    print('array has too many dimensions')
    return None
  if th_range is None:
    minval, maxval = np.inf, 0
    for i in arr:
      for j in i:
        if j < minval:
          minval = j
        if j > maxval:
          maxval = j
  elif len(th_range) == 2:
    minval, maxval = min(th_range), max(th_range)
  # linspace 9 thresholds through and find which is best
  threshes = np.linspace(minval, maxval, 9)
  fig = plt.figure()
  for t in range(len(threshes)):
    ax = fig.add_subplot(3, 3, t)
    img = show_t(arr, threshes[t])
    ax.imshow(img)
    ax.set_title('Threshold of %.2f' %threshes[t])
  plt.show()
  return



def import_images(folder, thresh=False, par=True):
  """
  This function loads images from a folder as PIL Image files and
  thresholds them with the given scalar. Par parallelizes this, ttime does
  both serial and parallel and times them. 
  *Must be 'tif' extension!*
  """
  fils = [f for f in os.listdir(folder) if f.split('.')[-1]=='tif']
  fils.sort()
  if folder[-1] != '/':
    folder = folder+'/'
  newtiflist = [folder+f for f in fils] # alphabetize
  # tifobjs = [load_img_array(f) for f in tiflist]
  
  def thresh_load(imfile):
    im_arr = np.array(Image.open(imfile))
    arr = np.zeros(im_arr.shape)
    for i in range(len(im_arr)):
      for j in range(len(im_arr[i])):
        if im_arr[i,j] >= thresh:
          arr[i,j] = 1
    return arr
  def no_thresh_load(imfile):
    im_arr = np.array(Image.open(imfile))
    return im_arr
    
  # here start parallel stuff
  if par:
    start_time_par = timer()
    pool = Pool(8)
    if thresh is not False:
      results_par = pool.map(thresh_load, newtiflist)
    else:
      results_par = pool.map(no_thresh_load, newtiflist)
    pool.close()
    pool.join()
    total_time_par = timer() - start_time_par
    print('Time for parallel: %.2f seconds' % total_time_par)
    return results_par
  # or non-parallel stuff
  else:
    start_time_nopar = timer()
    if thresh is not False:
      results_nopar = [thresh_load(f) for f in newtiflist]
    else:
      results_nopar = [no_thresh_load(f) for f in newtiflist]
    total_time_nopar = timer() - start_time_nopar
    print('Time for non-parallel: %.2f seconds' % total_time_nopar)
    return results_nopar



############################ CONTROL ##############################

# usage: python imageMatrix.py[0] tiffDirectory/[1] outFile[2]
if __name__ == '__main__':
  arguments = sys.argv
  directory = arguments[1]
  if len(arguments) > 2:
    outfile = arguments[2]
    
  print('command-line usage not finished yet!')

 







