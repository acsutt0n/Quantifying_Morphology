
import numpy as np
import math



def dist3(pt0, pt1):
  if len(pt0) == len(pt1) and len(pt0) == 3:
    return math.sqrt(sum([(pt0[i]-pt1[i])**2 for i in range(3)]))
  else:
    print('dimension mismatch')
    print(pt0, pt1)



def get_angle(pt0, midpt, pt1):
  if pt0 in [midpt, pt1] or pt1 in [midpt, pt0] or midpt in [pt0,pt1]:
    print('Some points are the same!')
    print(pt0, midpt, pt1)
  PT0 = dist3(pt1, midpt)
  PT1 = dist3(pt0, midpt)
  MIDPT = dist3(pt0, pt1)
  try:
    ang = math.acos( (MIDPT**2 - PT1**2 - PT0**2) / (2*PT1*PT0) )
    ang = ang*180/math.pi
  except:
    ang = 'nan'
  return ang



def find_points(seg0, seg1):
  seg0list, seg1list = [], []
  pt0where, pt1where, midwhere = None, None, None
  switchdict = {0: -1, -1: 0}
  # make a list of the node locations
  for n in seg0.nodes:
    seg0list.append([n.x,n.y,n.z])
  for n in seg1.nodes:
    seg1list.append([n.x,n.y,n.z])
    # find the common node, then use that to find the distinct ones
  for n in seg0list:
    if n in seg1list:
      midpt = n
  if seg0list.index(midpt) != 0:
    pt0where = 0
    pt0 = seg0list[0]
  else:
    pt0where = -1
    pt0 = seg0list[-1]
  if seg1list.index(midpt) != 0:
    pt1where = 0
    pt1 = seg1list[0]
  else:
    pt1where = -1
    pt1 = seg1list[-1]
  
  f = True
  if pt0 == pt1 or pt0==midpt:
    f = False
    if pt0where == 0:
      try:
        pt0=seg0list[1]
        f = True
      except:
        pass
    elif pt0where == -1:
      try:
        pt0=seglist[-2]
        f = True
      except:
        pass
  if pt0 == pt1 or pt1==midpt:
    if pt1where == 0:
      try:
        pt1=seg1list[1]
        f = True
      except:
        pass
    elif pt1where == -1:
      try:
        pt1=seg1list[-2]
        f = True
      except:
        pass
  if f == False:
    print('Tried to find new coordinates, but failed. Skipping')
  
  if pt0 in [midpt, pt1] or pt1 in [midpt, pt0] or midpt in [pt0,pt1]:
    print(seg0list, seg1list)
  #print('pt0 at %i, pt1 at %i' %(pt0where, pt1where))
  if pt1 and pt0 and midpt:
    return pt0, midpt, pt1
  else:
    print('could not figure out segments %s and %s' %(seg0.name, seg1.name))
    print(seg0list, seg1list)
    return [False]



def branch_angles(geo):
  angles = []
  for b in geo.branches:
    for n in b.neighbors:
      pts = find_points(n, b)
      if len(pts) == 3:
        pt0, midpt, pt1 = pts[0], pts[1], pts[2]
      angles.append(get_angle(pt0, midpt, pt1))
  angles = [a for a in angles if a!='nan']
  with open('temp_angles.txt', 'w') as fOut:
    for a in angles:
      fOut.write('%.10f, \n' %a)
  return angles
