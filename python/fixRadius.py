# correct illegal radius
# usage: python fixRadius.py hocFile newhocFile

import numpy as np
import sys



#########################################################################
# import contents of hoc file into segment dict and connections
def parse_pt(segments, line, getSeg):
  if getSeg not in segments.keys():
    segments[getSeg] = []
  csv = line.split(',')
  x = float(csv[0].split('(')[1])
  y = float(csv[1])
  z = float(csv[2])
  rad = float(csv[3].split(')')[0])
  segments[getSeg].append([x, y, z, rad])
  return segments



def parse_connection(connections, line):
  splitline = line.split(None)
  zero = int(splitline[2].split('_')[1].split('(')[0])
  one = int(splitline[1].split('_')[1].split('(')[0])
  curr_dict = {'0': zero, '1': one}
  connections.append(curr_dict)
  return connections



def load_hoc(hocFile):
  segments = {}
  connections = []
  getSeg = None
  with open(hocFile, 'r') as fIn:
    for line in fIn:
      if line:
        splitline = line.split(None)
        if splitline:
          if splitline[0] == '}':
            getSeg = None
          elif getSeg:
            if splitline[0].split('_')[0] == 'section':
              pass
            elif splitline[0].split('(')[0] == 'pt3dclear':
              pass
            elif splitline[0].split('(')[0] == 'pt3dadd':
              segments = parse_pt(segments, line, getSeg)
          if splitline[0] == 'connect':
            connections = parse_connection(connections, line)
          elif splitline[0] == 'create':
            getSeg = splitline[-1].split('_')[-1]
  print('Done reading %s' %hocFile)
  return segments, connections



########################################################################
# check and fix radii
def get_radius(pt, origSeg, segments, connections):
  # given a pt with a bad radius, find what it should be
  g_segs = []
  for con in connections:
    # print(con.values(), origSeg)
    if int(origSeg) in con.values():
      for h in con.values():
        # print(con.values())
        if h is not int(origSeg):
          g_segs.append(h)
  while int(origSeg) in g_segs:
    g_segs.pop(g_segs.index(int(origSeg)))
    
  # print(g_segs)
  # add the first and last points of each potential segment
  g_pts = []
  for g in g_segs:
    g_pts.append(segments[str(g)][0])
    g_pts.append(segments[str(g)][-1])
  rad = None
  for p in g_pts:
    if p[:3] == pt[:3]:
      print('Found correct radius; switching %.5f for %.5f'
            %(pt[-1],p[-1]))
      rad = p[-1]
  return rad



def gen_allpts(segments):
  # generate the master pt list for slow_rad
  allpts, pts_segs = [], []
  for k in segments.keys():
    for n in segments[k]:
      allpts.append(n)
      pts_segs.append(k)
  return allpts, pts_segs



def slow_rad(pt, allpts, pts_segs):
  # crawl through segpoints and find the one to match the rad
  # else, it's a tip and give it a rando radius (otherwise it must
  # match another x-y-z tuple or NeuronGeometry will freak out
  for p in allpts:
    possibles = []
    if p[:3] == pt[:3] and p != pt:
      possibles.append(p[-1])
  u_possibles = list(np.unique(possibles))
  if None in u_possibles:
    u_possibles.pop(u_possibles.index(None))
  if len(u_possibles) > 1:
    print('Multiple radii possible, choosing one')
    print(u_possibles)
  elif len(u_possibles) == 1:
    return u_possibles[0]
  elif len(u_possibles) == 0:
    print('No good radius found; using 0.09')
    return 0.09



def check_rads(segments, connections):
  # find segments with only 2 members -- these are the troublemakers
  counter = 0
  allpts, pts_segs = gen_allpts(segments)
  for k in segments.keys():
    if len(segments[k]) == 2:
      for node in range(len(segments[k])):
        counter = counter + 1
        newrad = get_radius(segments[k][node], k, segments, connections)
        if newrad is None:
          newrad = slow_rad(segments[k][node], allpts, pts_segs)
        if newrad == 0.:
          newrad = 0.09
        segments[k][node][-1] = newrad
  # now do slow check to fix the 
  print('Fixed %i radii' %counter)
  return segments, connections



########################################################################
# write new hoc file
def create_section(fOut, segments, seg):
  fOut.write('create section_%s\n' %seg)
  fOut.write('section_%s {\n' %seg)
  for p in segments[seg]:
    if type(p) is list:
      # print(p)
      fOut.write('pt3dadd(%f,%f,%f,%f)\n' %(p[0],p[1],p[2],p[3]))
  fOut.write('}\n\n')
  return



def create_connection(fOut, connections):
  # write connection matrix
  fOut.write('\n')
  for con in connections:
    fOut.write('connect section_%s(1), section_%s(0)\n'
               %(con['1'], con['0']))
  return



def write_file(outFile, segments, connections):
  #print(segments, connections)
  with open(outFile, 'w') as fOut:
    # first do segments
    for seg in segments.keys():
      create_section(fOut, segments, seg)
    create_connection(fOut, connections)
  print('Create new hoc file %s' %outFile)
  return




########################################################################
if __name__ == '__main__':
  args = sys.argv
  hocFile = args[1]
  if len(args) < 3:
    outFile = 'new_hoc_radii.hoc'
  elif len(args) >= 3:
    outFile = args[2]
  segments, connections = load_hoc(hocFile)
  segments, connections = check_rads(segments, connections)
  write_file(outFile, segments, connections)
  














