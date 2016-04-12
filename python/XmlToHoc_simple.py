# XmlToHoc.py - converts an xml file created in knossos to a hoc file
# usage: python XmlToHoc_simple.py skelfile.nml hocfile (opt)


import os, sys
import numpy as np
import networkx as nx
import copy


class SkelHoc():
  def __init__(self, xmlFile, hocFile='hoc_from_xml.hoc'):
    self.xmlfile = xmlFile
    self.hocfile = hocFile
    self.nodes = {'*': None}
    self.nodelist = []
    self.edges = []
    self.branchpoints = []
    self.branchnodes = []
    self.branchset = []
    self.segments = []
    self.connections = []
    self.connlist = []
    self.sources, self.targets = None, None
    self.dist = None
    self.nochange = 0
    self.getrad = False
    self.used_nodes = []
    self.connver = 'dic'
    
    # now do shit
    self.read_xml()
    self.q_distance()
    self.make_segments()
    #self.refine_connections_2()
    self.refine_connections_3()
    #self.refine_connections()
    self.write_hoc()
  
  
  
  ####### Read File #######
  def add_node(self, line):
    splits = line.split('"')
    ID, x, y, z = splits[1], float(splits[5]), float(splits[7]), \
                       float(splits[9])
    thing = [x,y,z]
    if self.getrad:
      thing.append(float(splits[3]))
    self.nodes[ID] = thing
    self.nodelist.append(int(ID))
    return self
  
  def add_edge(self, line):
    splits = line.split('"')
    e = {'source': int(splits[1]), 'target': int(splits[3])}
    self.edges.append(e)
    return self
  
  def add_branchpoint(self, line):
    splits = line.split('"')
    self.branchpoints.append(splits[1])
    return self
  
  def read_xml(self):
    with open(self.xmlfile, 'r') as fIn:
      lineNum = 0
      for line in fIn:
        lineNum = lineNum + 1
        splitLine = line.split(None)[0].split('<')[1]
        if splitLine == 'node':
          self.add_node(line)
        elif splitLine == 'edge':
          self.add_edge(line)
        elif splitLine == 'branchpoint':
          self.add_branchpoint(line)
        
    return self
  
  
  def q_distance(self):
    def dist(pt0, pt1):
      if len(pt0) != len(pt1):
        print('Dimension mismatch between pt0 and pt1')
        return
      return np.sqrt(sum([(pt1[i]-pt0[i])**2 for i in range(len(pt0))]))
    mindist = np.inf
    nlist = list(self.nodes.keys())
    nlist.pop(nlist.index('*'))
    for n in range(1,len(nlist)):
      curr_dist = dist(self.nodes[nlist[0]][:3], self.nodes[nlist[n]][:3])
      if curr_dist < mindist and curr_dist != mindist:
        mindist = curr_dist
    self.dist = mindist
    print('Quantal distance is %.5f ' %mindist)
    return self
  
  
  
  ####### Do analysis ########
  
  def make_segments(self):
    # Source is 0th end, target is 1th end
    self.sources = [i['source'] for i in self.edges]
    self.targets = [i['target'] for i in self.edges]
    # connections always 1th to 0th
    for s in range(len(self.sources)):
      if self.targets.count(self.sources[s]) > 0:
        inds = [i for i in range(len(self.targets)) if self.targets[i] == self.sources[s]]
        for ix in inds:
          con = [s, ix]
          if con not in self.connections and [con[1],con[0]] not in self.connections:
            self.connections.append(con)
    for t in range(len(self.targets)):
      if self.sources.count(self.targets[t]) > 0:
        inds = [i for i in range(len(self.sources)) if self.sources[i] == self.targets[t]]
        for ix in inds:
          con = [t, ix]
          if con not in self.connections and [con[1],con[0]] not in self.connections:
            self.connections.append(con)
    
    # eliminate redundancy
    newcons = []
    for con in self.connections:
      if con not in newcons or con.reverse() not in newcons:
        newcons.append(con)
    self.connections = newcons
    
    return self
  
  
  
  def add_connection(self, one, zero):
    self.connlist.append({'1': one, '0': zero})
    return self
  
  
  
  def refine_connections(self):
    # supposedly this finds which point 2 segments have in common
    # and creates a connection based on that, but it clearly doesn't work
    self.segments = list(zip(self.sources, self.targets))
    for con in self.connections:
      # target=0, source = 1
      if self.targets[con[1]] == self.sources[con[0]]:
        self.add_connection(con[1],con[0])
      elif self.targets[con[0]] == self.sources[con[1]]:
        self.add_connection(con[0],con[1])
      else:
        print('could not find connection points for [%i, %i]' \
              %(con[0],con[1]))
    
    for con in self.connlist:
      pass
    
    return self
  
  
  
  def refine_connections_2(self):
    # this one uses only overlapping nodes to connect segments
    # only the copies are popped, all saved stuff should come from self.segments
    self.segments = list(zip(self.sources, self.targets))
    self.nodes.pop('*')
    tempnodes = self.nodes
    tempsegs = copy.deepcopy(self.segments)
    conn_positions, conn_indices = [], []
    # for each node
    for n in tempnodes.keys():
      currnode = tempnodes[n]
      # print('currnode is %s' %n)
      self.used_nodes.append(n)
      # tempnodes.pop(n) ## can't change dict size during iteration
      # find which first (n) segment (sn=index)
      for sn in range(len(self.segments)):
        if int(n) in self.segments[sn]:
          # print('found matching segment')
          n_pos = self.segments[sn].index(int(n))
          # check if that node is somewhere else
          for m in tempnodes.keys():
            if tempnodes[m] == currnode and m not in used_nodes:
              print('Found matching node')
              # find the matching (m) segment(s) (sm=index)
              for sm in range(len(self.segments)):
                if int(m) in self.segments[sm]:
                  # found a possible match, make sure not already there in some order
                  if [sn, sm] not in conn_indices and [sm, sn] not in conn_indices and sm!=sn:
                    m_pos = self.segments[sm].index(int(m))
                    conn_indices.append([sn, sm])
                    conn_positions.append([n_pos, m_pos])
              # if node m was a match, pop it
              # tempnodes.pop(m) ## can't change dict size during iteration
              self.used_nodes.append(m)
    # now make connlist
    print(conn_indices)
    for k in range(len(conn_indices)):
      self.connlist.append([conn_indices[k][i] for i in conn_positions[k]])
  
    return self
  
  
  
  def refine_connections_3(self):
    # this does same as 2 but with segment-node index instead of nodes
    print('Called Refine #3')
    self.segments = list(zip(self.sources, self.targets))
    allnodes, connections = [], []
    # for each segment s
    for s in range(len(self.segments)):
      # for each node n in s
      for n in range(len(self.segments[s])):
        first = self.segments[s][n]
        # add it to used_nodes first
        # used_segs.append(self.segments[s][n])
        # now go through the possible segments
        for ss in range(len(self.segments)):
          for m in range(len(self.segments[ss])):
            second = self.segments[ss][m]
            # print('comparing %i to %i' %(self.segments[s][n], self.segments[ss][m]))
            if first == second:
              keylist = [list(k.keys()) for k in connections]
              if [str(ss),str(s)] not in keylist and [str(s),str(ss)] not in keylist:
                #print('Found match')
                connections.append({str(ss): m, 
                                    str(s): n})
                print('A redundant connection has been skipped')
              # used_segs.append(self.segments[ss][m])
    # show connections
    newconns = []
    for k in range(len(connections)):
      if len(connections[k]) == 2:
        newconns.append(connections[k])
    
    print(newconns)
    self.connlist = newconns
    
    return self
  
  
  
  
  
  ####### write the hoc file #####
  def write_hoc(self):
    with open(self.hocfile, 'w') as fOut:
      
      def start_hoc(total_segs):
        fOut.write('create filament_999[%i]\n' %total_segs)
        return
      
      def create_filament(sname):
        fOut.write('filament_999[%i] {\n' %sname)
        fOut.write('  pt3dclear()\n')
        return
      
      def pt3dadd(node):
        if self.getrad:
          rad = node[3]
        else:
          rad = 1.5
        fOut.write('  pt3dadd(%f, %f, %f, %f)\n'
                   %(node[0], node[1], node[2], rad))
        return
      
      def end_filament():
        fOut.write('}\n')
        return
      
      def connect_filaments(conns):
        if self.connver=='dic': # this means refine_2 was not used
          # conns is a dict of shape {'0': filA, '1': filB}
          keys, vals = list(conns.keys()), list(conns.values())
          fOut.write('connect filament_999[%i](%i), filament_999[%i](%i)\n'
                     %(int(keys[0]),vals[0],int(keys[1]),vals[1]))
        else:
          # conns is a list of shape [filA, filB] for 0th and 1th ends, respectively
          fOut.write('connect filament_999[%i](0.0), filament_999[%i](1.0)\n'
                     %(conns[0], conns[1]))
        return
      
      # start the hoc file
      start_hoc(len(self.segments))
      # create all the individual filaments
      seg_ind = 0
      for seg in self.segments:
        create_filament(seg_ind)
        for s in seg:
          pt3dadd(self.nodes[str(s)])
        end_filament()
        seg_ind = seg_ind + 1
      
      # create the connection matrix
      for con in self.connlist:
        connect_filaments(con)
      
    return
  




###### helper functions ########

def split_sources_targets(ldict):
  # splits a list of dicts{'1': #, '0': #} into sources and targets
  sources, targets = [], []
  for i in ldict:
    sources.append(i['1'])
    targets.append(i['0'])
  return sources, targets, list(zip(sources, targets))
  





##########################
if __name__ == '__main__':
  
  args = sys.argv
  if len(args) < 2:
    print('Need a skeleton file to analyze')
  skelfile = args[1]
  if len(args) > 2:
    hocfile = args[2]
    h = SkelHoc(skelfile, hocfile)
  else:
    h = SkelHoc(skelfile)




#########################################################
# Garbage from previous versions
"""
def create_segments(self):

#Run through edges, nodes and branchpoints to return the actual
#segments. 

# start with a random edge (or node)
self.sources = [i['source'] for i in self.edges]
self.targets = [i['target'] for i in self.edges]
for i in range(len(self.sources)):
c = self.sources.count(self.sources[i]) + \
self.targets.count(self.sources[i])
if c >= 3 or c == 1 and self.sources[i] not in self.branchnodes:
self.branchnodes.append(self.sources[i]) # add branchpoints and tips to branchnodes
for coo in range(c):
self.branchset.append(self.sources[i])
c = self.sources.count(self.targets[i]) + \
self.targets.count(self.targets[i])
if c >= 3 or c == 1 and self.targets[i] not in self.branchnodes:
self.branchnodes.append(self.targets[i])
for coo in range(c):
self.branchset.append(self.targets[i])
self.branchnodes = set(self.branchnodes)

print('Found %i branchnodes or tips' %len(self.branchnodes))
print('Branchset is %i points long' %len(self.branchset))

# print(self.branchset)
orig = int(len(self.branchset))
prev = int(len(self.branchset))
while len(self.branchset) > 0 and self.nochange < 1000:
#if len(self.branchset)%10 == 0:
#  print('%f percent branchpoints done.' %float(100-100*len(self.branchset)/orig))
print('%i branchpoints remaining ...' %len(self.branchset))
rand_ind = int(np.random.random(1)*len(self.branchset))
self.crawl(int(self.branchset[rand_ind]))
if int(len(self.branchset)) == prev:
self.nochange = self.nochange + 1
else:
prev = int(len(self.branchset))
return self



def crawl(self, b):
"""
#Current is always already added to node list.
"""
print('Crawling')
try:
nodes = [self.branchset.pop(b)] # only in scope of current function; else self.nodes
except:
return
go = False
try:
if self.branchset[b] in self.sources or self.branchset[b] in self.targets:
go = True
else:
print('%i not a source or a target!?' %self.branchset[b])
# already been popped
return
except:
return
current = b
count = 0
while go and count < 1000:
count = count + 1
if nodes[-1] in self.sources: # found match in targets
ind = self.sources.index(nodes[-1]) # find index
nodes.append(self.targets.pop(ind)) # pop it onto current nodes
self.sources.pop(ind) # remove the linker
elif nodes[-1] in self.targets:
ind = self.targets.index(nodes[-1])
nodes.append(self.sources.pop(ind))
self.targets.pop(ind)
if nodes[-1] in self.branchset:
go = False
try:
self.branchset.pop(self.branchset.index(nodes[-1]))
except:
print('Tried to pop %i from branchset but was already popped'
% current)
# end of while loop
#if nodes[0] not in self.branchpoints or nodes[-1] not in self.branchpoints:
#  print('Segment does not start and end with branchpoints:')
self.format_segment(nodes) # format the segment to prepare for hoc
# print('Finished segment for branchnode %i' %self.branchset[b])
print('done crawling')
return self



def format_segment(self, nodes):
newseg = len(self.segments.keys())
self.segments[str(newseg)] = {}
for n in range(len(nodes)):
self.segments[str(newseg)][str(n)] = \
self.nodes[str(nodes[n])]
return self



def create_connections(self):
end_segs, end0s, end1s, no_connections = [], [], [], []
if '*' in self.segments.keys():
self.segments.pop('*')
for c in self.segments.keys():
end_segs.append(c)
end0s.append(self.segments[c]['0'])
last = len(self.segments[c])-1
end1s.append(self.segments[c][str(last)])
print(len(end0s), len(end1s))
for e in range(len(end_segs)):
if end1s[e] in end0s:
# found a segment that connects to e's 1th end
f = end0s.index(end1s[e])
newcon = {'0': end_segs[f], '1': end_segs[e]}
if newcon not in self.connections:
self.connections.append(newcon)
if end0s[e] in end1s:
# found a segment that connects to e's 0th end
f = end1s.index(end0s[e])
newcon = {'0': end_segs[e], '1': end_segs[f]}
if newcon not in self.connections:
self.connections.append(newcon)
if end0s[e] not in end1s and end1s[e] not in end0s:
no_connections.append(e)
print('Total connections: %i. No connections found for %i segments.'
%(len(self.connections), len(no_connections)))
# end of for loop
return self

"""


























