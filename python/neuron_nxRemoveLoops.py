# neuron_nxRemoveLoops.py -- removes loops from hoc file based on 
#                            a networkx object and methods
# Much of this is in the ipython notebook networkx_loops.ipynb
# usage: python neuron_nxRemoveLoops.py hocIn hocOut (optional)


import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys



def as_nx_object(infile):
  """
  Get the 'nodes' and edges from a hoc file; return with nx.Graph
  """
  edges, nodes = [], []
  with open(infile, 'r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        if type(splitLine) is list and len(splitLine) > 0:
          if splitLine[0] == 'connect':
            e1 = splitLine[1].split('[')[1].split(']')[0]
            if e1 not in nodes:
              nodes.append(e1)
            e1 = e1 + '.' + splitLine[1].split('(')[1].split(')')[0]
            e2 = splitLine[2].split('[')[1].split(']')[0]
            if e2 not in nodes:
              nodes.append(e2)
            e2 = e2 + '.' + splitLine[2].split('(')[1].split(')')[0]
            edges.append([e1, e2])
  # Add the "intrinsic" edges that connect the 0th end to 1th 
  # end of each node
  for s in nodes:
    edges.append([s+'.'+'0', s+'.'+'1'])
  G = nx.Graph()
  G.add_edges_from(edges)
  return G, edges, nodes



def remove_loops(infile, G=None, edges=None, nodes=None):
  """
  Use the nx object to remove the loops from the Graph object. Returns
  nx.Graph object. Together with the segment node-list can create a 
  hoc file without any loops.
  """
  # Remove the first edge from each loop
  def remove(G, loops=None):
    if loops is None:
      loops = nx.cycle_basis(G)
    for l in loops:
      s1, s2 = l[0], l[1]
      if (s1, s2) in G.edges():
        G.remove_edge(s1, s2)
      elif (s2, s1) in G.edges():
        G.remove_edge(s2, s1)
    loops2 = nx.cycle_basis(G)
    print('New graph has %i loops' %len(loops2))
    return G
  if G is None:
    G, edges, nodes = as_nx_object(infile)
  # Check how many loops there are
  loops = nx.cycle_basis(G)
  # An undirected tree/graph should have n-1 edges; here there are
  print('%s has %i nodes and %i edges' %(infile, len(nodes)*2, len(edges)))
  print('Should have %i edges' %(len(nodes)*2-1))
  print('Predicts ~ %i loops' %(len(edges)-(len(nodes)*2-1)))
  print('Actually found %i loops' %len(loops))
  # Copy nodes and edges first
  edges2, nodes2 = [i for i in edges], [j for j in nodes]
  # Try to remove all the loops, but with some tolerance in case a wall
  same = 0
  while len(loops) > 0 and same < 100:
    G = remove(G, loops)
    loops2 = nx.cycle_basis(G)
    if len(loops) == len(loops2):
      same = same + 1
    else:
      loops = [l for l in loops2]
  return G



def rewrite_hoc(infile, outfile=None, G=None):
  """
  Use networkx to remove the loops and then save the new hoc file.
  """
  # Turn each edge into a hoc-friendly dictionary
  def add_connection(edge):
    s1, s2 = edge[0], edge[1]
    seg1, loc1 = s1.split('.')[0], s1.split('.')[1]
    seg2, loc2 = s2.split('.')[0], s2.split('.')[1]
    if seg1 != seg2:
      return {int(seg1): int(loc1), int(seg2): int(loc2)}
    else:
      return {}
  # Write the dictionary into the hocfile; keep the rest the same
  def write_hoc(infile, outfile, connections):
    with open(infile, 'r') as fIn:
      with open(outfile, 'w') as fOut:
        for line in fIn:
          if line:
            splitLine = line.split(None)
            if type(splitLine) is list and len(splitLine) > 0:
              if splitLine[0] == 'connect': # Ignore the connections
                pass 
              else:
                fOut.write(line)
            else:
              fOut.write(line)
    with open(outfile, 'a') as fOut: # Open the outfile again for appending
      fOut.write('\n')
      for c in connections:
        fOut.write('connect filament_999[%i](%i), filament_999[%i](%i)\n'
                   %(list(c.keys())[0], c[list(c.keys())[0]],
                     list(c.keys())[1], c[list(c.keys())[1]]))
    return
  #
  if G is None:
    G = remove_loops(infile)
  else:
    G = remove_loops(infile, G=G)
  if outfile is None:
    outfile = infile.split('.')[-2]+'_NoLoops.hoc'
  connections = []
  gedges = G.edges()
  print('Num edges: %i' %(len(gedges)))
  for e in gedges:
    E = add_connection(e)
    if len(E) == 2:
      connections.append(E)
  write_hoc(infile, outfile, connections)
  return G
  
        
  

########################################################################
if __name__ == "__main__":
  args = sys.argv
  hocIn = args[1]
  if len(args) > 2:
    hocOut = args[2]
  else:
    hocOut = None
  rewrite_hoc(hocIn, hocOut)


















