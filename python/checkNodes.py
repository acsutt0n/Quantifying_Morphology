## check if any nodes are backwards

def checkNodes(hocFile):
  
#  class ConnectionData:
#    def __init__(self):
#      self
  
  
  
  ################# parse hoc file
  def parseHoc(hocFile):
    
    
    # -------- parse connection data ----------
    def parseConnections(hocFile):
      # get connection data first
      with open(hocFile, 'r') as infile:
        zerothEnd = [];
        onethEnd  = [];  # set connection arrays to 0
        for line in infile:
          cols = line.strip().split('_')

          if cols[0] == 'connect section':
            current_zeroth = cols[1].split('(')
            zerothEnd.append(float(current_zeroth[0]))
            current_oneth = cols[2].split('(')
            onethEnd.append(float(current_oneth[0]))
            # print('connection read')
          else:
            continue
      return zerothEnd, onethEnd
    
    # zero_nodes and one_nodes
    def parseSec
    with open(hocFile, 'r') as infile:
      Xs = Ys = Zs = Rads = 0  # initialize 
      for line in infile:
        cols = line.strip().split('_') # start the same way
        if cols[0] == 'connect section':
          continue # skip this part now
          
        elif cols[0] == 'section':
          newcols = cols[1].split('{')
          current_section = int(newcols[0])
          openSection = 0 # keep this as 0 until the first point is found
          # print('found section')
        
        # if end of a section is found, create node index
        elif cols[0] == '}' and openSection:
          openSection = 0
          print(current_section)
          writeNodes(zero_node, one_node, zeros, ones)
          zero_node = one_node = 0 # then reset the nodes

        
        # parse pt3d points
        else:
          pt3dcols = cols[0].split('(')
          if pt3dcols[0] == 'pt3dclear': # skip it
            # print('pt3dclear skipped')
            continue
            
          # add points
          elif pt3dcols[0] == 'pt3dadd': 
            pointcols = pt3dcols[1].split(',')
            # also get radius
            radcols = pointcols[3].split(')')
            if not zero_node: # if zero node isn't defined, create it
              zero_node = [float(pointcols[0])]
              zero_node.append(float(pointcols[1]))
              zero_node.append(float(pointcols[2]))
              zero_node.append(float(radcols[0]))
            else: # zero node already defined, write/rewrite one node
              one_node = [float(pointcols[0])]
              one_node.append(float(pointcols[1]))
              one_node.append(float(pointcols[2]))
              one_node.append(float(radcols[0]))

          elif pt3dcols[0] == '\n':
            continue
          else:
            continue # print('line not recognized')
    
    return zeros, ones, zerothEnd, onethEnd
  
  
  ################# writeNodes to array
  def writeNodes(zero_node, one_node, zeros, ones):
    
    # zeros and ones
    if not zeros and not ones:
      zeros = {current_section: zero_node}
      ones = {current_section: one_node}
    else:
      zeros[current_section] = [zero_node]
      ones[current_section] = [one_node]
    return zeros, ones
  
  
  
  #################### control
  zeroes, ones, zerothEnd, onethEnd = parseHoc(hocFile)
  # cycle through the nodes
  
  # sanity test
  numNodes0 = len(zerothEnd)
  numNodes1 = len(onethEnd)
  numZeros = len(zeros)
  numOnes = len(ones)
  
  if numNodes0 != numNodes1:
    print('Connection matrix mismatch')
  elif numZeros != numNodes0:
    print('No. of zero nodes doesn"t match zeroth connection matrix')
  elif numOnes != numZeros:
    print('No. of nodes mismatch')
  elif numOnes != numNodes1:
    print('No. of one nodes doesn"t match oneth connection matrix')
  else:
    print('All matrices match')
  
  
  for node in range(len(zerothEnd)):
    
    goodConnect = 0 # initialize counter
    badConnect = 0
    
    # mismatch
    if zeros[zerothEnd[node]] == ones[onethEnd[node]]:
      goodConnect = goodConnect + 1
    else:
      badConnect = badConnect + 1
    
  print('Good segment connections: %i' % goodConnect)
  print('Bad segment connections: %i' % badConnect)
  
  
  
  
  
  
##################### start of control
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]
  checkNodes(hocFile)
  
  
