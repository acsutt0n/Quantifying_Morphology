## hocOrders.py -- for a neuron with max branch order N, create
#                  N/10 hoc files, each with 10 steps of increasing
#                  complexity
# usage: python hocOrder.py infile.hoc outfile.hoc (optional)


def get_seglist(geo):
  """
  Return the seglist and the PathDistanceFinder object.
  """
  pDF = PathDistanceFinder(geo, geo.soma)
  seglist = list(pDF.branchOrders.keys())
  return pDF, seglist



def get_ordern(hocfile, geo, n, outfile):
  """
  Only return segments of order n and below.
  """
  geo.calcBranchOrder()
  keepInds = [seg.filamentIndex for seg in geo.segments
              if seg.branchOrder <= n]
  print('Soma is %s (fil index: %i)' %(geo.soma.name, geo.soma.filamentIndex))
  
  # Now write only these segments (& soma)
  with open(hocfile, 'r') as fIn:
    with open(outfile, 'w') as fOut:
      for line in fIn:
        splitLine = line.split(None)
        
        




def write_hoc(geofile, segarray, outfile=None):
  """
  This writes the hoc file.
  """
  if outfile is None:
    outfile = geofiles.split('.')[-2]+'_order.hoc'
  with open(outfile, 'w') as fOut:
    
    prevseg = None
    for s in range(len(segarray)):
      elem = segarray[s]
      if elem[-1] != prevseg: # Start a new segment
      
  
  


def get_orders(geo, orders):
  """
  
  """
  
  pDF, seglist = get_seglist(geo)
  





















