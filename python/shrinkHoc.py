# shrinkHoc.py -- create new hocs that are different sizes, they "grow"
# usage: python shrinkHoc.py hocfile num_of_growths 
#                                         ^
#                  number of steps to go from soma only to full neuron


import numpy as np
from neuron_readExportedGeometry import *
import os



def order_hoc(geo, order, root_name=None):
  """
  Get all the branches from the soma to this order of branch.
  But am using segment orders (max(segOrder) = 1/2(max(branchOrder)))
  """
  def write_seg(seg, sec_names, filament_root, fOut):
    # Write a segment to a new hoc
    fOut.write('%s[%i] {\n' %(filament_root, sec_names.index(seg.name)))
    fOut.write('  pt3dclear() \n')
    for n in seg.nodes:
      fOut.write('  pt3dadd(%.7f, %.7f, %.7f, %.6f)\n' 
                 %(n.x, n.y, n.z, n.avgRadius))
    fOut.write('}\n\n')
    return
  #
  def write_connections(geo, sec_names, filament_root, fOut):
    # Write only the connections that are between included segments
    for con in geo.connections:
      if con['filament1'] in sec_names and con['filament2'] in sec_names:
        sec1 = str('%s[%i]' %(filament_root, sec_names.index(con['filament1'])))
        sec2 = str('%s[%i]' %(filament_root, sec_names.index(con['filament2'])))
        fOut.write('connect %s(%i), %s(%i) \n' 
                   %(sec1, int(con['location1']),
                     sec2, int(con['location2'])))
    return
  #
  # Get a list of segments to keep
  keep_inds = [s.filamentIndex for s in geo.segments
                                   if s.branchOrder <= order]
  sec_names = [s.name for s in geo.segments if s.branchOrder <= order]
  # Use the index of sec_names for the count
  # Get the filament name
  filament_root = sec_names[0].split('[')[0]
  # Create the new hocfile
  if root_name is None:
    root_name = geo.name
  fname = root_name + '_order-' + str(order) + '.hoc'
  # Write the hoc with all branches up to branchOrder=order
  with open(fname, 'w') as fOut:
    # If a header is needed place it here
    fOut.write('create %s[%i] \n' %(filament_root, len(sec_names)))
    count = 0
    for s in geo.segments:
      if s.name in sec_names:
        write_seg(s, sec_names, filament_root, fOut)
        count = count + 1
    write_connections(geo, sec_names, filament_root, fOut)
  print('%s (order %i) written!' %(fname, order))
  return



def get_orders(geo, orders=10, thing='segs'):
  """
  Get how many (and which) orders for a new series of hocfiles.
    If thing=='int' you'll get back that many hocs, i.e.: splits the
      branchOrders by 'orders'. 
    If thing=='segs' each subsequent hoc will add all the branches
      with orders >= previous_order + 'orders'; so will get back
      max(branchOrders)/'orders'
    If 'orders' is a list, it just uses those as the orders and thing is ignored.
  """
  maxOrder = max([s.branchOrder for s in geo.segments])
  if type(orders) is not list:
    # Get back 'orders' number of hocs
    if thing == 'int':
      orders = [int(i) for i in np.linspace(1,maxOrder, orders)]
    # Get back branchOrder/'orders' number of hocs
    elif thing == 'segs':
      orders = [int(i) for i in np.linspace(1, maxOrder, int(maxOrder/orders))]
    else:
      print('Orders (arg 2) should be int or list!'); return None
  # Create all the hoc files
  for o in orders:
    order_hoc(geo, o, geo.name)
  return




#########################################################################

if __name__ == "__main__":
  args = sys.argv
  geo = demoReadsilent(args[1])
  if len(args) == 3:
    orders = args[2]
  else:
    orders = 10
  if len(args) == 4:
    things = args[3]
  else:
    things = 'segs'
  get_orders(geo, orders, things)
  











