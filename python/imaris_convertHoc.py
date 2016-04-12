# 
# This program changes hocs from imaris to amira/home-brew format
# so they can be visualized in amira
# This current version retains only the first radius and eliminates the
# second/minor radius. Future versions may incorporate both rad1 and rad2
# into rad1, but currently amira can only handle rad1
#
# usage: python imaris_convertHoc.py oldHoc newHoc (op) filament(op)
# If filament should be an integer that indicates which filament to keep.
# e.g.: 100000007 for 837_123_LP


import numpy as np



def reduce_rads(hocfile, outfile=None):
  """
  Remove the second radius and write a new hoc file.
  """
  if outfile is None:
    outfile = hocfile.split('.')[0] + '_radius.hoc'
  #
  def write_onerad(line, fOut):
    splitLine = line.split(',')
    if len(splitLine) > 4: # Contains additional radii, will be omitted
      thing = splitLine[:4]
      thing[-1] = thing[-1]+')'
      fOut.write(', '.join(thing))
      fOut.write('\n')
    else:
      fOut.write(line)
    return
  #
  with open(hocfile, 'r') as fIn:
    with open(outfile, 'w') as fOut:
      for line in fIn:
        if line:
          splitLine = line.split(None)
          if len(splitLine) > 0:
            if splitLine[0].split('(')[0] == 'pt3dadd': # Found a point
              write_onerad(line, fOut)
            elif splitLine[0] == 'access': # This also throws an error
              pass
            #elif splitLine[0].split('_')[0] == 'define': #
            #  pass
            else:
              fOut.write(line)
          else:
            fOut.write(line)
  print('File %s written' %outfile)
  return


"""
def filament_reduce(hocfile, outfile, filament=100000007):
  #This deals with strange cases like filament_A[B] where both A and B
  #change, like filament_99[3], filament_99[4], filament_100[0], by 
  #keeping only the designated filament.
  def write_onerad(line, fOut):
    splitLine = line.split(',')
    if len(splitLine) > 4: # Contains additional radii, will be omitted
      thing = splitLine[:4]
      thing[-1] = thing[-1]+')'
      fOut.write(', '.join(thing))
      fOut.write('\n')
    else:
      fOut.write(line)
    return
  # If it's bad, return False, else write & return True
  def check_line(line, fOut, filament): # 
    splitLine = line.split(None)
    if len(splitLine) > 1:
      try: # Start of filament
        ch = int(splitLine[0].split('_')[1].split('[')[0])
      except:
        ch = None
      try:
        
        if ch == filament:
          fOut.write(line)
          return True
          
  # If filament matches, write it, else skip it
  go = True # Until we hit something that shouldn't be there
  with open(hocfile, 'r') as fIn:
    with open(outfile, 'w') as fOut:
      for line in fIn:
        if line:
"""
  


####################
if __name__ == "__main__":
  import sys
  args = sys.argv
  hocfile = args[1]
  if len(args) > 2:
    outfile = args[2]
    reduce_rads(hocfile, outfile)
  elif len(args) > 3:
    outfile = args[2]
    try:
      filament = int(args[3])
    except:
      print('Filament argument must be an int!')
    filament_reduce(hocfile, outfile, filament)
  else:
    reduce_rads(hocfile)
    
  
  









