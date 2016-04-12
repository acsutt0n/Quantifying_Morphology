# xmlwhich.py - give it two sections and it will remove them from the 
# hoc file
# usage: python xmlwhile.py hocfile seg1 seg2 (file, int, int)
# enter the numbers in the order NeuronGeometry returns them



def is_this_the_line(splitline, seg0, seg1):
  # check if this connect line is the right one to delete
  # returns true to keep, false to delete
  zer = int(splitline[1].split('[')[1].split(']')[0])
  one = int(splitline[2].split('[')[1].split(']')[0])
  if zer == seg0 or zero == seg1:
    if one == seg1 or one == seg0:
      print('Removed it')
      return False
  return True



def remove(hocfile, seg0, seg1):
  print('Reading file %s and removing connection between %i and %i'
        %(hocfile, seg1, seg0))
  lol = []
  with open(hocfile, 'r') as fIn:
    lineNum = 0
    for line in fIn:
      splitLine = line.split(None)
      try:
        if splitLine[0] == 'connect':
          if is_this_the_line(splitLine, seg0, seg1):
            lol.append(line)
        else:
          lol.append(line)
      except:
        lol.append(line)
  return lol





def write_new(listoflines, newfilename='new_hoc.hoc'):
  # write the lines list to a new file
  print('Writing new hoc file %s' %newfilename)
  with open(newfilename, 'w') as fOut:
    for l in listoflines:
      fOut.write(l)
      fOut.write('\n')
  return




def control(hocfile, seg0, seg1):
  lol = remove(hocfile, seg0, seg1)
  write_new(lol)
  return




#################
if __name__== '__main__':
  import sys
  hocfile = sys.argv[1]
  seg1, seg0 = int(sys.argv[2]), int(sys.argv[3]) # switch order here
  control(hocfile, seg0, seg1)
