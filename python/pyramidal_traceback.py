# this removes all of the redundant sections, catching tracebacks
# into a file and parsing them to figure out which sections to remove
# temp out file = temp_out.txt
# usage: python pyramidal_traceback.py hocfile 


import sys, traceback
from pyramidal_readExportedGeometry import *



########################################################################
# xmlwhich

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



def write_new(listoflines, newfilename='temphoc.hoc'):
  # write the lines list to a new file
  print('Writing new hoc file %s' %newfilename)
  with open(newfilename, 'w') as fOut:
    for l in listoflines:
      fOut.write(l)
      fOut.write('\n')
  return



def xmlwhich_control(hocfile, seg0, seg1):
  lol = remove(hocfile, seg0, seg1)
  write_new(lol)
  return


#######################################################################
# remove_carriage_return


def remove_carriage_returns(filein, fileout):
  lol = []
  with open(filein, 'r') as fIn:
    lineNum, carrcount = 0, 0
    for line in fIn:
      lineNum = lineNum + 1
      if line == '\n':
        carrcount = carrcount + 1
      else:
        carrcount = 0
      if carrcount < 2:
        lol.append(line)
  print('File %s was %i lines long' %(filein, lineNum))
  with open(fileout, 'w') as fOut:
    for l in lol:
      fOut.write(l)
  
  print('New file %s created is %i lines long' %(fileout, len(lol)))
  return


######################################################################
# traceback stuff


def call_demoRead(hocfile):
  try:
    geo = demoRead(hocfile)
    return True
  except:
    print('Logging traceback')
    tb = traceback.format_exc()
    return tb



######################################################################
# loop through removing redundant sections until none exist

def parse_traceback(tb):
  splitTb = [str(k) for k in tb.split(None)]
  if splitTb[0] == 'Traceback':
    if splitTb[-1] == 'connected':
      seg1 = int(splitTb[-7].split('[')[1].split(']')[0])
      seg0 = int(splitTb[-5].split('[')[1].split(']')[0])
      return [seg0, seg1]
    else:
      print('UNEXPECTED TRACEBACK!')
      print(tb)
      return [None]
  else:
    print(splitTb[0])
    return [None]


def remove_redundant(hocfile):
  """
  Hocfile is a str, not file or openfile obj.
  """
  tempout = 'temp_outfile.txt'
  # Run first pass
  tb = call_demoRead(hocfile)
  if tb is True or tb.split(None)[0] == 'None':
    return # file is already ready
  else:
    go = True
    
  # Loop until the file is ready  
  while go:
    print('Got new traceback')
    seglist = parse_traceback(tb)
    if len(seglist) == 1:
      return # file is ready
    elif len(seglist) == 2:
      # two segments returned -- remove them
      xmlwhich_control(hocfile, seglist[0], seglist[1])
      # remove carriage returns
      remove_carriage_returns('temphoc.hoc', hocfile)
      # check new hocfile
      tb = call_demoRead(hocfile)
    else:
      print('bad seglist:')
      print(seglist)
      return 
    if tb is True or tb.split(None)[0] == 'None':
      return # file is ready
    
  # end of while loop and function
  


if __name__ == '__main__':
  hocfile = sys.argv[1]
  remove_redundant(hocfile)
  #sys.exit(1)
    
