# create a list of points with all Hoc info for each point 
#  (for use in matlab)
#  format: segNum, x, y, z, hocRad
#######
# usage: python hocPointsSections.py hocFileName pointsFileName

def hocToPoints(hocFile):
  
  outfile = open(pointFile, 'w')
  current_section = -1
  
  # write line to point file
  def writeline(newPoint):

    if type(newPoint) is list: # if newPoint exists (if parseHoc returned something
      # other than None)
      print('writing newPoint')
      print(newPoint)
      ptString = ['0','0','0','0','0']
      
      for c in range(len(newPoint)):
        ptString[c] = str(newPoint[c])
      printPtLine = ' '.join(ptString)
      outfile.write(printPtLine)
      outfile.write('\n')
  
    else:
      pass


  # function order -- open the infile (hoc) and parse it
  with open(hocFile, 'r') as infile:
    for line in infile:
      newPoint = parseHoc(line, current_section)
      if type(newPoint) is int and newPoint >= 0:
        current_section = newPoint
      writeline(newPoint)
      print(newPoint)
      
      





def parseHoc(line, current_section):
  cols = line.strip().split('_')
  
  # part of connection matrix; skip
  if cols[0] == 'connect section': 
    current_section = -1
    return current_section
  
  # found new section
  elif cols[0] == 'section':
    newcols = cols[1].split('{')
    current_section = int(newcols[0])
    print('found section %i' % current_section)
    return current_section
  
  # if it's the end of the section, close it
  elif cols[0] == '}':
    current_section = -1
    return current_section
  
  # else it must be a pt3d line
  else:
    pt3dcols = cols[0].split('(')
    
    # ptd3dclear line
    if pt3dcols[0] == 'pt3dclear':
      return -1
    
    # pt3dadd line
    elif pt3dcols[0] == 'pt3dadd': # and current_section >= 0:
      print('adding new point')
      pointcols = pt3dcols[1].split(',')
      x = float(pointcols[0])
      y = float(pointcols[1])
      z = float(pointcols[2])
      radcols = pointcols[3].split(')')
      rad = float(radcols[0])
      newPoint = [current_section, x, y, z, rad]
      return newPoint
    
    elif pt3dcols[0] == '\n':
      return -1
    
    # line not recognized
    else:
      return -1
  



#####################################################
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]
  print('using %s' % hocFile)
  pointFile = sys.argv[2]
  hocToPoints(hocFile)

