# removeIsolatedSections.py -- remove isolated sections (connected to '-1')

def doShit(fName):
  with open(fName, 'r') as infile:
    pLines = [parse(line) for line in infile if not isStupid(line)]
    fLines = [line for line in infile]
  
  



if __name__ == '__main__':
  import sys
  fName = sys.argv[1]
  doShit(fName)




from sys import argv

script, in_file, out_file = argv

infile = open(in_file, 'r')

# isolated sections to be removed
iso_sections = [];

# write a line to the output file
outfile = open(out_file, 'w')

def lineout (line, omit):
  if omit == 1:
    omit = omit # do nothing
  else:
    outfile.write(line)
    print('line written')

# run through lines and focus only on the connection instructions
for line in infile:
  # print(line)
  cols = line.split('_')
  
  if cols[0] == 'connect section':
    isoSection = cols[2].split('(')
    #print(isoSection)
    if isoSection[0] == '-1':
      isoSecNum = cols[1].split('(')
      newIsoSecNum = int(isoSecNum[0])
      iso_sections.append(newIsoSecNum)

# jump to top of file
infile.seek(0)

# remove isolated sections and isolated connections from hoc file
for i, line in enumerate(infile):
  cols = line.split('_')
  print i
  # if a new section is being created, check if it's an isolated section
  if cols[0] == 'create section': 
    sec = int(cols[1])
    if sec in iso_sections:
      omit = 1
    else:
      omit = 0
      lineout(line, omit)
      
  elif cols[0] == 'section':
    lineout(line, omit)
  elif cols[0] == 'connect section':
    isoSection = cols[2].split('(')
    # if trying to connect an isolated section
    if isoSection[0] == '-1':
      lineout(line, 1) # omit = 1
    else:
      lineout(line, 0) # omit = 0
  elif cols[0] == '}':
    lineout(line, omit)
  elif cols[0] == '\n':
    outfile.write('\n')
    
  # else, must be a pt3dclear or pt3dadd
  else:
    pt3d = cols[0].split('(')
    if pt3d[0] == 'pt3dclear':
      lineout(line, omit)
      print('Found new segment')
    elif pt3d[0] == 'pt3dadd':
      lineout(line, omit)
    else:
      print('Error trying to write line %i ' %i)

# print(iso_sections)
  
infile.close()
outfile.close()

    
    




