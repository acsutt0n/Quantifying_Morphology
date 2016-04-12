# remove all multiple blank lines

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


###########
if __name__ == '__main__':
  import sys
  args = sys.argv
  remove_carriage_returns(args[1], args[2])

