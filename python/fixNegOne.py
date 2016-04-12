# fix -1 in connection matrix
# usage: python fixNegOne.py fileToFix newfileName
# this used to fix bad radii also but now it does not!!!!!!

import os, sys



def readfile(brokenFileName, newFileName):
  
  with open(newFileName, 'w') as fOut:
  
    with open(brokenFileName, 'r') as fIn:
      
      for line in fIn:
        splitLine = line.split(None)
        if type(splitLine) is list and len(splitLine) >= 1:
          #print(splitLine)
          #  and splitLine[2].split('_')[1].split('(')[0] == -1
          if splitLine[0] == 'connect':
            prevsec = int(splitLine[1].split('_')[1].split('(')[0])
            cursec = int(splitLine[2].split('_')[1].split('(')[0])
            if cursec == -1:
              cursec = prevsec + 1
            connect_string = 'connect section_%i(1), section_%i(0)' \
                             % (prevsec, cursec)
            fOut.write(connect_string)
            fOut.write('\n')
            #print('Fixed bad connection')
          else:
            fOut.write(line)
        else:
          fOut.write(line)
        """
        if type(splitLine) is list:
          print(splitLine)
          try:
            first = splitLine.split('(')[0]
          except:
            if splitLine[0] != 'connect':
              first = 'monkeyshit'
              fOut.write(line)
              fOut.write('\n')
          if first == 'pt3dadd':
            #splitLine = line.split('(') 
            #if splitLine[0] == 'pt3dadd':
            comsplit = splitLine[1].split(',')
            if float(comsplit[3].split(')')[0]) <= 0:
              #comsplit = splitLine[1].split(',')
              rad = float(comsplit[3].split(')')[0])
              comsplit[3]=str(0.1)
              newline = ','.join(comsplit)
              paren = str(')')
              newline.append(paren)
              newline = ''.join(newline)
              pt3d = 'pt3dadd('
              newline = pt3d.append(newline)
              fOut.write(newline)
              fOut.write('\n')
              #print('Normal connection')
        
        
        #if line.split(None) is None:
        #  fOut.write(line)
        #  fOut.write('\n')
        """
    


############################################################
if __name__ == '__main__':
  arguments = sys.argv
  brokenFileName = arguments[1]
  newFileName = arguments[2]
  readfile(brokenFileName, newFileName)
