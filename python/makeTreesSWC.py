# create a SWC file readable by TREES (N x 6)
### usage: python makeTreesSWC.py hocFile swcFile somaSeg
# [ node# regiontype  X  Y  Z  rad  parID ]



def makeTreesSWC(hocFile, swcFile, somaSeg):
  import numpy as np
  
  outfile = open(swcFile, 'w')
  openSection = 0 # this tells the function whether it should be adding
                  # the next point to the section and skip other parsing
  
  # write line
  def writeSWCline(current_section, iden, xmean, ymean, zmean, radmean):
    # find parent
    if current_section == int(somaSeg): # use manually identified root section
      ind = -1
    else:
      try:
        ind = zerothEnd.index(current_section)
      except Exception as e:
        print('section %i has no parent' %current_section)
        # print(e)
        ind = -2
    
    # identify soma
    if ind >= 0:
      parent_section = onethEnd[ind]
    # bad section, give it an artificial parent and tell the user
    elif ind == -2:
      print('presumed parent of section %i' %current_section)
      parent_section = current_section-1
    # if root (mark as soma)
    elif ind == -1:
      print('root section is: %i' %int(rootSeg))
      parent_section = -1
      iden = 1;
    
    # identify axon
    if current_section == int(axonSeg):
      iden = 2
      print('axon section is: %i' %int(axonSeg))
    
    parent_section = int(parent_section)
    currentSWCline = [current_section, iden, xmean, ymean, zmean, \
    radmean, parent_section]
    
    # LMeasure starts counting segments at ID = 1
    currentSWCline[0] = currentSWCline[0] +1
      
    SWCstring = ['0','0','0','0','0','0','0']
    for c in range(len(currentSWCline)):
      SWCstring[c] = str(currentSWCline[c])
    
    printSWCstring = ' '.join(SWCstring)
    # SWCstring = str(currentSWCline)
    outfile.write(printSWCstring)
    outfile.write('\n')
    # print('SWC line written')
    
    
  

                 
  # get connection data first
  with open(hocFile, 'r') as infile:
    zerothEnd = [];
    onethEnd  = [];  # set connection arrays to 0
    for line in infile:
      cols = line.strip().split('_')
      # turns 'connect section_355(0), section_301(1)'
      # into ['connect section', '352(0), section', '349(1)']

      if cols[0] == 'connect section':
        current_zeroth = cols[1].split('(')
        zerothEnd.append(float(current_zeroth[0]))
        current_oneth = cols[2].split('(')
        onethEnd.append(float(current_oneth[0]))
        # print('connection read')
      else:
        continue
    
    # once entire connection arrays are known, 
    # get SWC sections
  with open(hocFile, 'r') as infile:
    Xs = Ys = Zs = Rads = 0  # initialize 
    for line in infile:
      cols = line.strip().split('_') # start the same way
      if cols[0] == 'connect section':
        continue # skip this part now
        
      elif cols[0] == 'section':
        newcols = cols[1].split('{')
        current_section = int(newcols[0])
        openSection = 1
        # print('found section')
      
      # if end of a section is found, write section data
      elif cols[0] == '}' and openSection:
        openSection = 0
        # print(current_section)
        xmean = sum(Xs) / len(Xs)
        ymean = sum(Ys) / len(Ys)
        zmean = sum(Zs) / len(Zs)
        radmean = sum(Rads) / len(Rads)
        writeSWCline(current_section, 3, xmean, ymean, zmean, radmean)
        Xs = Ys = Zs = Rads = 0  # reset these to 0
      
      # parse lines to get points
      # turns 'pt3dadd(309.060028,174.105026,52.000000,0.090000)'
      # into ['pt3dadd', '309.060028,174.105026,52.000000,0.090000)']
      else:
        pt3dcols = cols[0].split('(')
        if pt3dcols[0] == 'pt3dclear': # skip it
          # print('pt3dclear skipped')
          continue
          
        # gotta add these points and rads if openSection==1
        elif pt3dcols[0] == 'pt3dadd' and openSection: 
          pointcols = pt3dcols[1].split(',')
          # into ['309.060028', '174.105026', '52.000000', '0.090000)']
          # coordinates
          if not Xs:
            Xs = [float(pointcols[0])]
          else:
            Xs.append(float(pointcols[0]))
          if not Ys:
            Ys = [float(pointcols[1])]
          else:
            Ys.append(float(pointcols[1]))
          if not Zs:
            Zs = [float(pointcols[2])]
          else:
            Zs.append(float(pointcols[2]))
          
          # for radius
          # and from ['0.090000)'] into ['0.090000']
          radcols = pointcols[3].split(')')
          if not Rads:
            Rads = [float(radcols[0])]
          else:
            Rads.append(float(radcols[0]))
          # print('point added')
        elif pt3dcols[0] == '\n':
          continue
        else:
          continue # print('line not recognized')







#####################################################
if __name__ == '__main__':
  import sys
  hocFile = sys.argv[1]
  swcFile = sys.argv[2]
  somaSeg = sys.argv[3]

  makeTreesSWC(hocFile, swcFile, somaSeg)

