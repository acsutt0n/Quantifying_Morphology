# fake_turn.py 
# usage: python fake_turn.py img_dir

# Fake-turn filament, a fake filament that turns in 3-D. 


import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os



def load_filament(img_dir, makebinary=False, keepfirst=True, switch=True):
  """
  Best if labels are like '01' instead of '1' so it can be sorted easier.
  makebinary: 0's and 1's
  switch: inverts colors
  """
  print('Loading images...')
  fils = os.listdir(img_dir)
  fils.sort()
  fils = [img_dir+f for f in fils if f.split('.')[-1].lower() == 'tif']
  arr = []
  for f in fils:
    if fils.index(f) % 20 == 0:
      print('%i (of %i) images loaded.' %(fils.index(f), len(fils)))
    img = Image.open(f)
    arr.append(np.array(img))
  if keepfirst is True:
    varr = keep_first_color(arr)
    arr = varr
  if switch is True: # This assumes range of 0-255
    print('Switch is True; assumes a range of 0-255!')
    sarr = []
    for a in arr:
      t = []
      for i in range(len(a)):
        j = [abs(x-255) for x in a[i]]
        t.append(j)
      sarr.append(t)
    arr = sarr
  return arr



def get_skeleton(arr, target=0):
  """
  Returns the skeleton of the filament, taken as the point where the
  row has the most 0's.
  """
  skelpts = []
  for a in range(len(arr)):
    ind, cnt, mid = None, 0, None
    for row in range(len(arr[a])):
      if list(arr[a][row]).count(target) > cnt:
        cnt = list(arr[a][row]).count(target)
        ind = row
        mid = list(arr[a][row]).index(target) + int(cnt/2) # Middle of thickest segment
    skelpts.append([a, ind, mid])
  # Now have all the skelpoints
  return skelpts
  


def keep_first_color(arr):
  """
  Sometimes tifs have 4 values (R,G,B,intensity). This keeps just first.
  """
  print('Romoving all but first color values...')
  varr = []
  for a in arr:
    t = []
    for i in range(len(a)):
      j = [x[0] for x in a[i]]
      t.append(j)
    varr.append(t)
  return varr


  
  
  
  








