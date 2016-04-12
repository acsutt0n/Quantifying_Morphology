# spiral.py -- this spirals in a space-filling way

import numpy as np
# from imageMatrix import *



class Spiral(): # cross-scetion
  def __init__(self, arr, start=None, tolerance=5000):
    self.arr = np.array(arr)
    self.shape = np.shape(self.arr)
    self.area_pts = []
    self.area = 0
    self.surface = None
    self.live_pts = []
    self.dead_pts = []
    self.tolerance = tolerance
    if start == None:
      self.pos = [int(self.shape[0]/2), int(self.shape[1]/2)]
    else:
      self.pos = start
    self.start = self.pos
    self.x, self.y = None, None
    self.prev = None
    self.target = self.arr[self.pos[0], self.pos[1]]
    self.no_change = {'prev_dead_pts': None, 'prev_live_pts': None,
                      'count': 0}
    # here do stuffs
    self.spiral_control()

  ####################
  # the first group of functions help spiral()
  def set_start(self):
    # assumes only one filament is present, white on dark
    _, self.pos = clean_filament(self.arr)
    


  
  def default_dead_pts(self):
    # make non-target points dead, requires center point to be in segment
    for i in range(self.shape[0]):
      for j in range(self.shape[1]):
        if self.arr[i,j] != self.target:
          self.dead_pts.append([i,j])
    #print('Starting with %i live points and %i dead points' 
    #      %(len(self.live_pts),len(self.dead_pts)))
    return self

  
  
  def no_changes(self):
    # assess whether spiral process is done
    if self.no_change['prev_dead_pts'] is None:  # set prev dead pts
      #print('Setting self.no_change[prev_dead_pts]')
      self.no_change['prev_dead_pts'] = len(self.dead_pts)
    elif self.no_change['prev_dead_pts'] == len(self.dead_pts): # increment count
      self.no_change['count'] = self.no_change['count'] + 0.5
      #print('Incremented delta to %.1f' %self.no_change['count'])
    else:
      self.no_change['prev_dead_pts'] = len(self.dead_pts) # set new prev dead pts
      #self.no_change['count'] = 0
    
    if self.no_change['prev_live_pts'] is None:  # set prev live pts
      self.no_change['prev_live_pts'] = len(self.live_pts)
    elif self.no_change['prev_dead_pts'] == len(self.live_pts): # increment count
      self.no_change['count'] = self.no_change['count'] + 0.5
      #print('Incremented delta to %.1f' %self.no_change['count'])
    else:
      self.no_change['prev_live_pts'] = len(self.live_pts) # set new prev live pts
      #self.no_change['count'] = 0
    #print(self.no_change)
    return self.no_change['count'] # 'delta'



  def check_next(self):
    check = {'right': move_down,      #    -------->
             'left': move_up,         #    ^   x-->
             'up': move_right,        #    |       |
             'down': move_left}       #    <-------  
    next_pt = check[self.prev](self.x, self.y) # get next pt
    if next_pt not in self.dead_pts and next_pt not in self.live_pts:
      if next_pt[0] in range(0,self.shape[0]) and next_pt[1] in range(0,self.shape[1]):
        #print(check[self.prev](self.x,self.y))
        #print('Next move OK')
        # print(self.live_pts)
        return [True, check[self.prev](self.x,self.y)]
      else:
        return [False]
    else:
      #print('Next move not OK')
      return [False]
  
    

  
  def check_live_pt(self, i,j):
    # determine if point should be dead or alive
    if self.arr[i,j] != self.target:
      if [i,j] not in self.dead_pts or [i,j] in self.live_pts:
        self.dead_pts.append([i.j])
      return False
    
    live = 0
    for fn in [move_down, move_up, move_right, move_left]:
      if fn(i,j) in self.live_pts or fn(i,j) in self.area_pts:
        live = live+1 
    if live >= 1 and [i,j] not in self.dead_pts:
      # if [i,j] is connected to an area pt, add it
      # self.area_pts.append([i,j]) # add it to dead points
      if [i,j] not in self.live_pts:
        #ind = self.live_pts.index([i,j])
        self.live_pts.append([i,j]) # remove it from live points
        #print('Added [%i,%i] to area_pts and removed it from live_pts'
        #      %(i,j))
      #else:
      #  self.live_pts.append([i,j])
      return False
    else:
      return True
        
    
  
  def log_pt(self): 
    # log and check if point is dead
    #print('Trying to add point %i %i' %(self.x, self.y))
    if self.arr[self.x,self.y] == self.target: # double check
      self.area = self.area + 1
      if self.check_live_pt(self.x, self.y): # if point is alive, add it
        if [self.x,self.y] not in self.live_pts:
          self.live_pts.append([self.x,self.y])
          #print('Added point %i %i. Live points: %i' %(self.x, self.y,
          #                                           len(self.live_pts)))
      # else it is automatically added to dead_pts
  
  
  ## Not used!
  def retire_pts(self):
    # if a point is surrounded on all sides by active or dead points
    # (or both), make it dead
    # print(self.shape)
    for p in self.live_pts:
      i,j = p[0], p[1]
      if not self.check_live_pt(i,j):
        return
        #print('Live points: %i, dead points: %i, area points: %i' 
        #      %(len(self.live_pts), len(self.dead_pts), len(self.area_pts)))
    
        
  
  def spiral(self):
    # continue to spiral until an edge is hit; then try other directions
    # when all directions hit an edge, move to another point
    from itertools import cycle
    self.area = self.area + 1
    self.x, self.y = self.pos[0], self.pos[1]
    self.prev = 'up' # set fake previous move
    self.log_pt()
    moves = [move_down, move_left, move_up, move_right]
    prevs = ['right','down','left','up']
    _moves = cycle(moves)
    _prevs = cycle(prevs)
    delta = self.no_changes()
    #print(self.no_change)
    # cycle through points
    while len(self.live_pts) > 0:
      count = 0
      # self.retire_pts()
      delta = self.no_changes()
      # print('Delta is %.1f' %delta)
      # take a point until exhausted
      self.x, self.y = self.live_pts[int(np.random.random(1)*len(self.live_pts))]
      # print('Starting new point %i %i' %(self.x, self.y))
      while count < 4: # if next exists, go to it
        if self.check_next()[0]:
          # print(self.check_next())
          try:
            self.x, self.y = self.check_next()[1][0], self.check_next()[1][1]
          except:
            self.x, self.y = self.check_next()[1][0][0], self.check_next()[1][0][1]
          self.log_pt()
          count = 0 # reset count, found a new point
          
          next_move = next(_moves)
          self.prev = next(_prevs)
        else:
          next_move = next(_moves)
          self.prev = next(_prevs)
          count = count + 1
          #print('Live points: %i, dead points: %i, area points: %i' 
          #      %(len(self.live_pts), len(self.dead_pts), len(self.area_pts)))
        if delta > self.tolerance:
          return
          #if count ==4:
          #  print(count)
          # otherwise, check the others and increase count
    
  
  
  ###########################
  # the next group of functions perform analyses
  def get_surface(self):
    if self.surface == None:
      self.surface = {'xs':0, 'ys':0}
    flip = lambda u: (u+1)%2 # switch between L-R and U-D
    # for each live_pt
    xs, ys = 0, 0
    for p in self.live_pts:
      i,j = p[0],p[1]
      f = 0
      for fn in [move_down, move_right, move_up, move_left]:
        f = flip(f) # f = 1 (up, down), f = 0 (left, right)
        if fn(i,j) in self.dead_pts:
          if f == 1:
            ys = ys + 1
          elif f == 0:
            xs = xs + 1
    self.surface['xs'], self.surface['ys'] = xs, ys
  
  
  
  def spiral_control(self):
    # control the class
    self.default_dead_pts()
    self.spiral()
    self.get_surface()
    # print(self.dead_pts)
    return self
        

# directional functions
move_up = lambda x,y: [x,y-1]
move_down = lambda x,y: [x,y+1]
move_right = lambda x,y: [x+1,y]
move_left = lambda x,y: [x-1,y]





# Generator and proof-reading functions

def gen_arr(choose=1):
  if choose==0:
    ary = [[0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 0, 0],
           [0, 1, 1, 1, 1, 0],
           [1, 1, 1, 1, 1, 0],
           [0, 0, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0]]
  elif choose==1:
    ary = [[0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,1,1,1,1,0,0,0],
           [0,0,1,1,1,1,1,1,1,1,0,0],
           [1,1,1,1,1,1,1,1,1,1,1,0],
           [1,1,1,1,1,1,1,1,1,1,1,0],
           [0,1,1,1,1,1,1,1,1,1,0,0],
           [0,0,1,1,0,1,1,1,1,1,1,0],
           [0,0,0,0,0,0,0,1,1,1,1,0],
           [0,0,0,0,0,0,0,0,1,1,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0]]
  
  return ary



def reproduce_matrix(pts):
  xs, ys = [p[0] for p in pts], [p[1] for p in pts]
  new_arr = np.zeros([max(xs)+1,max(ys)+1])
  for p in pts:
    new_arr[p[0],p[1]] = 1
  return new_arr



##############################
# no inputs, run a demo of the algorithm
if __name__ == '__main__':
  arr = gen_arr()
  g = Spiral(arr)











