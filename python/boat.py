class Boat:
  def __init__(self, settins):
    
    if settins == None:
      self.sail = [1,2,3]
      self.mast = {1:[1,2], 3:[3,4]}
      self.bow = 1
    else:
      self.sail = [s ** 2 for s in settins[s]]
      self.mast = {1:[0,0,0,0]}
      self.bow = settins
    

## control
def foo():
  bar = Boat
  attrs = vars(bar)
  print(', '.join('%s: %s' % item for item in attrs.items()))
  print(bar.sail, bar.mast, bar.bow)
  


##   MAIN 
if __name__ == '__main__':
  foo()
  
