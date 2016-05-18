# targets_Dendrograms.py

"""
Creates rectangular and circular dendrograms with and without axon
identifiers.
"""


from target_Paths import *


########################################################################
# Dendrogram from tips (square)

"""
Make a dendrogram using a geo object (known connectivity, branches)
"""
class ConnDendrogram:
  def __init__(self, geo):
    """ Requires a geo object as input. """
    self.geo = geo
    self.axon_ids = [] # 
    self.Z = [] # Linkage matrix
    
    self.make_linkage()
    # end init
  
  
  def getBranchTips(self):
    """ Get the branch tips (indices are returned as segments). """
    tipInds, tipLocs = self.geo.getTipIndices()
    # Populate the leaf dictionary -- easy access to branches
    tipSegs = [seg.name for seg in self.geo.segments if 
               seg.filamentIndex in tipInds]
    bTipsIdx = []
    for tip in tipSegs:
      bTipsIdx.append([self.geo.branches.index(b) for b in self.geo.branches
                       if tip in b.tags][0])
    leafDict = {}
    for tipI in bTipsIdx: # Leaf is key, branch obj is value
      leafDict[bTipsIdx.index(tipI)] = self.geo.branches[tipI]
    
    return leafDict
  
  
  def make_linkage(self):
    """
    Do everything to make the linkages. A cluster combines two branches
    or two branch sets.
    """
    Z, tolerance = [], 0
    # Get tips
    leafDict = self.getBranchTips()
    contains = [1 for i in leafDict.keys()] # Each tip contains 1 tip
    branchDict = {}
    
    for k in leafDict.keys(): # Branch obj is key, Leaf is value
      branchDict[leafDict[k]] = k
    activeLeafs = list(range(len(leafDict))) # At first, all are active
    newLeaf = len(activeLeafs) # This is _n_ at first - the number of observations
    # New leafs will be added starting at number newLeaf
    
    # Iterate through all active leafs
    while len(activeLeafs) > 0:
      to_pop = []
      to_push = []
      print('Iterating through %i active leaves...' %len(activeLeafs))
      for act in activeLeafs: # Try to add the leaves in activeLeafs
        print('Working on leaf %i' %act)
        nebs = leafDict[act].neighbors # Get the neighbors
        for ne in nebs:
          if ne in branchDict.keys(): # Make sure this key exists
            # If the branch is also active (its downstream
            if branchDict[ne] in activeLeafs and ne not in to_pop: 
                                              # branches and tips are already logged)
              Z.append([act, branchDict[ne], 1, (contains[act]+contains[branchDict[ne]])])
            
              # Pop them both (act, ne) and activate other nebs
              to_pop.append(act)
              to_pop.append(branchDict[ne])
            for new_neighbor in nebs:
              if new_neighbor not in to_pop and \
                 new_neighbor not in branchDict.keys(): # Update everthing!
                   to_push.append(newLeaf)
                   # activeLeafs.append(newLeaf) # Add this later (don't change while iter.)
                   branchDict[new_neighbor] = newLeaf
                   leafDict[newLeaf] = new_neighbor
                   contains.append(contains[act]+contains[branchDict[ne]])
                   newLeaf += 1 
      
      # All active leaves are iterated through once; pop and push
      activeLeafs = [a for a in activeLeafs if a not in to_pop]
      for pu in to_push:
        activeLeafs.append(pu)
      if to_pop == 0 and to_push == 0: # Tolerance
        tolerance += 1
      else:
        tolerance = 0
      if tolerance >= 10:
        print('Reached a wall'); break
    
    # In theory, everything has been updated and now activeLeafs is empty
    self.Z = Z
    return self
  
  
# def plot_dend
    


########################################################################
# Testing


for act in activeLeafs:
    nebs = leafDict[act].neighbors
    for ne in nebs:
        if ne in branchDict.keys():
            if branchDict[ne] in activeLeafs and ne not in to_pop:
                Z.append([act, branchDict[ne], 1, contains[act]+contains[branchDict[ne]]])
                to_pop.append(act)
                to_pop.append(branchDict[ne])
            for new_neighbor in nebs:
                if new_neighbor not in to_pop and new_neighbor not in branchDict.keys():
                  to_push.append(newLeaf)
                  branchDict[new_neighbor] = newLeaf
                  leafDict[newLeaf] = new_neighbor
                  contains.append(contains[act]+contains[branchDict[ne]])
                  newLeaf = newLeaf + 1














########################################################################

if __name__ == "__main__":
  print('Module is used interactively.')


