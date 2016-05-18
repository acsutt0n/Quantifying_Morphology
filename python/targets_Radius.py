# targets_Radius.py

"""
This includes all radius stuff -- taper, Rall, fitting, etc.
"""


from target_Paths import *




def load_params(addpath='/home/alex/code/morphology/python/build-morphology/'):
  """
  Load the parameters for the fits. These are derived in Cuntz... Segev 
  (2007). For each segment of normalized length ldend[i], P[i] are the
  parameters for  y = P[i,0]x^2 + P[i,1]x + P[i,2] = 
  """
  P = []
  with open(addpath+'P.txt','r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        P.append([float(i) for i in splitLine])
  
  # All ldend values are tab or space-separated on one line
  ldend = []
  with open(addpath+'ldend.txt', 'r') as fIn:
    for line in fIn:
      if line:
        splitLine = line.split(None)
        for s in splitLine:
          ldend.append(float(s))
  return P, ldend



def get_quad_targets(geo, diamfile, primNeurDiam=15.4, 
                     minTip=0.0001, show=False, ):
  """
  For a simpler analysis, this just fits the paths to a quad taper
  and returns what the diameter should be at the targets.
  Targets are (x,y,z) midpoints. Filaments are filInds (ints).
  This requires adriane2 to be already loaded!
  """
  def solve(params, x):
    return params[0]*x**2 + params[1]*x + params[2]
  
  # Diagnostics
  pDF = PathDistanceFinder(geo, geo.soma)
  tipInds, tipLocs = geo.getTipIndices()
  tipSegs = []
  for i in tipInds:
    tipSegs.append([s for s in geo.segments if s.filamentIndex==i][0])
  paths = [pDF.pathTo(tipSegs[u], tipLocs[u]) for u in range(len(tipInds))]
  pathlengths = [pDF.distanceTo(tipSegs[u], tipLocs[u])
                for u in range(len(tipInds))]
  # Look @ paths
  print('Found %i paths ranging from %.2f to %.2f' 
        %(len(pathlengths), min(pathlengths), max(pathlengths)))
    
  # Now check the targets
  targ_stuff = []
  temp_targ = getDiameters(geo, diamfile, show=None)
  for t in temp_targ.keys():
    # Find the closest segment first -- this is only efficient with small geofiles
    closestSeg = [s for s in geo.segments if 
                  s.filamentIndex==temp_targ[t]['closestFil']][0]
    targDist = temp_targ[t]['pathLength']
    try:
      whichPath = paths.index([p for p in paths if closestSeg in p][0])
    except: # The target is suspected to be on an axon!
      print('Target is on an axon; adding a new path.')
      paths.append(pDF.pathTo(closestSeg, 1.)) # Assume 1. Doesn't make a big difference w/ knossos (short segments)
      pathlengths.append(pDF.distanceTo(closestSeg, 1.)) 
      whichPath =  int(len(paths)-1) # Len-1 (newest added path)
    targ_stuff.append([closestSeg, whichPath, targDist])
  
  # Look @ targets
  print('Found %i targets ranging from %.2f to %.2f'
        %(len(targ_stuff), min([t[2] for t in targ_stuff]), 
          max([t[2] for t in targ_stuff])))
  
  # Fit the longest path to the longest ldend[-1] value 
  #   and the shortest path to ldend[0], interpolate between
  P, ldend = load_params()
  norm_pathlengths = [p/max(pathlengths) for p in pathlengths] # 0-1
  print([t[1] for t in targ_stuff])
  norm_targdists = [t[2]/pathlengths[t[1]] for t in targ_stuff] # 0-1 # Norm. by its respective path
  path_param_inds = [int(p*len(ldend)-1) for p in norm_pathlengths] # ldend[ind]
  targ_param_inds = [path_param_inds[u[1]] for u in targ_stuff] # ldend[ind]
  
  # Fit 100 pts for each target's path; divide by ldend[0]*P values
  #   and multiply by initial_xy/z
  
  path_with_targets, Xs = [], [] # Interpolate 100 points for the x
  for t in range(len(targ_param_inds)):
    temp_path = None
    temp_xs = np.linspace(0.,ldend[targ_param_inds[t]],100)
    temp_path = [solve(P[targ_param_inds[t]], x) 
                 for x in temp_xs]
    temp_path = [(p/temp_path[0])*primNeurDiam for p in temp_path]
    path_with_targets.append(temp_path)
    Xs.append([i/max(temp_xs)*pathlengths[targ_stuff[t][1]] for i in temp_xs])
  fit_diams = [path_with_targets[p][len([i for i in Xs[p] if i < targ_stuff[p][2]])]
               for p in range(len(path_with_targets))]
  
  # Show the target paths
  if show:
    for p in range(len(path_with_targets)):
      col = np.random.random(3)
      plt.plot(Xs[p], path_with_targets[p], color=col)
      plt.plot(targ_stuff[p][2], 
               path_with_targets[p][len([i for i in Xs[p] if i < targ_stuff[p][2]])],
               's', color=col)
      plt.xlabel('Distance (um)')
      plt.ylabel('Fitted XY-diameter (um)')
      plt.title(geo.name.split('GM')[0])
    plt.show()
  return fit_diams, [p[2] for p in targ_stuff], Xs, path_with_targets

#



def show_tapers(geo, diamfile, primNeurDiam=15.4):
  """
  """
  return







########################################################################

if __name__ == "__main__":
  print('Module is used interactively.')




