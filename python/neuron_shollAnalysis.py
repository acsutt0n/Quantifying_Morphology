# needs to be imported into python to run at the moment -- not stand-alone

  def shollAnalysis(geo, straightenNeurites=True):
    """
    Find the number of neurites that intersect a sphere of a given radius
    """
    # get the centroid of the soma, weighting each compartment's contribution
    # by volume
    
    # define how distance from centroid to compartment is measured
    if straightenNeurites:
      # get distance traveled along neurites (e.g. as though neuron was
      # straightened out
      
      centroid = geo.soma.centroidPosition(mandateTag='Soma')
      
      # compute distance from soma to each segment
      somaPaths = PathDistanceFinder(geo, geo.soma, centroid)
      
      # store results in an array
      distances = []
      for s in geo.segments:
        d0, d1 = somaPaths.distanceTo(s, 0.0), somaPaths.distanceTo(s, 1.0)
        if d1 < d0:
          d0, d1 = d1, d0
        
        distances.append((d0, 1))
        distances.append((d1, -1))
      
    else:
      # get euclidean distance from soma centroid to each compartment
      # (must be done compartment by compartment, because segments curve)

      centroid = geo.soma.centroid(mandateTag='Soma')
      
      # define distance from centroid to compartment
      def _centroidDist(c):
        def _tupleDist(_t1, _t2):
          return sqrt( (_t1[0] - _t2[0])**2 + \
                       (_t1[1] - _t2[1])**2 + \
                       (_t1[2] - _t2[2])**2 )

        d0 = _tupleDist(centroid, (c.x0, c.y0, c.z0))
        d1 = _tupleDist(centroid, (c.x1, c.y1, c.z1))
        if d0 <= d1:
          return d0, d1
        else:
          return d1, d0
      
      # compute the distance from the soma centroid to each compartment that's
      # not in the soma
      distances = []
      for c in geo.compartments:
        if 'Soma' in c.tags:
          continue
        d0, d1 = _centroidDist(c)
        distances.append((d0, 1))  
        distances.append((d1, -1))
    
    # sort the distances by increasing distance
    distances.sort(key=lambda x: x[0])
    return distances
