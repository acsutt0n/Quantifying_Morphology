 pyneuron (python2) test for figures for morphology grant

import neuron
import matplotlib.pyplot as plt
import numpy as np
# most NEURON specific commands start with `neuron.h.`


########################################################################
# Create geometry for model
def set_params(Ra=35, Gleak=0.01, Eleak=-50,finit=-50,
              stim_dur=5, stim_amp=-1, stim_delay=1, run_time=10):
  P = {'Ra':Ra, 'Gleak': Gleak, 'Eleak': Eleak, 'finit': finit,
     'stim_dur': stim_dur, 'stim_amp': stim_amp, 'stim_delay': stim_delay,
     'run_time': run_time}
  return P



def init_section(L, Ra, diam, nseg=None):
  dend = neuron.h.Section()
  dend.Ra = Ra
  dend.L = L
  dend.diam = diam
  if nseg is None:
    dend.nseg = int(dend.L/10)
  else:
    dend.nseg = nseg
  return dend



def createModelGeo(geo):
  """
  Create a model with the given geometry.
  """
  nrnSegs, nrnSegNames = [], []
  # gleak = 0.001 uS/nF for soma & dendrites, gleak = 0.2 for axons
  # Ra = 60 ohm*cm (Rabbah et al., 2005; Golowalsch et al., 2009)
  segNames = []
  # Create the neuron segment objects
  for s in geo.segments:
    L = s.length
    Ra = 60
    diam = s.avgRadius*2
    nrnSegs.append(init_section(L, Ra, diam))
    segNames.append(s.name)
  # Connect all the segments
  for s in range(len(nrnSegs)):
    for n, nLoc in zip(geo.segments[s].neighbors, geo.segments[s].neighborLocations):
      nrnSegs[s].connect(nrnSegs[segNames.index(n.name)], nLoc[0], nLoc[1])
  # Setting passive parameters
  for sec in neuron.h.allsec():
    # Do with the present `sec`
    sec.insert('pas')
    # sec.Ra = props['Ra'] # already set
    # Do for each segment within `sec`:
    for seg in sec: 
      # Do with the segment `seg`:
      seg.pas.g = 0.0015 # avg value per Taylor et al., 2009
      seg.pas.e = -18 # avg value per Taylor et al., 2009
  return nrnSegs



########################################################################
# Simulate with various currents
def stepCurrent(
  stim = neuron.h.IClamp(soma(0.5))
  # Setting recording paradigm
  stim.delay = props['stim_delay']
  stim.amp = props['stim_amp']
  stim.dur = props['stim_dur']



def simple_geo(props=None, retsoma=False):
  # simulate stimulations, returns time, dend2_volt +/- soma_volt
  P = set_props()
  if props is None:
    props=P
  else:
    for k in P.keys():
      if k not in props:
        props[k] = P[k]
    for k in props.keys():
      if k not in P.keys():
        print("Don't know what the fuck %s is. Options are: " %k)
        print(list(P.keys()))
    
  # Creating the morphology
  # soma (compartment: neuron.h.Section() )
  soma = init_section(100, props['Ra'], 80, 10)
  # dendrite0
  dend_0 = init_section(200, props['Ra'], bound0)
  # dendrite1, with taper
  dend_1 = init_section(200, props['Ra'], bound0)
  diams = np.linspace(bound0, bound1, dend_1.nseg)
  rad = -1
  for seg in dend_1:
    rad = rad+1
    seg.diam = diams[rad]
  # dendrite2
  dend_2 = init_section(200, props['Ra'], bound1)

  dend_0.connect(soma, 1, 0) # connect soma(1) with dend_0(0)
  dend_1.connect(dend_0, 1, 0) # connect dend_0(1) with dend_1(0)
  dend_2.connect(dend_1, 1, 0) # connect dend_1(1) with dend_2(0)

  # Implementing a current clamp electrode
  # Locate the electrode at the center of the soma
  stim = neuron.h.IClamp(soma(0.5))
  # Setting recording paradigm
  stim.delay = props['stim_delay']
  stim.amp = props['stim_amp']
  stim.dur = props['stim_dur']
  


  # Record Time from NEURON (neuron.h._ref_t)
  rec_t = neuron.h.Vector()
  rec_t.record(neuron.h._ref_t)
  # Record Voltage from the center of the soma and the end of dend_2
  rec_v = neuron.h.Vector()
  rec_v.record(soma(0.5)._ref_v)
  rec_2 = neuron.h.Vector()
  rec_2.record(dend_2(1)._ref_v)

  neuron.h.finitialize(-50)
  neuron.init()
  neuron.run(props['run_time'])
  
  if retsoma==True:
    return list(rec_t), list(rec_2), list(rec_v)
  else:
    return list(rec_t), list(rec_2)
