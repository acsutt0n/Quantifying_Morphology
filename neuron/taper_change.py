# pyneuron (python2) test for figures for morphology grant

import neuron
import matplotlib.pyplot as plt
import numpy as np
# most NEURON specific commands start with `neuron.h.`


def set_props(Ra=35, Gleak=0.01, Eleak=-50,finit=-50,
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



def middle_taper(bound0, bound1, props=None, retsoma=False):
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
  
  # Setting passive parameters
  for sec in neuron.h.allsec():
    # Do with the present `sec`
    sec.insert('pas')
    sec.Ra = props['Ra']
    # Do for each segment within `sec`:
    for seg in sec:
      # Do with the segment `seg`:
      seg.pas.g = 0.01
      seg.pas.e = -50

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



def scale_tapers(n=6, props=None):
  # default is from 20->20 to 20->2
  bounds = np.logspace(0, 5, n, base=2)
  #lower_bounds = 80
  lower_bounds = list(reversed(bounds))
  time, volt0, soma = middle_taper(lower_bounds[0], lower_bounds[0], props=props, retsoma=True)
  volts = [volt0]
  for i in lower_bounds[1:]:
    _, v = middle_taper(lower_bounds[0], i, props=props)
    volts.append(v)

  return time, soma, volts, lower_bounds



def plot_tapers(n=6):
  # Plot the stuff
  time, soma, volts, bounds = scale_tapers(n)
  cols = ['b','g','r','c','m','y','k']*4
  
  fig = plt.figure()
  ax = fig.add_subplot(222)
  ax.plot(time, soma, 'k', linewidth=4, alpha=0.4)
  for v in range(len(volts)):
    ax.plot(time, volts[v], cols[v], linewidth=2, alpha=0.7)

  plt.title('Taper vs. Response', fontsize=40)
  plt.xlabel('Time (s)', fontsize=25)
  plt.ylabel('Voltage (mV)', fontsize=25)
  traces = ['Soma - no taper', 'Dend (no taper)']
  for b in range(1, len(bounds)):
    traces.append('Dend (%.1f taper)' %(1-bounds[b]/bounds[0]))
  plt.legend(traces, bbox_to_anchor=(-.5,1), loc=1, borderaxespad=0.)
  plt.yticks([-50.2,-50.],fontsize=25)
  plt.xticks([0,10],fontsize=25)
  plt.show()



def middle_length(midlength=200, props=None, retsoma=False):
  # change the length of the middle neurite
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
  dend_0 = init_section(200, props['Ra'], 20)
  # dendrite1, with taper
  dend_1 = init_section(midlength, props['Ra'], 10)
  diams = list(reversed(np.linspace(10,20, dend_1.nseg)))
  rad = -1
  for seg in dend_1:
    rad = rad+1
    seg.diam = diams[rad]
  # dendrite2
  dend_2 = init_section(200, props['Ra'], 10)
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
  #print('Stimulating at %.1f nA for %.1f s' %(props['stim_amp'],props['stim_dur']))
  
  # Setting passive parameters
  for sec in neuron.h.allsec():
    # Do with the present `sec`
    sec.insert('pas')
    sec.Ra = props['Ra']
    # Do for each segment within `sec`:
    for seg in sec:
      # Do with the segment `seg`:
      seg.pas.g = 0.01
      seg.pas.e = -50

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



def scale_lengths(n=6, props=None):
  # get a response from a neuron of different lengths
  lengths = [200]
  for l in np.linspace(200,2000,n):
    lengths.append(l)
  if len(lengths) > 2:
    lengths.pop(0)
  time, volt0, soma = middle_length(lengths[0], props=props, retsoma=True)
  volts = [volt0]
  for l in lengths[1:]:
    _, v = middle_length(l, props=props)
    volts.append(v)
  
  # for IV curve:
  responses = []
  for v in volts:
    if min(v) < -50:
      responses.append(min(v))
    else:
      responses.append(max(v))
      
  return time, soma, volts, lengths, responses
  


def plot_lengths(n=6, props=None):
  # Plot the stuff
  time, soma, volts, lengths, _ = scale_lengths(n, props)
  cols = ['b','g','r','c','m','y','k']*4
  #print(lengths)
  fig = plt.figure()
  ax = fig.add_subplot(222)
  ax.plot(time, soma, 'k', linewidth=4, alpha=0.4)
  for v in range(len(volts)):
    ax.plot(time, volts[v], cols[v], linewidth=2, alpha=0.7)

  plt.title('Length vs. Response', fontsize=40)
  plt.xlabel('Time (s)', fontsize=25)
  plt.ylabel('Voltage (mV)', fontsize=25)
  traces = ['Soma (%i um)' %(lengths[0]), 'Dend (%i um)' %(lengths[0])]
  for b in range(1,len(lengths)):
    traces.append('Dend (%i um)' %(lengths[b]))
  plt.legend(traces, bbox_to_anchor=(-.5,1), loc=1, borderaxespad=0.)
  plt.yticks([-50.2,-50.],fontsize=25)
  plt.xticks([0,10],fontsize=25)
  plt.show()



def length_IV_curve(Isteps=[-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3], 
                    n=6):
  responses = []
  Isteps = [i*3 for i in Isteps]
  for i in Isteps:
    props = set_props(stim_amp=i)
    _, _, _, lengths, reps = scale_lengths(n, props=props)
    responses.append(reps)
  
  fig = plt.figure()
  ax = fig.add_subplot(222)
  for i in range(n):
    reps = [r[i]+50. for r in responses]
    ax.plot(Isteps, reps, 'k',linewidth=2, alpha=0.5)

  traces = ['%i um' %(l) for l in lengths]
  # plt.legend(traces, loc='best')
  ax.spines['left'].set_position('zero')
  ax.spines['bottom'].set_position('zero')
  plt.xticks([-10,10], fontsize=25)
  plt.yticks([-1,1], fontsize=25)
  plt.xlabel('Applied current (nA)', fontsize=25)
  plt.ylabel('Response (mV)', fontsize=25)
  ax.yaxis.set_label_coords(0,0.5)
  ax.xaxis.set_label_coords(0.5,0)
  plt.title('I-V curve for passive neuron', fontsize=40)
  plt.show()
  
  return


def IV(i=10):
  props=set_props(stim_amp=i)
  #responses = []
  time, soma, volts, lengths, response = scale_lengths(n=1, props=props)
  print(response)
  fig = plt.figure()
  plt.plot(time, soma, 'k')
  plt.xlabel('Applied current (nA)', fontsize=25)
  plt.ylabel('Response (mV)', fontsize=25)
  plt.title('I-V curve for passive neuron', fontsize=30)
  traces = ['%i um' %(l) for l in lengths]
  #plt.legend(traces, loc='best')
  plt.show()
  return


############
if __name__ == '__main__':
  #IV()
  #length_IV_curve()
  #plot_lengths()
  plot_tapers()
  




"""
# ############### ###### ######### ############### Old dev stuff
# Plot the recordings with matplotlib
# ===================================

import matplotlib.pyplot as plt

# get values from NEURON-vector format into Python format
times = [] # Use list to add another trace later.
voltages = []
times.append(list(rec_t)) # alternativ to `list(rec_t)`: `numpy.array(rec_t)`
voltages.append(list(rec_v))
voltages.append(list(rec_2))
# check types by:
# >>> type(rec_t)
# >>> type(time[0])

fig = plt.figure()
plt.plot(times[0], voltages[0], 'b')
plt.plot(times[0], voltages[1], 'r')
plt.title("sample plot")
plt.xlabel("time (ms)")
plt.ylabel("voltage (mV)")

plt.show() # If the interpreter stops now: close the figure.
# For interactive plotting, see `Part 1` -> `ipython`

# Does the voltage trace look as you expected?
"""




