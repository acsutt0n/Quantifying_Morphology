# fake neuron from p. 131 of NEURON book
# uses python interpreter, NOT hoc 


from neuron import *
from neuron import h
from nrn import *


############### topology ##############
"""
Could use h object, as below:
h('create soma, apical, basilar, axon')
h('connect apical(0), soma(1)')
h('connect basilar(0), soma(0)')
h('connect axon(0), soma(0)')
"""
soma = Section() # section class from nrn
apical = Section()
basilar = Section()
axon = Section()
apical.connect(soma) # child.connect(parent)
basilar.connect(soma, 0, 0)
axon.connect(soma, 0, 0)
sections = [soma, apical, basilar, axon]


############### geometry ##################
"""
Again could use lots of h objects but that's obnoxious.
"""
soma.L = 30
soma.diam = 30
soma.nseg = 1

apical.L = 600
apical.diam = 1
apical.nseg = 23

basilar.L = 200
basilar.diam = 2
basilar.nseg = 5

axon.L = 1000
axon.diam = 1
axon.nseg = 37


############### biophysics ###############

for sec in sections:
  sec.Ra = 100
  sec.cm = 1

soma.insert('hh')
soma.insert('pas')
soma().g_pas = 0.0002
soma().e_pas = -65

basilar.insert('pas')
basilar.g_pas = 0.0002
basilar.e_pas = -65

axon.insert('hh')


################# testing ###################

# don't know how to avoid using the hoc interpreter for this part...
h('topology()') # prints the crude version of the topology
h('forall psection()') # prints geometry and biophysical properties


################# instrumentation ##################

# synaptic input
syn = h.AlphaSynapse(0.5, sec=soma)
syn.onset = 0.5
syn.tau = 0.1
syn.gmax = 0.05
syn.e = 0

# graphical display
g = h.Graph()
g.size(0,5,-80,40)

















