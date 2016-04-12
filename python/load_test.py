from accurateHoc_parallel import *

def load_sample():
  
  from timeit import default_timer as timer

  start_time = timer()
  g = Hoc('home/marderlab/acsutton/code/morphology/sample_data/fake.hoc',
          'home/marderlab/acsutton/code/morphology/sample_data/',
          [1,1,1], 8) # runs on 8 cores
  print('Time taken: %f sec ' %(timer() - start_time))
  return g
