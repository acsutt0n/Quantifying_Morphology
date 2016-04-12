#!/usr/bin/python


from multiprocessing import Process, Queue, cpu_count
import tempfile
import os
from operator import xor
import pickle
import textProgress
import sys
if sys.version_info[0] == 3:
  from io import StringIO
else:
  from cStringIO import StringIO
from time import sleep
import functools


###############################################################################
def robust_map(f, inputList, args=tuple(), kwargs=dict(), numProcesses=-1,
               initFunc=None):
  rHash = hash(f.__name__)
  rHash ^= functools.reduce(xor, map(hash, inputList))
  resumeFile = os.path.join(tempfile.gettempdir(), 'resume_' + str(rHash))
  logFile = 'map_' + str(rHash) + '.log'
  logOut = None
  inputQueue = Queue()
  outputQueue = Queue()

  if os.access(resumeFile, os.R_OK):
    print('Resuming previously interrupted map (%s)' % resumeFile)
    # resume interrupted job
    computedIndices = []
    with open(resumeFile, 'rb') as fIn:
      fIn.seek(0, 2)
      endInd = fIn.tell()
      fIn.seek(0)
      while fIn.tell() < endInd:
        # load pickled calculated obj (index, value, textOut) from resume file
        computedVal = pickle.load(open('fIn','rb'))
        # add computed index to computed indices
        computedIndices.append(computedVal[0])
        # add computed obj to outputQueue
        outputQueue.put(computedVal + (None,))
        textOut = computedVal[-1]
        if textOut:
          if logOut is None:
            print('Logging output to %s' % logFile)
            logOut = open(logFile, 'w')
          logOut.write(textOut)
    for ind, inputVal in enumerate(inputList):
      if ind not in computedIndices:
        inputQueue.put((ind, inputVal))
    print('Computing new values')
  else:
    # start new job
    for inputVal in enumerate(inputList):
      inputQueue.put(inputVal)
  for n in range(numProcesses):
    inputQueue.put((None, None))
  
  if numProcesses <= 0:
    numProcesses += cpu_count()
  
  sendArgs = (f, inputQueue, outputQueue, initFunc) + args
  procs = [Process(target=_mapFunc, args=sendArgs, kwargs=kwargs)
           for n in range(numProcesses)]

  numInput = len(inputList)  
  outputList = [None] * numInput
  textOutList = [''] * numInput
  
  # start the processes
  for p in procs:
    p.start()

  numOutput = 0
  numPreCalc = outputQueue.qsize()
  textProgress.startProgress(numInput - numPreCalc)
  while numOutput < numInput:
    try:
      ind, val, textOut, err = outputQueue.get()
    except BaseException as err:
      pass
    
    if textOut:
      if logOut is None:
        print('Logging output to %s' % logFile)
        logOut = open(logFile, 'w')
      logOut.write(textOut)
    
    if err is None:
      numOutput += 1
      if numOutput > numPreCalc:
        textProgress.updateProgress()
      outputList[ind] = val
      textOutList[ind] = textOut
    else:
      # there was a problem, save current progress
      with open(resumeFile, 'wb') as fOut:
        for ind, (val, textOut) in enumerate(zip(outputList, textOutList)):
          if val is not None:
            pickle.dump((ind, val, textOut), open(fOut, 'rb'))
      # stop all the functioning Processes
      for p in procs:
        if p.is_alive():
          p.terminate()
      # clean up log
      if logOut is not None:
        logOut.close()
      # re-raise the error
      raise err

  # wait until all procs finished
  for p in procs:
    p.join()

  assert outputQueue.empty()
  
  if logOut is not None:
    logOut.close()
  
  if os.access(resumeFile, os.F_OK):
    # resume file exists, but map concluded successfully
    os.remove(resumeFile)
  
  return outputList


def _mapFunc(f, inputQueue, outputQueue, initFunc, *args, **kwargs):
  import sys
  import traceback
  sys.stdout = textOutput = StringIO()
  
  if initFunc is not None:
    initFunc()
  try:
    while not inputQueue.empty():
      ind, inputVal = inputQueue.get()
      if inputVal is None:
        break
      textOutput.seek(0)
      outputVal = f(inputVal, *args, **kwargs)
      outputQueue.put((ind, outputVal, textOutput.getvalue(), None))
  except BaseException as err:
    if err.args and err.args[0]:
      if not isinstance(err.args[0], str):
        err.args = (str(err.args[0]) + '\nProcess traceback:\n' +
                  traceback.format_exc(),) + err.args[1:]
      else:
        err.args = (err.args[0] + '\nProcess traceback:\n' +
                  traceback.format_exc(),) + err.args[1:]
    else:
      err.args = ('\nProcess traceback:\n' + traceback.format_exc(),) + \
                 err.args[1:]
    try:
      outputQueue.put((ind, None, textOutput.getvalue(), err))
    except:
      pass
    raise err


def test_robust_map():
  
  from time import sleep
  from scipy import linspace
  import random
  
  def _testFunc(x, myStr='stupid'):
    sleep(1)
    r = random.uniform(0.0, 100.0)
    print(myStr)
    ans = x * x
    if r < 2.0:
      raise RuntimeError('Oh Snap! %g' % r)
    return ans, myStr

  inputList = linspace(0, 10, 100)
  output = robust_map(_testFunc, inputList, numProcesses=-1,
                      kwargs={'myStr':'called _testFunc()'},
                      initFunc=random.seed)
  vals, strings = zip(*output)
  print(vals)


if __name__ == "__main__":
  test_robust_map()
