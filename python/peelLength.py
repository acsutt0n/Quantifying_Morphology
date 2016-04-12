#!/usr/bin/python


import sys
write = sys.stdout.write
import scipy
from scipy import optimize, diff
from math import exp, log, isinf, sqrt
from matplotlib import pyplot


###############################################################################
def expSum(t, model, vErr=0, t0=None):
  if t0 is None:
    v = [vErr + sum(dV * (1.0 - exp(-tn/tau)) for tau,dV in model) for tn in t]
  else:
    def _f(_dV, _tn, _t0, _tau):
      if _tn <= _t0:
        return 0.0
      else:
        return _dV * (1.0 - exp((_t0 - _tn)/_tau))
    v = [vErr + sum(_f(dV, tn, t0, tau) for tau,dV in model) for tn in t]
  
  return v


###############################################################################
def expSumParams(t, *params):
  offset = params[0]
  if any(p <= 0 for p in params[1::2]):
    return scipy.array([float('inf')] * len(t))
  try:
    v = scipy.array([offset + sum(dV * (1.0 - exp(-tn/tau)) for tau, dV in
                                  zip(params[1::2], params[2::2]))
                     for tn in t])
  except OverflowError:
    print(params)
    raise
  return v


###############################################################################
def _linearFit(x, y, startInd, stopInd):
  x = x[startInd:stopInd]
  y = y[startInd:stopInd]
  p, residuals, rank, singular_values, rcond = scipy.polyfit(x, y, 1,full=True)
  linearErr = sum(e*e for e in residuals)
  return p[0], p[1], linearErr, startInd, stopInd


###############################################################################
def estimateVInf(v, vNoise=1e-10):
  stopInd, maxV = max(enumerate(v), key=lambda x: x[1])
  
  im = int(stopInd * 0.9)
  i1 = int(stopInd * 0.95)
  ip = stopInd
  
  vm, v1, vp = v[im], v[i1], v[ip]

  dVm = v1 - vm
  dVp = vp - v1
  if dVp < vNoise or dVm < vNoise or dVm < dVp:
    vInf = vp
  else:
    vInf = v1 + dVm * dVp / ( dVm - dVp )
  return vInf, stopInd


###############################################################################
def peelLength(t, v, startInd=0, debugPlots=False):
  # Model a voltage deflection as v(t) = sum_n dV_n * (1.0 - exp(-t/tau_n))
  # Find the model by looking at the late time period to find the slowest
  #  exponential, then subtracting that exponential to find the next slowest,
  #  etc
  
  # convert v to a scipy array
  v = scipy.array(v)
  
  # get tStart, vStart
  tStart, vStart = t[startInd], v[startInd]
  fitStartInd = startInd
  
  # estimate the final voltage after infinite time
  vInf, stopInd = estimateVInf(v)
  
  # convert v to exponential decay towards vInf
  v = vInf - v
  
  # compute and save unexplained voltage
  vStart = v[startInd]
  
  # start with empty model
  model = []

  # loop until there's essentially no deltaV remaining
  while v[startInd] > 0.01 * vStart:
    # find the slowest remaining time constant, and it's associated dV
    
    # fix the duration of the search window
    dI = max(10, (stopInd - startInd) / 10)
    if dI >= stopInd - startInd:
      # not enough room to search, so give up
      break
    
    oldStopInd = stopInd
    oldFitStartInd = fitStartInd
    # find the best subregion/linear fit to log v in that subregion
    logV = [log(vn) for vn in v[:stopInd]]
    slope, offset, linearErr, fitStartInd, stopInd = \
      min((_linearFit(t, logV, i1, i1 + dI)
           for i1 in range(startInd, stopInd - dI)),
          key=lambda x: x[2])
    
    tau = -1.0 / slope
    dV = exp(offset)
    
    if debugPlots:
      # compute the change in log of the slowest part of the voltage decay
      linear = [offset + slope * tn for tn in t[:oldStopInd]]
      fig = pyplot.figure()
      ax = fig.add_subplot(111)
      ax.plot(t[:oldStopInd] - t[0], logV)
      ax.plot(t[:oldStopInd] - t[0], linear, 'r--', linewidth=2)
      ax.plot(t[fitStartInd:stopInd] - t[0], linear[fitStartInd:stopInd], 'g-',
              linewidth=2)
      pyplot.xlim(0.0, t[oldStopInd] - t[0])
      pyplot.ylim(min(min(linear[:oldStopInd]), min(logV[:oldStopInd])),
                  max(linear[0], logV[0]))

    if tau < 0 or (len(model) > 0 and tau > model[-1][0]):
      # the exponential we've found is garbage
      break
    # the exponential is okay, add it to the model
    model.append((tau, dV))
    
    # remove the newly modeled exponential from the voltage
    def _expFun(_t, _dV, _slope):
      if _t > 0.0:
        return _dV * exp(_slope * _t)
      else:
        return _dV
    v -= [_expFun(tn, dV, slope) for tn in t]

    # refine the stopInd for next time
    try:
      stopInd = min(stopInd,
                    next(vn for vn in enumerate(v) if vn[1] <= 0)[0] - 1)
    except StopIteration:
      pass
  
  vErr = v[startInd]
  numV = len(v) - startInd
  vResid = sqrt(sum(vn*vn for vn in v[startInd:]) / numV)
  return model, vErr, vResid, oldFitStartInd


###############################################################################
def fitLength(t, v, startInd, startModel, vErr=0):
  
  expInd = max(enumerate(diff(v)), key=lambda x:x[1])[0] + 1
  startInd = max(expInd, startInd)
  
  # only fit the exponential part of the traces
  t, v = t[startInd:], v[startInd:]
  
  # start with initial guess and fit parameters of model
  startParams = [vErr] + [p for pair in startModel for p in pair]
  
  try:
    params, pCov = optimize.curve_fit(expSumParams, t, v, p0=startParams,
                                      maxfev=500)
  except RuntimeError as err:
    if 'Number of calls to function has reached maxfev' in err.message:
      print(err.message)
      return [], float('Inf'), float('inf')
    else:
      raise
    
  fitModel = [(tau, dV) for tau, dV in zip(params[1::2], params[2::2])]
  vErr = params[0]
  fitV = expSum(t, fitModel, vErr=vErr)
  vResid = sqrt(sum((vn - fitVn)**2 for vn, fitVn in zip(v, fitV)) / len(fitV))
  return fitModel, vErr, vResid


###############################################################################
def visualizeFit(t, v, fitModel, vErr=0):
  fitV = expSum(t, fitModel, vErr=vErr)
  fig = pyplot.figure()
  ax = fig.add_subplot(111)
  ax.plot(t, v, 'b-')
  ax.plot(t, fitV, 'r--')
  pyplot.ylim(0, max(max(v), max(fitV)))


###############################################################################
def printModel(model, fit=True, vErr=None, vResid=None):
  numExp = len(model)
  if fit:
    preamble = 'Found %d exponentials:   ' % numExp
  else:
    preamble = 'Created %d exponentials: ' % numExp
  if numExp == 0:
    print(preamble)
    return
  
  print(preamble + 'tau = %6.2f ms, dV = %5.2f mV' % model[0])
  indent = ' ' * 24
  for expVals in model[1:]:
    print(indent + 'tau = %6.2f ms, dV = %5.2f mV' % expVals)
  if vErr is not None:
    print(indent + 'Remaining unexplained delta v = %5.2f mV' % vErr)
  if vResid is not None:
    print(indent + 'Residual rms voltage error = %5.2f mV' % vResid)


###############################################################################
def getStepWindow(t, v):
  # return time and voltage vectors during the stimulus period only
  
  # find the point of maximum voltage, and cut off everything afterwards
  maxInd, maxV = max(enumerate(v), key=lambda x: x[1])
  minInd, minV = min(enumerate(v), key=lambda x: x[1])
  if maxV - v[0] > v[0] - minV:
    # this is a positive step
    t = t[:maxInd]
    v = scipy.array(v[:maxInd])
  else:
    # this is a negative step, flip it for now
    t = t[:minInd]
    v = v[0] - scipy.array(v[:minInd])
  
  # re-center time to start at the point of maximum voltage change
  diffV = diff(v)
  dVInd, maxDV = max(enumerate(diffV), key=lambda x: x[1])
  dVInd -= 1
  while diffV[dVInd] > 0:
    dVInd -= 1
  dVInd += 1
  
  t -= t[dVInd]
  v -= v[dVInd]
  
  return t, v, dVInd
  

###############################################################################
def modelResponse(t, v, verbose=False, plotFit=False, debugPlots=False,
                  findStepWindow=False, displayModel=None):
  # Model voltage response as sum of decaying exponentials
  if findStepWindow:
    t, v, startInd = getStepWindow(t, v)
  else:
    startInd = 0
  
  if displayModel is None:
    displayModel = verbose
  
  if verbose:
    print('Modeling voltage by peeling exponentials')
  fitModel1, vErr1, vResid1, fitStartInd \
    = peelLength(t, v, startInd, debugPlots=debugPlots)
  if displayModel:
    printModel(fitModel1, vErr=vErr1, vResid=vResid1)
  if verbose:
    print('Performing nonlinear fit to voltage data')
  fitModel2, vErr2, vResid2 = fitLength(t, v, fitStartInd, fitModel1, vErr=vErr1)
  if displayModel:
    printModel(fitModel2, vErr=vErr2, vResid=vResid1)
  
  if vResid1 < vResid2:
    if verbose:
      print('Peeling exponentials gives best result')
    fitModel, vErr, vResid = fitModel1, vErr1, vResid1
  else:
    if verbose:
      print('Nonlinear fit gives best result')
    fitModel, vErr, vResid = fitModel2, vErr2, vResid2
  
  if plotFit:
    visualizeFit(t, v, fitModel, vErr=vErr)

  return fitModel, vErr, vResid


###############################################################################
def evaluateFit(model, fitModel, vErr, vResid):
  def _percErr(modelExp, fitExp):
    return (100 * (fitExp[0] / modelExp[0] - 1),
            100 * (fitExp[1] / modelExp[1] - 1))
  print('Steady state voltage error = %.2g mV' % vErr)
  print('RMS voltage error = %.2g mv' % vResid)
  print('Percent model errors: ' + 
        'tau = %4.2g%%, dV = %4.2g%%' % _percErr(model[0], fitModel[0]))
  indent = ' ' * 22
  for modelExp, fitExp in zip(model[1:], fitModel[1:]):
    print(indent + 'tau = %4.2g%%, dV = %4.2g%%' % _percErr(modelExp, fitExp))


###############################################################################
def getDemoTrace(dataDir):
  import os
  tList = [f for f in os.listdir(dataDir) if f.startswith('time')]
  
  import cPickle
  for tFile in tList:
    vFile = tFile.replace('time', 'voltage')
    
    with open(os.path.join(dataDir, tFile), 'r') as fIn:
      t = cPickle.loads(fIn.read())
    with open(os.path.join(dataDir, vFile), 'r') as fIn:
      v = cPickle.loads(fIn.read())
    yield t, v


###############################################################################
def demo(plotFit=False, debugPlots=False, verbose=True):
  #demoDir = 'data/first broken'
  #demoDir = 'data/saverangeparameters'
  #demoDir = 'data/High Residual Error'
  demoDir = 'data/High Residual Error 2'
  for t, v in getDemoTrace(demoDir):
    model = modelResponse(t, v, verbose=verbose, plotFit=plotFit,
                          debugPlots=debugPlots)
  return model


###############################################################################
def makeSyntheticData(deltaVTotal=10.0, tFinal=2000, dT=0.25, numLengths=(2,5),
                      noiseAmp=0.0, upSwing=False):
  import random
  random.seed()
  numLengths = random.randint(*numLengths)
  tau = tFinal * 0.5 * \
    scipy.cumprod([random.uniform(0.1, 0.5) for n in range(numLengths)])
  dV = scipy.array([random.uniform(0.01, 1)**3 for n in range(numLengths)])
  dV *= (deltaVTotal / sum(dV))  
  model = [(tau_n, dV_n) for tau_n, dV_n in zip(tau, dV)]
  
  t = scipy.linspace(0, tFinal, int(tFinal/dT) + 1)
  v = expSum(t, model)
  
  if upSwing:
    tauUp = min(0.5 * min(tau), 10.0)  
    for n, tn in enumerate(t):
      v[n] *= (1.0 - exp(-tn / tauUp))
    
  if noiseAmp > 0:
    v = [vn + random.normalvariate(0, noiseAmp) for vn in v]
  return t, v, model


###############################################################################
def testPeelLength(plotFit=False, debugPlots=False, noiseAmp=0.0,
                   verbose=True):
  if verbose:
    print('Creating synthetic data')
  t, v, model = makeSyntheticData(noiseAmp=noiseAmp)
  if verbose:
    printModel(model, fit=False)
  fitModel, vErr, vResid = modelResponse(t, v, plotFit=plotFit,
                                         debugPlots=debugPlots,
                                         verbose=verbose)
  evaluateFit(model, fitModel, vErr, vResid)


###############################################################################
def _parseArguments():
  import argparse
  parser = argparse.ArgumentParser(description=
    "Simulate a neuron with geometry exported in a .hoc file, and passive "
    + "properties specified in a separate json .txt file",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("--demo", action='store_true',
                      help="test previously-generated data in demo dir, "
                          +"instead of synthetically generating new data")
  parser.add_argument("--plotFit", "-f", action='store_true',
                      help="visualize the original data and fit model")
  parser.add_argument("--debugPlots", "-d", action="store_true",
                      help="make additional plots demonstrating fit process")
  parser.add_argument("--quiet", "-q", action="store_false", dest="verbose",
                      help="suppress printing information")
  parser.add_argument("--noiseAmp", "-n", default=0.0, type=float,
                      help="specify noise level for synthetic data")
              
  return parser.parse_args()


###############################################################################
if __name__ == "__main__":
  options = _parseArguments()
  if options.demo:
    demo(plotFit=options.plotFit, debugPlots=options.debugPlots,
         verbose=options.verbose)
  else:
    testPeelLength(plotFit=options.plotFit, debugPlots=options.debugPlots,
                   noiseAmp=options.noiseAmp, verbose=options.verbose)
  
  if options.plotFit or options.debugPlots:
    pyplot.show()
