#!/usr/bin/python


import sys
import time

_totalTasks = None
_finishedTasks = None
_startTime = None
_message = 'Completed'

def startProgress(totalTasks, message=None):
  global _totalTasks
  global _finishedTasks
  global _startTime
  global _message
  
  _totalTasks = totalTasks
  _finishedTasks = 0
  _startTime = time.time()
  if message is not None:
    _message = message


def updateProgress(tasks=1):
  global _finishedTasks
  _finishedTasks += tasks

  elapsedTime = time.time() - _startTime

  if _finishedTasks == _totalTasks:
    elapsedStr = getTimeString(elapsedTime)
    
    sys.stdout.write('\r%s %d/%d = %.1f%%. %s elapsed time.\n' % \
      (_message, _finishedTasks, _totalTasks, 100.0, elapsedStr))
  else:
    ratio = _totalTasks / float(_finishedTasks)
    remainTime = (ratio - 1.0) * elapsedTime
    remainStr = getTimeString(remainTime)
    
    sys.stdout.write('\r%s %d/%d = %.1f%%. %s remaining.' % \
      (_message, _finishedTasks, _totalTasks, 100 / ratio, remainStr))
  
  sys.stdout.flush()


def getTimeString(timeInterval):
  # return a human-friendly string describing a time interval in seconds
  m = int(timeInterval / 60.0)
  if m == 0:
    return '%.3fs' % timeInterval
  else:
    s = int(round(timeInterval % 60.0))
    h = m / 60
    if h == 0:
      return '%dm %ds' % (m, s)
    else:
      m = m % 60
      d = h / 24
      if d == 0:
        return '%dh %dm %ds' % (h, m, s)
      else:
        h = h % 24
        return '%dd %dh %dm %ds' % (d, h, m, s)


def testProgress(numTasks = 101, totalTime = 15.0, message='Testing progress'):
  delayTime = totalTime / numTasks
  startProgress(numTasks, message)
  for n in range(numTasks):
    time.sleep(delayTime)
    updateProgress()


if __name__ == "__main__":
  testProgress()
  
  
