# python test return



def mainfun():
  
  def submainfun(submainin):
    
    if submainin == 1:
      print('submain input = 1')
    else:
      print('submain input != 1')
      
  secondfuninput1 = 1
  name1 = secondfun(secondfuninput1)
  submainfun(name1)
  
  secondfuninput2 = 0
  name2 = secondfun(secondfuninput2)
  submainfun(name2)


def secondfun(secondfunin):
  
  if secondfunin == 1:
    print('secondfun input = 1')
    return 0
  elif secondfunin == 0:
    print('secondfun input = 0')
    return [
  else:
    print('bad input to secondfun')

if __name__ == '__main__':
  mainfun()
  
