import numpy as np
import time
import sys
from methods import *

if __name__ == "__main__":

  fn = 'a1data1.txt'
  #fn = 'a1data2.txt'
  #fn = 'a1data3.txt'
  #fn = 'a1data4.txt'
  #fn = 'a1data6.txt'

  snp = getData_simple(fn)

  if ((fn == 'a1data3.txt') or (fn == 'a1data4.txt')):
    # for 'a1data3.txt' and 'a1data4.txt'
    sys.setrecursionlimit(2000)
  if (fn == 'a1data6.txt'):
    # for 'a1data6.txt'
    sys.setrecursionlimit(20000)

  #print snp
  #swap_columns(snp,0,1)
  #print

  start = time.time()

  #sort_snp(snp)
  ans = is_pp(snp)
  print 'is perfect phylogeny: ' , ans

  end = time.time()

  print 'time: ', end - start

  #print snp

  outfn = 'output_' + fn
  np.savetxt(outfn,snp,fmt='%d',delimiter='')
