import numpy as np
from methods import *

if __name__ == "__main__":

  fn = 'a1globedata.txt'
  snp = getData_table(fn)
  print snp

  print

  #make_zero_major(snp)
  #print snp

  #print 

  sort_snp(snp)
  print snp
  #ans = is_pp(snp)
  #print 'is pp: ' , ans
  #print snp
  #print snp.shape


