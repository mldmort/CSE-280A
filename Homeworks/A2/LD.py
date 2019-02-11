from methods import getData_simple
import numpy as np

def computeLD(snp):
   
  n,m = snp.shape

  LD = np.zeros_like(snp)

  for i in xrange(n):
    for j in xrange(i+1,n):

      A0s = (snp[:,i]==0)
      As0 = (snp[:,j]==0)
      A00 = np.logical_and(A0s, As0)

      P0s = np.sum(A0s)/float(n)
      Ps0 = np.sum(As0)/float(n)
      P00 = np.sum(A00)/float(n)

      LD[i,j] = P00 - P0s*Ps0

      #print 'P0s: ' , P0s
      #print 'Ps0: ' , Ps0
      #print 'P00: ' , P00
      #print 'P00 - P0s*Ps0: ' , P00 - P0s*Ps0
  return LD

if __name__ == "__main__":

  fn = "pop2_mod.txt"
  snp = getData_simple(fn)
  print snp
  print snp.shape

  LD = computeLD(snp)
  print LD
