from methods import getData_simple
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def computeLD(snp):
   
  n,m = snp.shape

  LD = np.zeros((m,m))
  Pval = np.ones((m,m))

  for i in range(m):
    for j in range(0,i):

      A0s = (snp[:,i]==0)
      As0 = (snp[:,j]==0)
      A1s = (snp[:,i]==1)
      As1 = (snp[:,j]==1)
      A00 = np.logical_and(A0s, As0)

      P0s = float(np.sum(A0s))/float(n)
      Ps0 = float(np.sum(As0))/float(n)
      P1s = float(np.sum(A1s))/float(n)
      Ps1 = float(np.sum(As1))/float(n)
      P00 = float(np.sum(A00))/float(n)

      D = P00 - P0s*Ps0
      if (D>0):
        denom = min(P0s*Ps1, P1s*Ps0)
        LD[i,j] = D/denom
      else: 
        denom = max(-P0s*Ps0, -P1s*Ps1)
        LD[i,j] = D/denom
  
      Pval[i,j] = 1.0 - stats.chi2.cdf(D*D*n/(P0s*Ps0*P1s*Ps1),1)
        

      #print 'P0s: ' , P0s
      #print 'Ps0: ' , Ps0
      #print 'P00: ' , P00
      #print 'P00 - P0s*Ps0: ' , P00 - P0s*Ps0
  return LD , Pval

if __name__ == "__main__":

  fn = "pop2_mod.txt"
  snp = getData_simple(fn)
  #print snp
  #print snp.shape

  LD , Pval = computeLD(snp)
  #LD , Pval = computeLD(snp[:,0:50])
  #print "LD.shape: ", LD.shape
  #print(LD[2,1])

  mask = np.zeros_like(LD)
  mask[np.triu_indices_from(mask)] = True
  with sns.axes_style("white"):
    plt.figure()
    ax = sns.heatmap(LD, mask=mask, vmax=.5, square=True,  cmap="YlGnBu")
    plt.savefig('LD_heatmap.png')
  with sns.axes_style("white"):
    plt.figure()
    ax = sns.heatmap(Pval, mask=mask, square=True,  cmap="YlGnBu")
    plt.savefig('Pval_heatmap.png')
  with sns.axes_style("white"):
    plt.figure()
    ax = sns.heatmap(-np.log10(Pval), mask=mask, vmin=0, vmax=3., square=True,  cmap="YlGnBu")
    plt.savefig('mlogPval_heatmap.png')
    plt.show()

  #print(Pval)
  #print(-np.log10(Pval))
  #plt.figure()
  #plt.imshow(LD, cmap='hot', interpolation='nearest')
  #plt.show()
  #print LD
