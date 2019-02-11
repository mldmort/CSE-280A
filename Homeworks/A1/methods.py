import numpy as np

def getData_simple(fn):
  fl = open(fn,'r')
  lines = fl.readlines()
  n = len(lines)
  m = len(lines[0]) - 1
  M = np.zeros((n,m),int)
  for i in range(n):
    for j in range(m):
      M[i][j] = int(lines[i][j])
  return M

  fl.close()

# be careful: this method adds a first row to track the location of each column when sorted
# in method make_zero_major we hacked the code so not to consider the first row
# when sorting it goes throught the first row but doesn't change the result
def getData_table(fn):
  M = np.loadtxt(fn, skiprows = 1, dtype = int)
  col_ids = range(1,M.shape[1])
  return np.vstack([col_ids,M[:,1:]])
  
def swap_columns(M,c1,c2):
  col = list(M[:,c2])
  M[:,c2] = list(M[:,c1])
  M[:,c1] = col

def sort_row(M,c1,c2,row):
  #print ' > sorting columns: c1: ', c1, ' c2: ' , c2, ' row: ' , row
  if ((c1>=c2) or (row>=M.shape[0])):
    return
  swap_loc = c1
  for ci in range(c1,c2+1):
    if (M[row][ci] == 0):
      break
    swap_loc += 1
  #print 'swap_loc: ' , swap_loc
  for ci in range(swap_loc+1,c2+1):
    if (M[row][ci] == 1):
      swap_columns(M,ci,swap_loc)
      swap_loc += 1

  #print M
  if (c1 < swap_loc-1):
    sort_row(M,c1,swap_loc-1,row+1)
  if (swap_loc < c2):
    sort_row(M,swap_loc,c2,row+1)

def sort_snp(M):
  row, col = M.shape
  sort_row(M,0,col-1,0)
  
# checks perfect phylogeny
def is_pp(M):
  sort_snp(M)
  row, col = M.shape
  for ic in range(col-1):
    if ( not is_pp_columns(M,ic,ic+1)):
      return False

  return True

def is_pp_columns(M,c1,c2):
  row, col = M.shape
  col1 = M[:,c1]
  col2 = M[:,c2]
  col3 = col1 * col2
  if( np.array_equal(col3,col1) or np.array_equal(col3,col2) or np.array_equal(col3,np.zeros_like(col3)) ):
    return True
  else:
    return False

def make_zero_major(M):
  row, col = M.shape
  for ic in range(col):
    if (sum(M[1:,ic]) >= row/2):
      M[1:,ic] = 1 - M[1:,ic]

