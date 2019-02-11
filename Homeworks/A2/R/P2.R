n = 10
theta = 40
#N = 1e6
N = 1000
#mu = theta/(4*N)
mu = 1e-1
l = 1e6

tree = list(list())

for (i in 1:n) {
  tree[[1]][i] = list(i)
}

# making a random tree
for (k in n:2) {
  tn = n+1-k
  rand2 = sample.int(k,2)
  tree[[tn+1]] = list()

  merge_set = unlist(union(tree[[tn]][rand2[1]], tree[[tn]][rand2[2]]))
  rem_set = tree[[tn]]
  rem_set = rem_set[-c(rand2[1],rand2[2])]

  tree[[tn+1]][[1]] = (merge_set)
  tree[[tn+1]] = c(tree[[tn+1]],(rem_set))

  #print('tn+1: ')
  #print(tn+1)
  #print(merge_set)
  #print(rem_set)
  #print(tree[[tn+1]])
}

# find the mutations
snp = matrix(0,nrow=n,ncol=1)

#n = 10
#N = 1000
#mu = 1e-1
mut_nums = c()
for (k in n:2) {
  k_ch_2 = k*(k-1)/2.0
  Tk = N/(k_ch_2)
  #print("Tk:")
  #print(Tk)
  p_mean = Tk*mu
  #print('p_mean:')
  #print(p_mean)
  rand_vec = rpois(k,p_mean)
  #print('rand_vec:')
  #print(rand_vec)
  mut_nums = c(mut_nums,rand_vec)
}

mut_sum = sum(mut_nums)
mut_locs = sample.int(l,mut_sum)

snp = matrix(0,nrow=n,ncol=mut_sum)
icol = 1
ibranch = 1
for (k in n:2) {
  tn = n+1-k
  for (i in 1:k) {
    if (mut_nums[ibranch] >= 1) {
      for (imu in 1:mut_nums[ibranch]) {
        snp[c(unlist(tree[[tn]][i])),icol] = 1
        icol = icol+1
      }
    }
    ibranch = ibranch+1
  }
  print('icol:')
  print(icol)
}

colnames(snp) = mut_locs




