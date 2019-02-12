
n = 100             # Sample size
N0 = 1000000       # Population size
mu = 1e-8          # number of mutations per generation per bp
theta = 0.04       # Scaled mutation rate
l = 10000          # number of base pairs
alph = 7/(4*N0)
#alph = 0

nloops = 200
xsi = matrix(0, ncol = (n-1))

for (il in 1:nloops) {
  gen = 0
  tree = list(list())

  for (i in 1:n) {
    tree[[1]][i] = list(i)
  }

  # Making a random tree
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

  # Find the mutations
  snp = matrix(0,nrow=n,ncol=1)
  mut_nums = c()

  for (k in n:2) {
    k_ch_2 = k*(k-1)/2
    Nt = exp(-alph*gen)*N0
    Tk = rgeom(1, k_ch_2/(2*Nt))
    #Tk = N0/(k_ch_2)
    gen = gen + Tk
  #  Tk = N/(k_ch_2)
    #print("Tk:")
    #print(Tk)
    p_mean = l*Tk*mu
    #print('p_mean:')
    #print(p_mean)
    rand_vec = rpois(k,p_mean)
    #print('rand_vec:')
    #print(rand_vec)
    mut_nums = c(mut_nums,rand_vec)
  }
  
  # mut_sum is the total number of mutations
  mut_sum = sum(mut_nums)
  mut_locs = sample.int(l,mut_sum)
  
  snp = matrix(0, nrow = n, ncol = mut_sum)
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
      ibranch = ibranch + 1
    }
  #  print('icol:')
  #  print(icol)
  }
  
  colnames(snp) = mut_locs
  
  Num = matrix(0, ncol = dim(snp)[2])
  for (i in 1:dim(snp)[2]) {
    Num[i] = sum(snp[,i])
  }
  
  for (i in 1:(n-1)) {
    xsi[i] = xsi[i] + length(which(Num == i))
  }
  
}

xsi = xsi / nloops
xsi = xsi * c(1:(n-1)) / l

png("P2_alpha=7.png")
plot(1:(n-1), xsi, type="b", xlab="i", ylab="i*xsi", main='alpha=7/(4*N)')
#lines(c(1,(n-1)),c(0.04,0.04), type="l")
dev.off()

