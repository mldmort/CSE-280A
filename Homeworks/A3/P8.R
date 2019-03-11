file_name = '/Users/miladmortazavi/Documents/Dev/BIOINF/CSE-280A/Homeworks/A3/a3data2.txt'

data = as.matrix(read.table(text=gsub("",' ',readLines(file_name))))
#n = dim(data)[1]
n = 20
m = dim(data)[2]

q = array(0,dim=c(n,n,m))

#q = numeric(n*n*m)

for (i in 1:n) {
  for (j in 1:n) {
    for (k in 1:m) {
      if ((data[i,k]==1) && (data[j,k]==1)) {
        q[i,j,k] = 1
      }
    }
  }
}

f = 0.1

C = matrix(0,nrow = m+n,ncol = 1)
for (i in 1:m) {
  C[i] = 1
}

A = matrix(0,nrow = (n*(n-1)/2*m)+1,ncol = m+n)
rhs = matrix(0,nrow = (n*(n-1)/2*m)+1,ncol = 1)

row = 1
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    #if (i) {
      for (k in 1:m){
        A[row,(i+m)] = 1
        A[row,(j+m)] = 1
        A[row,k] = 1
        rhs[row] = 2+q[i,j,k]
        row = row + 1
      }
    #}
  }
}
for (i in (m+1):(m+n)) {
  A[row,i] = 1
}
rhs[row] = f*n

sense <- rep(c('<='),((n*(n-1)/2*m)+1))
sense[n*(n-1)/2*m+1] = '>='

library(slam)
library(gurobi)

model <- list()

model$A          <- A
model$obj        <- C
model$modelsense <- 'max'
model$rhs        <- rhs
model$sense      <- sense
model$vtype      <- 'B'

params <- list(OutputFlag=0)

result <- gurobi(model, params)

print('Solution:')
print(result$objval)
print(result$x)
#print(sum(result$x)/length(result$x)*100.)