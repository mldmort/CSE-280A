file_name = '/Users/miladmortazavi/Documents/Dev/BIOINF/CSE-280A/Homeworks/A3/recombination_simulation/snp.r.40_mod.dat'

data = as.matrix(read.table(text=gsub("",' ',readLines(file_name))))
n = dim(data)[1]
m = dim(data)[2]

q = matrix(1, nrow=m, ncol=m)

for (i in 1:m) {
  for (j in i:m) {
    temp = data[,i]*data[,j]
    if ((sum(temp)==0) || all(temp==data[,i]) || all(temp==data[,j])) {
      q[i,j] = 1
      q[j,i] = 1
    } else {
      q[i,j] = -20000
      q[j,i] = -20000
    }
  }
}

#q[,5] = -q[,5]
#q[5,] = -q[5,]

mult = matrix(0,nrow=m,ncol=1)
const = matrix(0,nrow=m,ncol=m)
for (i in 1:m) {
  mult[i] = sum(q[i,]) - 1
  const[i,i] = 1
}



library(slam)
library(gurobi)

model <- list()

model$A          <- const
model$obj        <- t(mult)
model$modelsense <- 'max'
model$rhs        <- numeric(m)+1
model$sense      <- rep(c('<='),m)
model$vtype      <- 'B'

params <- list(OutputFlag=0)

result <- gurobi(model, params)

print('Solution:')
print(result$objval)
print(result$x)
print(sum(result$x)/length(result$x)*100.)