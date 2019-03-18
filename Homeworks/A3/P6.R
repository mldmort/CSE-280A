
n = 10                               # Size of adjacency matrix
Node_List = list()
for (i in 1:n) {
  Node_List <- c(Node_List,i)
}

# ---------------------------------------------------------------------------------
# Building the adjacency matrix:

Adj_Matrix <- matrix(0, nrow = n, ncol = n)
k = min(n, ceiling(2*log(n, base = 2)))

for (i in 1:n) {
  for (j in i:n) {
    if (i == j) {
      Adj_Matrix[i,j] = 0
    } else if((j - i) == 1 & j<= k) {
      Adj_Matrix[i,j] = sample(c(-3,10), 1)
    } else {
      Adj_Matrix[i,j] = sample(c(-10:3), 1)
    }
  }
}

Adj_Matrix[lower.tri(Adj_Matrix)] = t(Adj_Matrix)[lower.tri(Adj_Matrix)]

# cat("Is adjacency matrix symmetric?", isSymmetric(Adj_Matrix))

# ---------------------------------------------------------------------------------
# Path_Length function calculates length of a path:

Path_Length <- function(T) {
  Len <- 0
  if (length(T) <= 1) {
    return (0)
  }
  for (i in 1:(length(T)-1)) {
    Len <- Len + Adj_Matrix[as.integer(T[i]), as.integer(T[i+1])]
  }
  return (Len)
}

# ---------------------------------------------------------------------------------

Neighbor_List <- function(Path, Node_List) {
  Neighbor <- list()
  Set_Diff <- setdiff(Node_List, Path)
  
  # reduce one path
  for (i in 1:length(Path)) {
    Neighbor <- c(Neighbor, list(Path[-i]))
  }
  
  if (length(Set_Diff)==0) {
    return (Neighbor)
  }
  
  # add one path
  for (j in 1:length(Set_Diff)) {
    Neighbor <- c(Neighbor, list(c(Set_Diff[j],Path)))
  }
  for (i in 1:(length(Path)-1)) {
    for (j in 1:length(Set_Diff)) {
      Neighbor <- c(Neighbor, list(c(Path[1:i],Set_Diff[j],Path[(i+1):length(Path)])))
    }
  }
  for (j in 1:length(Set_Diff)) {
    Neighbor <- c(Neighbor, list(c(Path,Set_Diff[j])))
  }
  
  return(Neighbor)
}

############################## GRIDY ######################################
Path = list(5,  6,  8,)
max_iter = 20

max_cost = -1e20
for (iter in 1:max_iter) {
  my_neighbors = Neighbor_List(Path,Node_List)
  for (i in 1:length(my_neighbors)) {
    cost = Path_Length(my_neighbors[[i]])
    if (cost >= max_cost) {
      max_cost = cost
      Path <- my_neighbors[[i]]
    }
  }
  cat('iteration: ',iter, ' ')
  for (i in 1:length(Path)) {
    cat(Path[[i]],' ')
  }
  cat('\n')
  cat('max cost: ',max_cost,'\n')
}
###########################################################################

############################## Simulation Annealing #######################
Path = list(6,  5,  2,  1,  10,  9,  4,  3)
max_iter = 300
Temp = 1000.

max_cost = -1e20
for (iter in 1:max_iter) {
  my_neighbors = Neighbor_List(Path,Node_List)
  rand_ind = sample(c(1:length(my_neighbors)),1)
  rand_neighbor = my_neighbors[[rand_ind]]
  my_cost = Path_Length(rand_neighbor)
  if (my_cost >= max_cost) {
    max_cost = my_cost
    Path <- my_neighbors[[rand_ind]]
  } else {
     rand_val = runif(1,0,1)
     prob = exp(-(max_cost-my_cost)/Temp)
     cat('prob: ', prob,'\n')
     if (rand_val < prob) {
       max_cost = my_cost
       Path <- my_neighbors[[rand_ind]]
     }
     Temp = Temp*0.9
  }

  cat('iteration: ',iter, ' ')
  for (i in 1:length(Path)) {
    cat(Path[[i]],' ')
  }
  cat('\n')
  cat('max cost: ',max_cost,'\n')
}
###########################################################################