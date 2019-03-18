
Input_Data <- as.matrix(read.table(text = gsub("", ' ',
              readLines("/Users/miladmortazavi/Documents/Dev/BIOINF/CSE-280A/Homeworks/A3/a3data2.txt"))))

f = 0.7

colnames(Input_Data) <- c(1:dim(Input_Data)[2])
Sum_Index <- matrix(NA, nrow = dim(Input_Data)[2], ncol = 2)

for (i in 1:dim(Input_Data)[2]) {
  Sum_Index[i,1] = sum(Input_Data[,i])
  Sum_Index[i,2] = i
}

Sorted_Index <- Sum_Index[order(Sum_Index[,1], decreasing = TRUE), ]
NCol = which(Sorted_Index[,1] >= nrow(Input_Data)*f)

Mat <- Input_Data[,Sorted_Index[NCol,2]]

TMP <- list()
Sol <- which(Mat[,1] == 1)
Max <- nrow(Input_Data)*f

i = 1
iter = 0
Track <- list()
Updated_Mat <- Mat[,-1]

for (j in 1:(dim(Mat)[2]-2)) {
  i = 1
  Max_Local <- 0
  ideal <- list()
  while (i < dim(Updated_Mat)[2]) {
    TMP <- which(Updated_Mat[,i] == 1)
    Num <- length(intersect(TMP, Sol))
    if (Num > Max & Num > Max_Local) {
      Max_Local <- Num
      ideal <- TMP
      index_del = i
      i = i + 1
    }
    i = i + 1
  }
  
  if (Max_Local > 0) {
    iter = iter + 1
    Track <- c(Track, colnames(Updated_Mat)[index_del])
    cat("iter:", iter, "\n")
  }
  
  if (length(ideal) > 0) {
    Sol <- intersect(Sol, ideal)
  }
  
  Updated_Mat <- Updated_Mat[,-index_del]
}

Out_Mat = Mat[Sol,unlist(Track)]

cat("Number of columns :", dim(Out_Mat)[2], "\n", 
    "Set of columns :", colnames(Out_Mat), "\n")

cat("Sum of output matrix:", sum(Out_Mat), "\n",
    "Product of nrows and ncols of output matrix:", nrow(Out_Mat)*ncol(Out_Mat), "\n")
