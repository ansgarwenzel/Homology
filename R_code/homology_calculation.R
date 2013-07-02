#first load the packages necessary for third-party functions.
library(Matrix)
library(schoolmath)
library(MASS)
library(numbers)

#first define the subfunctions

#this is used to calculate the left elementary matrix X
findX <- function(A){
  X <- diag(nrow(A))    #create Identity matrix
  G_X <- cbind(A,X)     #combine identity matrix with original matrix for Gaussian elimination
  GX_res <- GaussianElimination(G_X)  #Gaussian elimination
  X <- GX_res[, (ncol(GX_res)-nrow(A)+1):ncol(GX_res)]  #extract the left matrix, X
  return(X)
}

#this calculates the smith normal form based on the calculation of the Hermite Normal Form.
#In particular, HNF((HNF(A))^T)
smith <- function(A){
  HF <- hermiteNF(A)$H
  S <- hermiteNF(t(HF))$H
  return(S)
}
#here is the main function to calculate the homology
homology <- function(degree, k, degenerate){
  #boundary_F <- matrix(nrow=12,ncol=6,byrow=T,data=c(1,-1,0,0,1,0,1,-1,-1,0,0,0,-1,1,1,0,0,0,-1,1,0,0,-1,0,0,0,1,-1,0,1,-1,0,1,-1,0,0,0,0,-1,1,0,-1,1,0,-1,1,0,0,0,-1,0,0,1,-1,0,0,0,1,1,-1,0,0,0,-1,-1,1,0,1,0,0,-1,1))#a matrix - Matrix(,sparse=TRUE)
  #boundary_G <- matrix(c(-1,0,1,-1,1,0,0,-1,1,1,-1,0,0,1,-1,1,0,-1),ncol=3,nrow=6,byrow=T)#another matrix
  boundary_F <- boundary_matrix(degree + 1, k, degenerate)
  boundary_G <- boundary_matrix(degree, k, degenerate)
  rho <- rankMatrix(boundary_G)[1]  #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  X <- findX(boundary_G)
  Z <- X[(rho + 1):q, ]       #only take the rows that map to zero (i.e. only the boundaries)
  B <- GaussianElimination(boundary_F)  #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  #B <- rref(boundary_F)
  N <- round(Z %*% ginv(B)) #calculate N. Details in documentation
  S <- smith(N)
  Delta <- diag(S)  #Extract the values necessary for the output.
  s <- length(Delta)
  l <- 0
  ones <- 0
  output <- c()
  for(i in 1:s){
    if(Delta[i]==1){
      ones <- ones + 1    #count number of ones
    } else if(Delta[i]!=0){
      l <- l + 1
      output <- append(output,Delta[i])    #count and extract nonzero and non-one values in the diagonal
    } else{
      break
    }
  }
  
  #the following is the output depending on the number of zeroes.
  if(s>l+ones){#check if there are any values not equal to one or zero
      print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), ifelse(degenerate," quandle"," rack"), " homology group of R_",k," is isomorphic to Z^", s-(l+ones)," plus the following:"))
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), ifelse(degenerate," quandle"," rack"), " homology group of R_",k," is isomorphic to the following:"))
  }
  if(l>0){ #if so, print out the resulting Z_n groups
    for(i in 1:l){
      print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
    }
  } else{
    print("0")
  }
}