library(Matrix)
library(schoolmath)
library(MASS)
library(numbers)

homology <- function(degree,k,degenerate){
  boundary_F <- matrix(nrow=12,ncol=6,byrow=T,data=c(1,-1,0,0,1,0,1,-1,-1,0,0,0,-1,1,1,0,0,0,-1,1,0,0,-1,0,0,0,1,-1,0,1,-1,0,1,-1,0,0,0,0,-1,1,0,-1,1,0,-1,1,0,0,0,-1,0,0,1,-1,0,0,0,1,1,-1,0,0,0,-1,-1,1,0,1,0,0,-1,1))#a matrix - Matrix(,sparse=TRUE)
  boundary_G <- matrix(c(-1,0,1,-1,1,0,0,-1,1,1,-1,0,0,1,-1,1,0,-1),ncol=3,nrow=6,byrow=T)#another matrix
  rho <- rankMatrix(boundary_G)[1]
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  X <- findX(boundary_G)
  Z <- X[(rho + 1):q, ]
  B <- GaussianElimination(boundary_F)
  B <- B[which(rowSums(abs(B))!=0), ]
  N <- round(Z %*% ginv(B))
  S <- smith(N)
  Delta <- diag(S)
  s <- length(Delta)
  l <- 0
  ones <- 0
  output <- c()
  for(i in 1:s){
    if(Delta[i]==1){
      ones <- ones + 1
    }
    else if(Delta[i]!=0){
      l <- l + 1
      output <- append(output,Delta[i])
    }
    else{
      break
    }
  }
  
  #the following is the output depending on the number of zeroes.
  if(s>l+ones){
      print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th")))," homology group of R_",k," is isomorphic to Z^",s-(k+ones)," plus the following:"))
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th")))," homology group of R_",k," is isomorphic to the following:"))
  }
  if(l>0){
    for(i in 1:l){
      print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
    }
  }
  else{
    print("0")
  }
}

findX <- function(A){
  X <- diag(nrow(A))
  G_X <- cbind(A,X)
  GX_res <- GaussianElimination(G_X)
  X <- GX_res[, (ncol(GX_res)-nrow(A)+1):ncol(GX_res)]
  return(X)
}

smith <- function(A){
  HF <- hermiteNF(A)$H
  S <- hermiteNF(t(HF))$H
  return(S)
}