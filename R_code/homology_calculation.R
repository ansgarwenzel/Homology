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
  S <- t(hermiteNF(t(HF))$H)
  D <- diag(S)
  for(i in 1:(length(D) - 1)){
      a <- D[i]
      b <- D[i + 1]
      if(a!=1||b!=1){
        d <- gcd(a,b)
        if(a==0 || b==0){
          d <- ifelse(a==0,abs(b),abs(a))
        }
        if(d!=0){
          alpha <- a/d
          D[i] <- d;
          D[i + 1] <- -(b * alpha)
        } else{
          D[i] <- 0
          D[i+1] <- 0
        }
      }
    }
  for(i in 1:length(D)){
    D[i] <- ifelse(D[i]<0, -D[i], D[i])
  }
  diag(S) <- D
  return(S)
}

i <- i+1

matrix_rank <- function(A){
  A <- GaussianElimination(A)
  A <- unique(A)
  return(nrow(A) - 1)
}

row_space <- function(B){
  B <- t(hermiteNF(t(B))$H)
  return(B)
}


#here is the main function to calculate the homology
homology <- function(degree, k, degenerate){
  #boundary_F <- amatrix(nrow=12,ncol=6,byrow=T,data=c(1,-1,0,0,1,0,1,-1,-1,0,0,0,-1,1,1,0,0,0,-1,1,0,0,-1,0,0,0,1,-1,0,1,-1,0,1,-1,0,0,0,0,-1,1,0,-1,1,0,-1,1,0,0,0,-1,0,0,1,-1,0,0,0,1,1,-1,0,0,0,-1,-1,1,0,1,0,0,-1,1))#a matrix - Matrix(,sparse=TRUE)
  #boundary_G <- matrix(c(-1,0,1,-1,1,0,0,-1,1,1,-1,0,0,1,-1,1,0,-1),ncol=3,nrow=6,byrow=T)#another matrix
  boundary_F <- boundary_matrix(degree + 1, k, degenerate)
  boundary_G <- boundary_matrix(degree, k, degenerate)
  rho <- matrix_rank(boundary_G)  #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  X <- findX(boundary_G)
  Z <- X[(rho + 1):q, ]       #only take the rows that map to zero (i.e. only the boundaries)
  B <- row_space(boundary_F)  #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  N <- round(B %*% ginv(Z)) #calculate N. Details in documentation
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