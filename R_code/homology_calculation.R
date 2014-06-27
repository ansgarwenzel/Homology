#first load the packages necessary for third-party functions.
library(Matrix)

library(MASS)
library(numbers)
library(compiler)


#first define the subfunctions

#this is used to calculate the left elementary matrix X
findXc <- function(A){
  X <- diag(nrow(A)) #create Identity matrix
  G_X <- cbind(A,X) #combine identity matrix with original matrix for Gaussian elimination
  GX_res <- GaussianElimination(G_X) #Gaussian elimination
  X <- GX_res[, (ncol(GX_res)-nrow(A)+1):ncol(GX_res)] #extract the left matrix, X
  return(X)
}
findX <- cmpfun(findXc)

push_downc <- function(D){
  for(i in 1:(length(D) - 1)){
    a <- D[i]
    b <- D[i + 1]
    if(a!=1||b!=1){
      d <- GCD(a,b)
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
  return(D)
}

push_down <- cmpfun(push_downc)

check_more_pushc <- function(D){
  for(i in 2:length(D)){
    if(D[i] < D[i - 1]){
      if(D[i]!=0){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
  }
  return(FALSE)
}

check_more_push <- cmpfun(check_more_pushc)

#this calculates the smith normal form based on the calculation of the Hermite Normal Form.
#In particular, HNF((HNF(A))^T)
smithc <- function(S){
  S <- hermiteNF(S)$H
  S <- t(hermiteNF(t(S))$H)
  D <- diag(S)
  D <- push_down(D)
  for(i in 1:length(D)){
    D[i] <- ifelse(D[i]<0, -D[i], D[i])
  }
  j <- 1
  more_push <- check_more_push(D)
  while(more_push){
    D <- push_down(D)
    for(i in 1:length(D)){
      D[i] <- ifelse(D[i]<0, -D[i], D[i])
    }
    print(j)
    j <- j+1
    more_push <- check_more_push(D)
  }
  for(i in 1:length(D)){
    D[i] <- ifelse(D[i]<0, -D[i], D[i])
  }
  diag(S) <- D
  return(S)
}

smith <- cmpfun(smithc)

matrix_rankc <- function(A){
  A <- GaussianElimination(A)
  A <- unique(A)
  return(nrow(A) - 1)
}

matrix_rank <- cmpfun(matrix_rankc)

row_spacec <- function(B){
  B <- t(hermiteNF(t(B))$H)
  return(B)
}

row_space <- cmpfun(row_spacec)



#here is the main function to calculate the homology
homologyc <- function(degree, k, quandle=TRUE){
  boundary_F <- boundary_matrix(degree + 1, k, quandle)
  boundary_G <- boundary_matrix(degree, k, quandle)
  rho <- matrix_rank(boundary_G) #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  q_rho <- q - rho
  boundary_G <- findX(boundary_G)
  boundary_G <- boundary_G[(rho + 1):q, ] #only take the rows that map to zero (i.e. only the boundaries)
  boundary_F <- row_space(boundary_F) #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  boundary_G <- round(boundary_F %*% ginv(boundary_G)) #calculate N. Details in documentation
  boundary_G <- smith(boundary_G)
  Delta <- diag(boundary_G) #Extract the values necessary for the output.
  s <- length(Delta)
  l <- 0
  ones <- 0
  output <- c()
  for(i in 1:s){
    if(Delta[i]==1){
      ones <- ones + 1 #count number of ones
    } else if(Delta[i]!=0){
      l <- l + 1
      output <- append(output,Delta[i]) #count and extract nonzero and non-one values in the diagonal
    } else{
      break
    }
  }
  #the following is the output depending on the number of zeroes.
  
  if(s>l+ones){#check if there are any values not equal to one or zero
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), ifelse(quandle," quandle"," rack"), " homology group of R_",k," is isomorphic to Z^", s-(l+ones)," plus the following:"))
  }
  if(s > (l + ones)){#check if there are any values not equal to one or zero
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), ifelse(quandle," quandle"," rack"), " homology group of R_",k," is isomorphic to Z^", s-(l+ones)," plus the following:"))
    
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), ifelse(quandle," quandle"," rack"), " homology group of R_",k," is isomorphic to the following:"))
  }
  if(l>0){ #if so, print out the resulting Z_n groups
    for(i in 1:l){
      print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
    }
  } else{
    print("0")
  }
}

homology <- cmpfun(homologyc)


degenerate_homologyc <- function(degree, k){
  boundary_F <- boundary_matrix_degenerate(degree + 1, k)
  boundary_G <- boundary_matrix_degenerate(degree, k)
  rho <- matrix_rank(boundary_G) #first, this calculates the rank of the matrix G. This removes the need to calculate D and Y later.
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  X <- findX(boundary_G)
  Z <- X[(rho + 1):q, ] #only take the rows that map to zero (i.e. only the boundaries)
  B <- row_space(boundary_F) #identify the cycles, i.e. remove the boundaries via Gaussian elimination
  N <- round(B %*% ginv(Z)) #calculate N. Details in documentation
  S <- smith(N)
  Delta <- diag(S) #Extract the values necessary for the output.
  s <- length(Delta)
  l <- 0
  ones <- 0
  output <- c()
  for(i in 1:s){
    if(Delta[i]==1){
      ones <- ones + 1 #count number of ones
    } else if(Delta[i]!=0){
      l <- l + 1
      output <- append(output,Delta[i]) #count and extract nonzero and non-one values in the diagonal
    } else{
      break
    }
  }
  
  #the following is the output depending on the number of zeroes.
  if(s>l+ones){#check if there are any values not equal to one or zero
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), " degenerate", " homology group of R_",k," is isomorphic to Z^", s-(l+ones)," plus the following:"))
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))), " degenerate", " homology group of R_",k," is isomorphic to the following:"))
  }
  if(l>0){ #if so, print out the resulting Z_n groups
    for(i in 1:l){
      print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
    }
  } else{
    print("0")
  }
}

degenerate_homology <- cmpfun(degenerate_homologyc)