library(Matrix)

homology <- function(degree,k,degenerate){
  boundary_F <- matrix(nrow=12,ncol=6,byrow=T,data=c(1,-1,0,0,1,0,1,-1,-1,0,0,0,-1,1,1,0,0,0,-1,1,0,0,-1,0,0,0,1,-1,0,1,-1,0,1,-1,0,0,0,0,-1,1,0,-1,1,0,-1,1,0,0,0,-1,0,0,1,-1,0,0,0,1,1,-1,0,0,0,-1,-1,1,0,1,0,0,-1,1))#a matrix - Matrix(,sparse=TRUE)
  boundary_G <- matrix(c(-1,0,1,-1,1,0,0,-1,1,1,-1,0,0,1,-1,1,0,-1),ncol=3,nrow=6,byrow=T)#another matrix
  rho <- rankMatrix(boundary_G)[1]
  q <- nrow(boundary_G)
  r <- ncol(boundary_G)
  q_rho <- q - rho
  GX <- find_G_and_X(boundary_G)
  D <- GX[[1]][(q_rho + 1):q,(r - rho + 1):r]
  Z <- matrix(nrow=q,ncol=q,rep.int(0,q*q))
  Z <- GX[[2]][1:(q_rho),1:(q_rho)]
  B <- find_B(boundary_F)
  N <- B*solve(Z)
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
  if(s>l+ones){
      print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th")))," homology group of R_",k," is isomorphic to Z^",s-(k+ones)," plus the following:"))
  } else{
    print(paste0("The ",degree,ifelse((degree%%10)==1,"st",ifelse((degree%%10)==2,"nd",ifelse((degree%%10)==3,"rd","th"))),"-th homology group of R_",k," is isomorphic the following:\n"))
  }
  for(i in 1:l){
    print(paste0("Z_",Delta[ones+i],ifelse(i!=l," plus","")))
  }
}

find_G_and_X <- function(G)
{
  
}

find_B <- function(M)
{
  
}

smith <- function(A)  #calculates smith form of a matrix
{
  
}