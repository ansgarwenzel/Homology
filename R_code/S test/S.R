library(compiler)

all_combinationsc <- function(k){
  output <- t(combn(rep(0:(k-1), 2), 2))
  output <- unique(output)
}

all_combinations <- cmpfun(all_combinationsc)


check_permutationsc <- function(B){
  C <- unique(B)
  if(C!=B){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
check_permutations <- cmpfun(check_permutationsc)

check_YBc <- function(B,k){
  for(i in 1:nrow(S_X)){
    
  }
  LHS <- 
  RHS <- 
    
  return_value <- (LHS==RHS)
  return(return_value)
}
check_YB <- cmpfun(check_YBc)

k <- 3#define k for Z_k (or X_k) where S lives in...
X_squared <- all_combinations(k)
S_X <- X_squared
#here, you have to change S_X to define S

#then, check that permutations hold:
permutations <- check_permutations(S_X)

#and check that Yang-Baxter holds
Yang_Baxter <- check_YB(S_X,k)