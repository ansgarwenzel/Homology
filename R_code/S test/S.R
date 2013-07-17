all_combinationsc <- function(k){
  output <- t(combn(rep(0:(k-1), 2), 2))
  output <- unique(output)
}

all_combinations <- cmpfun(all_combinationsc)

X_squared <- all_combinations(k)
S_X <- X_squared
#here, you have to change S_X to define S

check_permutationsc <- function(B){
  C <- unique(B)
  if(C!=B){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
check_permutations <- cmpfun(check_permutationsc)

