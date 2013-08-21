library(compiler)

all_combinationsc <- function(k){
  output <- t(combn(rep(0:(k-1), 2), 2))
  output <- unique(output)
}

all_combinations <- cmpfun(all_combinationsc)

check_permutationsc <- function(B){
  C <- unique(B)
  if(all(C!=B)){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
check_permutations <- cmpfun(check_permutationsc)

find_S_resultc <- function(x,X,S){
  for(i in 1:nrow(X)){
    if(all(X[i, ]==x)){
      match.id <- i
      break
    }
  }
  return(S[match.id, ])
}
find_S_result <- cmpfun(find_S_resultc)

check_YBc <- function(S,k,X){
  for(i in 0:(k - 1)){
    for(j in 1:nrow(S)){
      LHS <- c(i, S[j, ])
      RHS <- LHS
      LHS[1:2] <- find_S_result(LHS[1:2],X,S)
      LHS[2:3] <- find_S_result(LHS[2:3],X,S)
      LHS[1:2] <- find_S_result(LHS[1:2],X,S)
      RHS[2:3] <- find_S_result(RHS[2:3],X,S)
      RHS[1:2] <- find_S_result(RHS[1:2],X,S)
      RHS[2:3] <- find_S_result(RHS[2:3],X,S)
      if(!all(LHS == RHS)){
        return(FALSE)
      }
    }
  }
    return(TRUE)
}
check_YB <- cmpfun(check_YBc)

check_fc <- function(S_X,k,X_squared){
  S <- unique(cbind(X_squared[, 1], S_X[, 1:2]))
    if(nrow(S)!=(k*k)){
      return(FALSE)
    } else{
      return(TRUE)
    }
}

check_f <- cmpfun(check_fc)


k <- 3#define k for Z_k (or X_k) where S lives in...
X_squared <- all_combinations(k)
S_X <- X_squared[, 2:1]
#here, you have to change S_X to define S
S_X[, 2] <- up_action(X_squared[, 1], S_X[, 1], k)
X_squared[, 2] <- down_action(S_X[, 1], X_squared[, 1], k)


#then, check that permutations hold:
permutations_S <- check_permutations(S_X)
permutations_f <- check_f(S_X, k, X_squared)
permutations_g <- check_f(X_squared, k, S_X) #ignore if not a biquandle!

#and check that Yang-Baxter holds, based on S and f operation
Yang_Baxter <- check_YB(S_X,k,X_squared)


print(paste0("The permutation checks hold that S is ", permutations_S, ", f is ",permutations_f," and g is ", permutations_g, " and that the Yang-Baxter check holds ", Yang_Baxter, "."))