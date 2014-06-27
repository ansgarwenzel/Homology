library(compiler)


up_actionc <- function(a, b, k){
  
  #result <- 2 * b - a ###for dihedral quandle
  #binv <- k - b
  #result <- binv * a * b
  #result <- result %% k
  
  result <- (2 * b - a) %% k ###dihedral quandle
  #result <- sth ### alexander quandle
  ####commutative quandle
  #action_matrix <- rbind(c(0,0,0,0,0,0),c(1,1,5,5,2,2),c(2,5,2,1,5,1),c(3,4,4,3,4,4),c(4,3,3,3,4,3),c(5,2,1,2,1,5))
  #result <-action_matrix[a + 1, b + 1]
  ############
  
  return(as.integer(result))
}

up_action <- cmpfun(up_actionc)

down_actionc <- function(a, b, k){
  return(as.integer(a))
}

down_action <- cmpfun(down_actionc)


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


S_test <- function(k){
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
}