#this collects all the functions necessary for the calculation of boundary matrices of (bi)quandles

up_action <- function(a, b, k){
  result <- (2 * b - a) %% k
  return(as.integer(result))
}

down_action <- function(a, b, k){
  return(as.integer(a))
}

boundary_names <- function(degree,k,degenerate){
  output <- t(combn(rep(0:(k-1), degree), degree))
  output <- unique(output)
  if(degenerate&&ncol(output)>1){
    for(i in nrow(output):1){
      for(j in 2:ncol(output)){
        if(output[i, j]==output[i, j - 1]){
          output <- output[-i, ]
          break         
        }
      }
    }
  }
  return(output)
}

boundary_matrix <- function(degree, k, degenerate=FALSE){
  if(degenerate){
    m <- k*((k-1)^(degree-1))
    n <- k*((k-1)^(degree-2))
  } else{
    m <- k^(degree)
    n <- k^(degree-1)
  }
  M <- matrix(ncol=n,nrow=m,0)
  column_names <- boundary_names(degree - 1,k,degenerate)
  row_names <- boundary_names(degree,k,degenerate)
  
  for(i in 1:nrow(row_names)){
    result_vector <- rep(0,n)
    for(j in 1:ncol(row_names)){
      name_row <- row_names[i, ]
      b <- name_row[j]
      name_row <- name_row[-j]
      no_double_names <- TRUE
      
      if(degenerate&&length(name_row) > 1){
        for(l in 2:length(name_row)){
          if(name_row[l] == name_row[l - 1]){
            no_double_names <- FALSE
            break
          }
        }
      }
      
      if(no_double_names){
        row.is.a.match <- apply(column_names, 1, identical, name_row)
        match.id <- which(row.is.a.match)
        
        if(j %% 2 == 0){
          result_vector[match.id] <- result_vector[match.id] + 1
        } else{
          result_vector[match.id] <- result_vector[match.id] - 1
        }
      }
      
      no_double_names <- TRUE
      if(j > 1){
        name_row[1:(j - 1)] <- up_action(name_row[1:(j - 1)], b, k)
      }
      
      if(j < ncol(column_names)){
        name_row[j:length(name_row)] <- down_action(name_row[j:length(name_row)], b, k)
      }

      if(degenerate&&length(name_row) > 1){
        for(l in 2:length(name_row)){
          if(name_row[l] == name_row[l - 1]){
            no_double_names <- FALSE
            break
          }
        }
      }
      
      if(no_double_names){
        row.is.a.match <- apply(column_names, 1, identical, name_row)
        match.id <- which(row.is.a.match)
        
        if(j %% 2 == 0){
          result_vector[match.id] <- result_vector[match.id] - 1
        } else{
          result_vector[match.id] <- result_vector[match.id] + 1
        }
      }
    }
    M[i, ] <- result_vector
  }
  return(M)
}