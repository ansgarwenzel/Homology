#this collects all the functions necessary for the calculation of boundary matrices of (bi)quandles

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

boundary_namesc <- function(degree,k,degenerate){
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

boundary_names <- cmpfun(boundary_namesc)

boundary_matrixc <- function(degree, k, degenerate=FALSE){
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

boundary_matrix <- cmpfun(boundary_matrixc)

boundary_matrix_degeneratec <- function(degree, k){
  
  m <- k^(degree) - k*((k-1)^(degree-1))
  n <- k^(degree-1) - k*((k-1)^(degree-2))
  
  M <- matrix(ncol=n,nrow=m,0)
  column_names <- boundary_names_degenerate(degree - 1,k)
  row_names <- boundary_names_degenerate(degree,k)
  
  for(i in 1:nrow(row_names)){
    result_vector <- rep(0,n)
    for(j in 1:ncol(row_names)){
      name_row <- row_names[i, ]
      b <- name_row[j]
      name_row <- name_row[-j]
      double_names <- FALSE
      
      if(length(name_row) > 1){
        for(l in 2:length(name_row)){
          if(name_row[l] == name_row[l - 1]){
            double_names <- TRUE
            break
          }
        }
      }
      
      if(double_names){
        row.is.a.match <- apply(column_names, 1, identical, name_row)
        match.id <- which(row.is.a.match)
        
        if(j %% 2 == 0){
          result_vector[match.id] <- result_vector[match.id] + 1
        } else{
          result_vector[match.id] <- result_vector[match.id] - 1
        }
      }
      
      double_names <- FALSE
      if(j > 1){
        name_row[1:(j - 1)] <- up_action(name_row[1:(j - 1)], b, k)
      }
      
      if(j < ncol(column_names)){
        name_row[j:length(name_row)] <- down_action(name_row[j:length(name_row)], b, k)
      }
      
      if(length(name_row) > 1){
        for(l in 2:length(name_row)){
          if(name_row[l] == name_row[l - 1]){
            double_names <- TRUE
            break
          }
        }
      }
      
      if(double_names){
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

boundary_matrix_degenerate <- cmpfun(boundary_matrix_degeneratec)

boundary_names_degeneratec <- function(degree,k){
  output <- t(combn(rep(0:(k-1), degree), degree))
  output <- unique(output)
  
  for(i in nrow(output):1){
    keep <- FALSE
    for(j in 2:ncol(output)){
      if(output[i, j]==output[i, j - 1]){
        keep <- TRUE
        break
      }
    }
    if(!keep){
      output <- output[-i, ]
    }
  }
  return(output)
}

boundary_names_degenerate <- cmpfun(boundary_names_degeneratec)