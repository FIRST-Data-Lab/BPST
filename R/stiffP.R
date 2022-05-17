stiffP <- function(V, Tr, d){
  V = as.matrix(V)
  Tr = as.matrix(Tr)
  n = nrow(Tr)
  m = (d + 1) * (d + 2) / 2
  msqr <- m * m
  Mat <- build(d - 1)
  Indx1 <- rep(0, n * msqr)
  Indx2 <- Indx1
  S = Indx1
  place <- 1
  for(k in 1:n){
    LocS <- locStiff(V[Tr[k, 1], ], V[Tr[k, 2], ], V[Tr[k, 3], ], Mat, d)
    row <- c()
    for(j in 1:ncol(LocS)){
      for(i in 1:nrow(LocS)){
        if(LocS[i, j] != 0){
          row <- c(row, i)
        }
      }
    }
    
    col <- c()
    for(j in 1:ncol(LocS)){
      for(i in 1:nrow(LocS)){
        if(LocS[i, j] != 0){
          col <- c(col, j)
        }
      }
    }
    
    val <- c()
    for(j in 1:ncol(LocS)){
      for(i in 1:nrow(LocS)){
        if(LocS[i, j] != 0){
          val <- c(val, LocS[i, j])
        }
      }
    }
    
    L <- length(row)
    Indx1[place:(place + (L - 1))] = (k - 1) * m + row
    Indx2[place:(place + (L - 1))] = (k - 1) * m + col
    S[place:(place + (L - 1))] = val
    place = place + L
  }
  
  K <- sparseMatrix(i = Indx1[1:(place - 1)],
                             j = Indx2[1:(place - 1)], 
                             x = S[1:(place - 1)], 
                             dims = c(n*m, n*m))
  return(K)
}
