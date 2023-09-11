ring2 <- function(Tr0, V, Tr, TV, nl = 1){
  I.tr = c(Tr0)
  V1 = c(); T1 = c();
  if(nl == 0) {
    V1 = V[Tr0, ]; 
    T1 = matrix(1:3, ncol = 3);
    I.tr = row.match(Tr0, Tr);
    T.sub = matrix(Tr0, ncol = 3);
    Tr0 = T1; 
  } else { 
    for(i in 1:nl){
      k = length(I.tr)
      J = c()
      for(j in 1:k){
        J1 = which(TV[, I.tr[j]] != 0)
        J = c(J, J1)
      }
      V2 = Tr[J, ]
      I.tr = c(V2)
      V1 = c(V1, I.tr)
      T1 = c(T1, J)
    }
    V1 = sort(unique(V1))
    I.tr = V1
    V1 = V[I.tr, ]
    T1 = Tr[sort(unique(T1)), ]
    T.sub = T1
    nt = nrow(T1)
    for(i in 1:nt){
      j = T1[i, 1]; nj = which(I.tr == j); T1[i, 1] = nj
      j = T1[i, 2]; nj = which(I.tr == j); T1[i, 2] = nj
      j = T1[i, 3]; nj = which(I.tr == j); T1[i, 3] = nj
    }
    j = Tr0[1]; nj = which(I.tr == j); Tr0[1] = nj
    j = Tr0[2]; nj = which(I.tr == j); Tr0[2] = nj
    j = Tr0[3]; nj = which(I.tr == j); Tr0[3] = nj
  }
  
  ring.tri = list(V1 = V1, T1 = T1, Tr0 = Tr0, I.tr = I.tr, T.sub = T.sub)
}