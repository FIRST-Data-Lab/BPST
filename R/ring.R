ring <- function(Tr0, V, Tr, TV, n.layer = 2){
  I.tr = c(Tr0)
  V1 = c(); T1 = c();
  for(i in 1:n.layer){
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
  
  ring.tri = list(V1 = V1, T1 = T1, Tr0 = Tr0, I.tr = I.tr, T.sub = T.sub)
}