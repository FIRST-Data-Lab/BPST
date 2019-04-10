tria <- function(Tr0,V0){
  VT <- V0[Tr0,]
  V1 <- VT[1,]; V2 <- VT[2,]; V3 <- VT[3,];
  triall <- triarea(V1,V2,V3)
  return(triall)
}

TArea <- function(V,Tr,Z){
  row.names(Tr) <- 1:nrow(Tr)
  Z <- matrix(Z,ncol=2)
  n <- nrow(Z)
  J <- apply(Tr,1,abary,V0=V,xx=Z[,1],yy=Z[,2])
  ind2 <- as.numeric(names(J))
  l <- sapply(J,length)
  x2 <- rep(ind2,l)
  triall <- apply(Tr,1,tria,V0=V)
  tria.all <- triall[x2]
  return(tria.all)
}
