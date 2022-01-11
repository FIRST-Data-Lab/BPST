abary <- function(Tr0,V0,xx,yy){
  xx <- as.vector(xx); yy <- as.vector(yy);
  n1 <- length(xx); n2 <- length(yy);
  if(n1!=n2) warning("The length of the coordinates do not match!")
  n <- min(n1,n2)
  ind <- 1:n
  tol <- -1e-12
  V1 <- V0[Tr0[1],]
  V2 <- V0[Tr0[2],]
  V3 <- V0[Tr0[3],]
  lam <- bary(V1,V2,V3,xx,yy)
  J <- ind[lam$lam1>tol & lam$lam2>tol & lam$lam3>tol]
  return(J)
}

B0_Generator <- function(V0,Tr0,Z){
  row.names(Tr0) <- 1:nrow(Tr0)
  Z <- matrix(Z,ncol=2)
  n <- nrow(Z)
  J <- apply(Tr0,1,abary,V0=V0,xx=Z[,1],yy=Z[,2])
  x1 <- unlist(J)
  ind2 <- as.numeric(names(J))
  l <- sapply(J,length)
  x2 <- rep(ind2,l)
  B <- sparseMatrix(i=x1,j=x2,x=rep(1,length(x1)),dims=c(n,nrow(Tr0)))
  B <- B[apply(B,1,function(x) !all(x==0)),]
  
  list(B=B,ind.inside=sort(unique(x1)))
}
