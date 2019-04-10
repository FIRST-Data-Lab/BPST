seval <- function(V,Tr,c,sz,IndP,cnt,Lam){
  k <- dim(Tr)[1]
  Z <- rep(0,sz)
  m <- length(c)/k
  d <- degree(m)
  Ind1 <- matrix(0,nr=1/6*d*(d+1)*(d+2),nc=1)
  Ind2 <- Ind1
  Ind3 <- Ind1
  cnt1 <- matrix(0,nr=d+1,nc=1)
  for (i in 1:d){
    cnt1[i+1] <- cnt1[i]+(d-i+1)*(d-i+2)/2
    IJK <- indices(d-i+1)
    I <- IJK$I
    J <- IJK$J
    K <- IJK$K
    I1J1K1 <- indices(d-i)
    I1 <- I1J1K1$I
    J1 <- I1J1K1$J
    K1 <- I1J1K1$K
    Ind1[(cnt1[i]+1):cnt1[i+1]] <- locate(I1+1,J1,K1,I,J,K)
    Ind2[(cnt1[i]+1):cnt1[i+1]] <- locate(I1,J1+1,K1,I,J,K)
    Ind3[(cnt1[i]+1):cnt1[i+1]] <- locate(I1,J1,K1+1,I,J,K)
  }
  C <- matrix(c,nr=m,nc=k)
  cnt <- as.matrix(cnt)
  Lam <- matrix(Lam,nc=3)
  B <- apply(Lam,2,function(X){X=as.matrix(X); HQblkdiag(X,cnt)})
  B1 <- matrix(B[,1],nr=sz)
  B2 <- matrix(B[,2],nr=sz)
  B3 <- matrix(B[,3],nr=sz)
  C <- C[Ind1[(cnt1[1]+1):cnt1[1+1]],]%*%t(B1)+
    C[Ind2[(cnt1[1]+1):cnt1[1+1]],]%*%t(B2)+
    C[Ind3[(cnt1[1]+1):cnt1[1+1]],]%*%t(B3)
  if(nrow(Lam)==1){
    s1 <- as.matrix(Lam[,1])
    s2 <- as.matrix(Lam[,2])
    s3 <- as.matrix(Lam[,3])
  }else{
    s1 <- diag(Lam[,1])
    s2 <- diag(Lam[,2])
    s3 <- diag(Lam[,3])
  }
  
  s1 <- Matrix(s1,sparse=TRUE)
  s2 <- Matrix(s2,sparse=TRUE)
  s3 <- Matrix(s3,sparse=TRUE)
  if (d>=2){
    for (i in 2:d){
      C <- C[Ind1[(cnt1[i]+1):cnt1[i+1]],]%*%s1+
        C[Ind2[(cnt1[i]+1):cnt1[i+1]],]%*%s2+
        C[Ind3[(cnt1[i]+1):cnt1[i+1]],]%*%s3
    }
  }
  C <- as.matrix(C)
  # IndP <- as.matrix(IndP)
  Z[IndP] <- C
  return(Z)
}
