plagCV <- function(B,H,K,lambda,Y){
  n <- length(Y)
  BB <- crossprod(B,B)
  BY <- crossprod(B,Y)
  p <- nrow(H); q <- ncol(H);
  
  nl <- length(lambda)
  gamma_all <- c(); beta_all <- c()
  res_all <- c(); sse_all <- c(); 
  lhs2 <- rbind(cbind(Matrix(0,nrow=p,ncol=p,sparse=TRUE),H),
             Matrix(0,nrow=1,ncol=(p+q),sparse=TRUE))
  rhs <- rbind(BY,Matrix(0,nrow=p+1,ncol=1,sparse=TRUE))
  
  nfold <- 10; sfold <- round(n/nfold); 
  indp <- sample(n,n,replace=FALSE);
  for(ii in 1:nfold){
    ssei <- c()
    if(ii<nfold){
      indi <- sort(indp[((ii-1)*sfold+1):(ii*sfold)])
    }
    if(ii==nfold){
      indi <- sort(indp[((ii-1)*sfold+1):n])
    }
    indn <- setdiff(1:n,indi)
    Bi <- B[indi,]; Bn <- B[indn,];
    Yi <- Y[indi]; Yn <- Y[indn]; 
    BBi <- crossprod(Bi,Bi); BBn <- crossprod(Bn,Bn);
    BYi <- crossprod(Bi,Yi); BYn <- crossprod(Bn,Yn);
    rhsn <- rbind(BYn,Matrix(0,nrow=p+1,ncol=1,sparse=TRUE))
    for(il in 1:nl){
      Lam <- lambda[il]
      Dlam <- Lam*K
      lhs <- rbind(cbind(t(H),BBn+Dlam),lhs2)
      asvd <- svd(lhs)
      tmp1 <- 1/asvd$d
      tmp1[tmp1>1e12] <- 0
      adiag <- diag(tmp1)
      alpha <- asvd$v %*% adiag %*% t(asvd$u) %*% rhsn
      # check <- lhs %*% alpha
      # max(abs(check-rhs))
      gamma <- alpha[-(1:p)]
      beta <- Bi%*%gamma
      sse <- sum((Yi-beta)^2)
      ssei <- c(ssei,sse)
    }
    sse_all <- rbind(sse_all,ssei)
  }
  sse_all <- apply(sse_all,2,sum)
  j <- which.min(sse_all)
  gcv <- sse_all[j]
  lambdac <- lambda[j]
  
  Dlam <- lambdac*K
  lhs <- rbind(cbind(t(H),BB+Dlam),lhs2)
  asvd <- svd(lhs)
  tmp1 <- 1/asvd$d
  tmp1[tmp1>1e12] <- 0
  adiag <- diag(tmp1)
  alpha <- asvd$v %*% adiag %*% t(asvd$u) %*% rhs
  # check <- lhs %*% alpha
  # max(abs(check-rhs))
  gamma <- alpha[-(1:p)]
  beta <- B%*%gamma
  res <- Y-beta
  sse <- sum(res^2)
  
  list(beta_hat=beta,gamma_hat=gamma,theta_hat=NA,
       sse=sse,gcv=gcv,bic=NA,lamc=lambdac,Yhat=beta,
       res=res,edf=NA)
}
