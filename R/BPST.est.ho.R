BPST.est.ho <- function(B,Q2,K,lambda,Y){
  n <- length(Y)
  J <- ncol(Q2)
  
  W <- as.matrix(B%*%Q2)
  WW <- crossprod(W,W)
  rhs <- crossprod(W,Y)
  D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW)<J)
  if(!flag){
    Ainv <- chol(WW,pivot=TRUE)
    A <- solve(t(Ainv))
    ADA <- A%*%D%*%t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  theta_all <- c(); gamma_all <- c(); beta_all <- c()
  res_all <- c(); sse_all <- c(); df_all <- c(); gcv_all <- c(); bic_all <- c()
  for(il in 1:nl){
    Lam <- lambda[il]
    Dlam <- Lam*D
    lhs <- WW+Dlam
    # chol2inv(chol()) < solve() < qr.solve(); 
    lhs.inv <- chol2inv(chol(lhs));
    # crossprod() < %*%;
    theta <- crossprod(t(lhs.inv),rhs)
    theta_all <- cbind(theta_all,theta)
    gamma <- crossprod(t(Q2),theta)
    gamma_all <- cbind(gamma_all,gamma)
    beta <- crossprod(t(W),theta)
    beta_all <- cbind(beta_all,beta)
    res <- Y-beta
    res_all <- cbind(res_all,res)
    sse <- sum(res^2)
    sse_all <- c(sse_all,sse)
    if(!flag){
      df <- sum(1/(1+Cval*Lam))
    }
    if(flag){
      Hmtx <- crossprod(t(crossprod(t(W),lhs.inv)),t(W))
      df <- sum(diag(Hmtx))
    }
    df_all <- c(df_all,df)
    gcv <- n*sse/(n-df)^2
    gcv_all <- c(gcv_all,gcv)
    bic <- log(sse/n)+df*log(n)/n
    bic_all <- c(bic_all,bic)
  }
  j <- which.min(gcv_all)
  
  lambdac <- lambda[j]
  theta <- theta_all[,j]
  gamma <- gamma_all[,j]
  beta <- beta_all[,j]
  df <- df_all[j]
  sse <- sse_all[j]
  gcv <- gcv_all[j]
  bic <- bic_all[j]
  list(beta_hat=beta,gamma_hat=gamma,theta_hat=theta,
       sse=sse,gcv=gcv,bic=bic,lamc=lambdac,Yhat=beta,
       res=res,edf=df)
}
