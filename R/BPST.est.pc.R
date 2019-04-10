BPST.est.pc <- function(B,Y){
  n <- length(Y)
  B <- as.matrix(B)
  J <- ncol(B)
  WW <- crossprod(B,B)
  rhs <- crossprod(B,Y)
  theta <- solve(WW,rhs)
  gamma <- theta
  beta <- B%*%gamma
  res <- Y-beta
  list(beta_hat=beta,gamma_hat=gamma,theta_hat=theta,
       sse=sse,gcv=NULL,lamc=NULL,Yhat=beta,
       res=res,lamc=NULL,edf=NULL)
}
