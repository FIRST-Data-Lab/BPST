qrH <- function(H){
  qrout <- qr(t(H),tol=1e-10)
  QH <- qr.Q(qrout,complete=TRUE)
  RH <- qr.R(qrout,complete=TRUE)
  r1 <- qrout$rank
  if(r1>0 & r1<ncol(QH)){
    Q2 <- QH[,-(1:r1)]
  }else{
    Q2 <- c()
  }
  return(Q2)
}
