locStiff <- function(V1, V2, V3, Mat, d){
  m = (d + 1) * (d + 2) / 2
  Id <- diag(m)
  vx <- c(1, 0)
  vy <- c(0, 1)
  lam1x <- tcord(V1, V2, V3, vx)[1]
  lam2x <- tcord(V1, V2, V3, vx)[2]
  lam3x <- tcord(V1, V2, V3, vx)[3]
  lam1y <- tcord(V1, V2, V3, vy)[1]
  lam2y <- tcord(V1, V2, V3, vy)[2]
  lam3y <- tcord(V1, V2, V3, vy)[3]
  Dx <- dirder(Id, lam1x, lam2x, lam3x)
  Dy <- dirder(Id, lam1y, lam2y, lam3y)
  K <- abs(triarea(V1, V2, V3)) * (t(Dx) %*% Mat %*% Dx + t(Dy) %*% Mat %*% Dy)
  return(K)
}
