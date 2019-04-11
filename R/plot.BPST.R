#' Produces the contour plot for the estimated surface of a fitted "BPST" object.
#'
#' This function produces the contour map for the estimated surface of a fitted "BPST" object.
#'
#' @importFrom Matrix Matrix
#' @importFrom pracma isempty
#' @importFrom Rcpp evalCpp
#' 
#' @param mfit Fitted ``BPST" object.
#' \cr
#' @param Zgrid The grid points used to construct the contour plot.
#' \cr
#' @return None
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @export
#'
plot.BPST <- function(mfit,Zgrid=NULL){
  if(isempty(Zgrid)){
    z1.min <- min(mfit$Z[,1])-0.1; z1.max <- max(mfit$Z[,1])+0.1;
    z2.min <- min(mfit$Z[,2])-0.1; z2.max <- max(mfit$Z[,2])+0.1;
    n1 <- 50; n2 <- 50;
    u1 <- seq(z1.min,z1.max,length.out=n1)
    v1 <- seq(z2.min,z2.max,length.out=n2)
    uu <- rep(u1,each=n2)
    vv <- rep(v1,times=n1)
    uu.mtx <- matrix(uu,n2,n1)
    vv.mtx <- matrix(vv,n2,n1)
    Zgrid <- as.matrix(cbind(uu,vv))
  }else{
    u1 <- sort(unique(Zgrid[,1]))
    v1 <- sort(unique(Zgrid[,2]))
    n1 <- length(u1)
    n2 <- length(v1)
  }
  Bfull <- basis(mfit$V,mfit$Tr,mfit$d,mfit$r,Zgrid,FALSE,FALSE,FALSE,FALSE)
  B <- Bfull$B
  Ind.inside <- Bfull$Ind.inside
  Ygrid <- rep(NA,nrow(Zgrid))
  Ygrid[Ind.inside] <- B%*%mfit$gamma_hat
  Ygrid.mtx <- matrix(Ygrid,ncol=n1,nrow=n2)
  image(v1,u1,Ygrid.mtx)
  contour(Ygrid.mtx,add=TRUE,method="edge",vfont=c("sans serif","plain"))
}
