#' Model Fitting using Bivariate Penalized Spline over Triangulation
#'
#' This function conducts the model fitting via bivariate penlized spline over triangulation.
#'
#' @importFrom Matrix Matrix
#' 
#' @param Y The response variable observed over the domain.
#' \cr
#' @param Z The cooridinates of dimension \code{n} by two. Each row is the coordinates of a point.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 5, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param lambda The tuning parameter -- default is \eqn{10^(-6,-5.5,-5,\ldots,5,5.5,6)}.
#' \cr
#' @param Hmtx The indicator of whether the smoothness matrix \code{H} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param Kmtx The indicator of whether the energy matrix \code{K} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param QR The indicator of whether a QR decomposition need to be performed on the smoothness matrix -- default is \code{TRUE}.
#' \cr
#' @param TA The indicator of whether the area of the triangles need to be calculated -- default is \code{TRUE}.
#' \cr
#' @return A list of vectors and matrice, including:
#' \item{gamma_hat}{The estimated spline coefficients.}
#' \item{lamc}{The tuning parameter selected by Generalized Cross Validation (GCV).}
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+1)(d+2)/2}}, where \code{n} is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. The length of points means the length of ordering indices of observation points. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.}
#' \item{Ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{H}{The smoothness matrix.}
#' \item{Q2}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H}.}
#' \item{K}{The thin-plate energy function.}
#' \item{tria.all}{The area of each triangle within the given triangulation.}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia ang Li Wang from the Iowa State University.
#'
#' @examples
#' # Triangulation
#' # Option 1;
#' # data(V1); data(Tr1); d=5; r=1; V=V1; Tr=Tr1;
#' # Option 2
#' data(V2); data(Tr2); d=5; r=1; V=V2; Tr=Tr2;
#' d=5; r=1;
#' # Grid Points
#' n1.grid=101; n2.grid=101; n.grid=n1.grid*n2.grid;
#' u.grid=seq(0,1,length.out=n1.grid)
#' v.grid=seq(0,1,length.out=n2.grid)
#' uu.grid=rep(u.grid,each=n2.grid)
#' vv.grid=rep(v.grid,times=n1.grid)
#' Z.grid=as.matrix(cbind(uu.grid,vv.grid))
#' func=1; sigma=0.1;
#' gridpoints=DataGenerator(Z.grid,V,Tr,func,sigma,2019)
#' Y.grid=gridpoints$Y; mu.grid=gridpoints$mu;
#' ind=gridpoints$ind; ind.grid=(1:n.grid)[ind==1];
#' # Simulation parameters
#' n=2000;
#' ind.sam=sort(sample(ind.grid,n))
#' Y=as.matrix(gridpoints$Y[ind.sam]); Z=as.matrix(gridpoints$Z[ind.sam,]);
#' mu=as.matrix(gridpoints$mu[ind.sam]);
#' mfit=fit.BPST(Y,Z,V,Tr,d,r,lambda=10^seq(-6,6,by=0.5))
#' rmse=sqrt(mean((Y-mfit$Yhat)^2,na.rm=TRUE))
#' mpred=pred.BPST(mfit,Z.grid)
#' rmspe=sqrt(mean((Y.grid-mpred$Ypred)^2,na.rm=TRUE))
#' cat("rmse =",rmse,"and rmspe =",rmspe,"\n")
#' plot(mfit,Z.grid)
#' @export

fit.BPST <- function(Y,Z,V,Tr,d=5,r=1,lambda=10^seq(-6,6,by=0.5),
                     Hmtx=TRUE,Kmtx=TRUE,QR=TRUE,TA=TRUE){
  this.call <- match.call()
  
  n <- length(Y)
  ind.nna <- (1:n)[!is.na(Y)]
  Z.nna <- Z[ind.nna,]
  Y.nna <- Y[ind.nna]
  Bfull <- basis(V,Tr,d,r,Z.nna,Hmtx,Kmtx,QR,TA)
  B <- Bfull$B
  H <- Bfull$H
  Q2 <- Bfull$Q2
  K <- Bfull$K
  Ind.inside <- Bfull$Ind.inside
  Ind.inside <- ind.nna[Ind.inside]
  tria.all <- Bfull$tria.all
  Yi <- Y[Ind.inside]
  if(d>1){
    if(Hmtx && QR){
      mfit <- BPST.est.ho(B,Q2,K,lambda,Yi)
    }
    if(!Hmtx || !QR){
      mfit <- plagCV(B,H,K,lambda,Yi)
    }
  }
  if(d== -1){
    mfit <- BPST.est.pc(B,Yi)
  }
  mfit$mse <- mean((Yi-mfit$Yhat)^2);
  Yhat <- rep(NA,n)
  Yhat[Ind.inside] <- mfit$Yhat
  mfit$Yhat <- Yhat
  mfit$V <- V; mfit$Tr <- Tr; mfit$d <- d; mfit$r <- r; 
  mfit$B <- B; mfit$Q2 <- Q2; mfit$K <- K; mfit$Ind.inside <- Ind.inside;
  mfit$tria.all <- tria.all; mfit$Y <- Y; mfit$Z <- Z;
  mfit$call <- this.call;
  class(mfit) <- "BPST"
  return(mfit)
}