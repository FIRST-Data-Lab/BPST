#' Cross-validation using Bivariate Penalized Spline over Triangulation
#'
#' This function implements k-fold cross-validation via bivariate penlized spline over triangulation, and returns the mean squared prediction error.
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
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' 
#' @param Hmtx The indicator of whether the smoothness matrix \code{H} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param Kmtx The indicator of whether the energy matrix \code{K} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param QR The indicator of whether a QR decomposition need to be performed on the smoothness matrix -- default is \code{TRUE}.
#' \cr
#' @param TA The indicator of whether the area of the triangles need to be calculated -- default is \code{TRUE}.
#' \cr
#' @return 
#' \item{lamc}{The tuning parameter selected by k-fold cross validation (CV).}
#' \item{mspe}{The mean squared prediction error calculated by k-fold cross validation (CV).}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @examples
#' # Triangulation
#' # Option 1;
#' data(V1); data(Tr1); d=5; r=1; V=V1; Tr=Tr1;
#' # Option 2
#' # data(V2); data(Tr2); d=5; r=1; V=V2; Tr=Tr2;
#' d=-1; r=1;
#' # Grid Points
#' n1.grid=101; n2.grid=101; n.grid=n1.grid*n2.grid;
#' u.grid=seq(0,1,length.out=n1.grid)
#' v.grid=seq(0,1,length.out=n2.grid)
#' uu.grid=rep(u.grid,each=n2.grid)
#' vv.grid=rep(v.grid,times=n1.grid)
#' Z.grid=as.matrix(cbind(uu.grid,vv.grid))
#' func=1; sigma=0.1;
#' gridpoints=data.BPST(Z.grid,V,Tr,func,sigma,2019)
#' Y.grid=gridpoints$Y; mu.grid=gridpoints$mu;
#' ind=gridpoints$ind; ind.grid=(1:n.grid)[ind==1];
#' # Simulation parameters
#' n=2000;
#' ind.sam=sort(sample(ind.grid,n))
#' Y=as.matrix(gridpoints$Y[ind.sam]); Z=as.matrix(gridpoints$Z[ind.sam,]);
#' cv.BPST(Y,Z,V,Tr,d,r,lambda=10^seq(-6,6,by=0.5),nfold=10)
#' @export

cv.BPST <- function(Y,Z,V,Tr,d=5,r=1,lambda=10^seq(-6,6,by=0.5),nfold=10,
                     Hmtx=TRUE,Kmtx=TRUE,QR=TRUE,TA=TRUE){
  if(nfold<3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold=10
  }
  if(!(is.matrix(Y) & ncol(Y)==1) & !is.vector(Y)){
    warning("The response variable, Y, should be a vector.")
    Y=as.vector(Y)
  }
  if(!is.matrix(Z)){
    warning("The coordinates of each point, Z, should be a n by two matrix.")
    Z=matrix(Z,ncol=2)
  }
  if(!is.matrix(V) | ncol(V)!=2){
    warning("The vertices matrix, V, should be a N by two matrix, with each row represents the coordinates of a vertex.")
    V=matrix(V,ncol=2)
  }
  if(!is.matrix(Tr) | ncol(Tr)!=3){
    warning("The triangulation matrix, Tr, should be a nT by three integer matrix, with each row represents the indexes of vertices of a triangle.")
    Tr=matrix(Tr,ncol=3)
  }
  n <- length(Y)
  ind.nna <- (1:n)[!is.na(Y)]
  Z <- Z[ind.nna,]
  Y <- Y[ind.nna]
  Bfull <- basis(V,Tr,d,r,Z,Hmtx,Kmtx,QR,TA)
  B <- Bfull$B
  H <- Bfull$H
  Q2 <- Bfull$Q2
  K <- Bfull$K
  Ind.inside <- Bfull$Ind.inside
  tria.all <- Bfull$tria.all
  Y <- Y[Ind.inside]
  n <- length(Y)
  
  if(d>1 & QR==TRUE){
    W <- as.matrix(B%*%Q2)
    D <- crossprod(t(crossprod(Q2,as.matrix(K))),Q2)
    D <- as.matrix(D)
    nl <- length(lambda)
  }
  if(d== -1){W <- as.matrix(B)}
  
  # k-fold Cross Validation
  sfold <- round(n/nfold); 
  indp <- sample(n,n,replace=FALSE);
  sspe_all <- c();
  for(ii in 1:nfold){
    sspei <- c()
    if(ii<nfold){
      indi <- sort(indp[((ii-1)*sfold+1):(ii*sfold)])
    }
    if(ii==nfold){
      indi <- sort(indp[((ii-1)*sfold+1):n])
    }
    indn <- setdiff(1:n,indi)
    
    Wi <- W[indi,]; Yi <- Y[indi];
    Wn <- W[indn,]; Yn <- Y[indn];
    WW <- crossprod(Wn,Wn)
    rhs <- crossprod(Wn,Yn)
    if(d>1){
      for(il in 1:nl){
        Lam <- lambda[il]
        Dlam <- Lam*D
        lhs <- WW+Dlam
        theta <- solve(lhs,rhs)
        gamma <- crossprod(t(Q2),theta)
        Yihat <- crossprod(t(Wi),theta)
        sspei <- c(sspei,sum((Yi-Yihat)^2))
      }
    }
    if(d== -1){
      lhs <- WW
      theta <- solve(lhs,rhs)
      Yihat <- crossprod(t(Wi),theta)
      sspei <- c(ssei,sum((Yi-Yihat)^2))
    }
    sspe_all <- rbind(sspe_all,sspei)
  } 
  sspe <- apply(sspe_all,2,sum)
  j <- which.min(sspe)
  if(d>1){lamc <- lambda[j]}
  if(d== -1){lamc <- NA}
  mspe <- min(sspe)/n
  list(lamc=lamc,mspe=mspe)
}