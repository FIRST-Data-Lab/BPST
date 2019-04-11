#' Make predictions from a fitted BPST object.
#'
#' This function is used to make predictions of a fitted BPST object.
#'
#' @importFrom Matrix Matrix
#' @importFrom pracma isempty
#' @importFrom Rcpp evalCpp
#'
#' @param mfit Fitted ``BPST" object.
#' \cr
#' @param Zpred The cooridinates for prediction -- default is the observed coordinates, \code{Z}.
#' \cr
#' @return A vector of predicted values is returned.
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
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
#' gridpoints=data.BPST(Z.grid,V,Tr,func,sigma,2019)
#' Y.grid=gridpoints$Y; mu.grid=gridpoints$mu;
#' ind=gridpoints$ind; ind.grid=(1:n.grid)[ind==1];
#' # Simulation parameters
#' n=2000;
#' ind.sam=sort(sample(ind.grid,n))
#' Y=as.matrix(gridpoints$Y[ind.sam]); Z=as.matrix(gridpoints$Z[ind.sam,]);
#' mfit=fit.BPST(Y,Z,V,Tr,d,r,lambda=10^seq(-6,6,by=0.5))
#' rmse=sqrt(mean((Y-mfit$Yhat)^2,na.rm=TRUE))
#' mpred=predict(mfit,Z.grid)
#' rmspe=sqrt(mean((Y.grid-mpred$Ypred)^2,na.rm=TRUE))
#' cat("rmse =",rmse,"and rmspe =",rmspe,"\n")
#' @export
#'
predict.BPST <- function(mfit,Zpred=NULL){
  if(identical(Zpred,mfit$Z) | isempty(Zpred)){
    Ypred <- mfit$beta_hat
    Ind.inside <- mfit$Ind.inside
  }else{
    Bfull <- basis(mfit$V,mfit$Tr,mfit$d,mfit$r,Zpred,FALSE,FALSE,FALSE,FALSE)
    B <- Bfull$B
    Ind.inside <- Bfull$Ind.inside
    Ypred <- rep(NA,nrow(Zpred))
    Ypred[Ind.inside] <- B%*%mfit$gamma_hat
  }
  list(Ypred=Ypred,Ind.Inside=Ind.inside)
}