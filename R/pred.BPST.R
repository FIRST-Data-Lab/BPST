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
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia ang Li Wang from the Iowa State University.
#'
#' @export
#'
pred.BPST <- function(mfit,Zpred=NULL){
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