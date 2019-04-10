#' Model Fitting using Bivariate Penalized Spline over Triangulation
#'
#' This function conducts the model fitting via bivariate penlized spline over triangulation.
#'
#' @importFrom Matrix Matrix
#' 
#' @param Z The cooridinates of dimension \code{n} by two. Each row is the coordinates of a point.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param func The choice of test function -- default is 1. Possible choices include 1,2,3,4,5,6,7,8.
#' \cr
#' @param sigma The standard deviation of the white noise --  default is 0.1.
#' \cr 
#' @return A list of vectors and matrice, including:
#' \item{Y}{The response variable.}
#' \item{mu}{The mean function.}
#' \item{Z}{The coordinates.}
#' \item{ind}{A vector contains the indicators whether the point is inside the given triangulation.}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia ang Li Wang from the Iowa State University.
#' 
#' @export
DataGenerator <- function(Z,V,Tr,func=1,sigma=0.1,iter=2019){
  set.seed(iter)
  # location information
  if(isempty(Z)){
    stop("The location information is required for data generation.")
  }
  Z <- matrix(Z,ncol=2); 
  xx <- Z[,1]; yy <- Z[,2]; n <- nrow(Z);
  ind <- inVT(V,Tr,xx,yy)$ind
    
  # triangulation information
  if(isempty(V) | isempty(Tr)){
    stop("The triangulation information is required for data generation.")
  }
  V <- matrix(V,ncol=2); Tr <- matrix(Tr,ncol=3);
    
  # test functions
  if(!(func %in% (1:8))){
    stop("Test function can only be integers between 1 and 8.")
  }
  if(func==1){
    mu <- xx^2+3*yy^2+4*xx*yy # Quadratic
  }
  if(func==2){
    mu <- xx^3+yy^3 # Cubic
  }
  if(func==3){
    mu <- xx^4+yy^4 # Quadruplicate
  }
  if(func==4){
    mu <- sin(pi*xx)*sin(pi*yy) #Sine
  }
  if(func==5){
    mu <- (1-cos(2*pi*xx))*(1-cos(2*pi*yy)) # Cosine
  }
  if(func==6){
    mu <- exp(-50*((xx-0.5)^2+(yy-0.5)^2)) # Bump
  }
  if(func==7){
    mu <- 1/(1+exp(-10*(xx+yy)+10)) # Logit
  }
  if(func==8){
    mu <- atan((8*xx-4)^2-(8*yy-4)^2) # arctan
  }
    
  eps <- rnorm(n,mean=0,sd=sigma)
  mu[ind==0] <- NA; eps[ind==0] <- NA;
  Y <- mu+eps
  list(Y=Y,mu=mu,Z=Z,ind=ind)
}
