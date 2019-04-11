#' Decide whether a point is inside of a given triangulation.
#'
#' This function is used to decided whether a point is inside of a given triangulation.
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @param V0 The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr0 The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param xx The x-cooridinate of points of dimension \code{n} by one. 
#' \cr
#' @param yy The y-cooridinate of points of dimension \code{n} by one.
#' \cr
#' @return A list of vectors, including:
#' \item{ind}{A vector of dimension \code{n} by one matrix that lists whether the points are inside of a given triangulation. 0 -- represents outside the triangulation, while 1 -- represents inside the triangulation.}
#' \item{ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' 
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @examples
#' xx=c(-0.25,0.75,0.25,1.25)
#' yy=c(-0.25,0.25,0.75,1.25)
#' V0=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
#' Tr0=rbind(c(1,2,3),c(1,3,4))
#' inVT(V0,Tr0,xx,yy)
#' 
#' @export

inVT <- function(V0,Tr0,xx,yy){
  xx <- as.vector(xx); yy <- as.vector(yy);
  n1 <- length(xx); n2 <- length(yy);
  if(n1!=n2) warning("The length of the coordinates do not match!")
  n <- min(n1,n2)
  ind <- as.numeric(insideVT(V0,Tr0,xx,yy))
  ind.inside <- (1:n)[ind==1]
  list(ind=ind,ind.inside=ind.inside)
}