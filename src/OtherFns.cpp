#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix mvrbind(NumericMatrix a, NumericVector b) {
  int nra=a.nrow();
  int nca=a.ncol();
  int nb=b.length();
  NumericMatrix out((nra+1),std::max(nca,nb));
  for(int ii=0; ii<(nra+1); ii++){
    if(ii<nra){
      out(ii,_)=a(ii,_);
    }else{
      out(ii,_)=b;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix mvcbind(NumericMatrix a, NumericVector b) {
  int nra=a.nrow();
  int nca=a.ncol();
  int nb=b.length();
  NumericMatrix out(std::max(nra,nb),(nca+1));
  for(int ii=0; ii<(nca+1); ii++){
    if(ii<nca){
      out(_,ii)=a(_,ii);
    }else{
      out(_,ii)=b;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix newcol(Rcpp::Nullable<Rcpp::NumericMatrix> B0, Rcpp::Nullable<Rcpp::NumericVector> c0) {
  Function mvrbind("mvrbind");
  Function mvcbind("mvcbind"); 
  NumericMatrix A;
  if(B0.isNotNull() && c0.isNotNull()){
    NumericMatrix B(B0.get());
    NumericVector c(c0.get());
    int nrB=B.nrow();
    int ncB=B.ncol();
    int nc=c.length();
    A=B;
    if(nc>nrB){
      for(int ii=0; ii<(nc-nrB); ii++){
        NumericVector tmp1(ncB);
        std::fill(tmp1.begin(),tmp1.end(),0);
        A=mvrbind(A,tmp1);
      }
      A=mvcbind(A,c);
    }
    if(nc==nrB){
      A=mvcbind(A,c);
    }
    if(nc<nrB){
      NumericVector tmp2(nrB);
      std::fill(tmp2.begin(),tmp2.end(),0);
      IntegerVector ind=seq(0,(nc-1));
      tmp2[ind]=c;
      A=mvcbind(A,tmp2);
    }  
  }
  
  if(B0.isNull() && c0.isNotNull()){
    NumericVector c(c0.get());
    int nc=c.length();
    A=NumericMatrix(nc,1,c.begin());
  }
  
  if(B0.isNotNull() && c0.isNull()){
    NumericMatrix B(B0.get());
    A=B;
  }
  return A;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# B=rbind(c(1,3,5,7),c(2,4,6,8)); C=c(1:3);
# newcol(B,C)
# B=rbind(c(1,3,5,7),c(2,4,6,8)); C=c(1);
# newcol(B,C)
# B=rbind(c(1,3,5,7),c(2,4,6,8)); C=c(1:2);
# newcol(B,C)
# B=c(); C=c(1:3)
# newcol(B,C)
# B=rbind(c(1,3,5,7),c(2,4,6,8)); C=c();
# newcol(B,C)
# B=c(); C=c();
# newcol(B,C)
*/
