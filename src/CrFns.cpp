#include <Rcpp.h>
#include <math.h>
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

// Maybe not necessary
// // [[Rcpp::export]]
// NumericVector vbind(NumericVector a, NumericVector b) {
//   int na=a.length();
//   int nb=b.length();
//   NumericVector out(na+nb);
//   if(na>0){
//     IntegerVector inda=seq(0,(na-1));
//     out[inda]=a;
//   }
//   if(nb>0){
//     IntegerVector indb=seq(na,(na+nb-1));
//     out[indb]=b;
//   }
//   return out;
// }
// 
// // [[Rcpp::export]]
// NumericMatrix mtxcbind(NumericMatrix a, NumericMatrix b) {
//   int nca=a.ncol();
//   int ncb=b.ncol();
//   NumericMatrix out(std::max(a.nrow(),b.nrow()),nca+ncb);
//   for(int ii=0; ii<nca+ncb; ii++){
//     if(ii<nca){
//       out(_,ii)=a(_,ii);
//     }else{
//       out(_,ii)=b(_,ii-nca);
//     }
//   }
//   return out;
// }

//////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List CrIndices(int d, int r) {
  Function mtxcbind("mtxcbind") ;
  
  NumericMatrix I1;
  int start=1;
  int tmp1=d+1;
  for(int ii=0; ii<=r; ii++){
    for(int jj=0; jj<=(r-ii); jj++){
      IntegerVector newcol=seq((start+jj),(start+jj+d-r));
      NumericMatrix tmp2(newcol.length(),1);
      tmp2(_,0)=newcol;
      I1=mtxcbind(I1,tmp2);
    }
    start=start+tmp1;
    tmp1--;
  }
  int tmp3=-r*r/2+r*(d+3/2)+1; int tmp4=-r*r/2+r*(d+1/2)+d+1;
  IntegerVector tmp5=seq(tmp3,tmp4);
  NumericMatrix I2(tmp5.length(),1);
  I2(_,0)=as<NumericVector>(tmp5);
  return List::create(Named("I1")=I1,
                      Named("I2")=I2);
}

// [[Rcpp::export]]
List CrCellArrays(int d, int r) {
  List Jfit=CrIndices(d,r);
  NumericMatrix J1=as<NumericMatrix>(Jfit["I1"]);
  NumericMatrix J2=as<NumericMatrix>(Jfit["I2"]);
  NumericVector D1(d+1), D2(d+1);
  int s1=d+1; int s2=1;
  for(int ii=0; ii<(d+1); ii++){
    D1(ii)=s1;
    s1=s1+d-ii;
    D2(ii)=s2;
    s2=s2+d+1-ii;
  }
  NumericMatrix I21=clone<NumericMatrix>(J2);
  std::reverse(I21.begin(),I21.end());
  NumericVector tmp1=D1-r;
  IntegerVector tmp2=seq(0,d-r);
  NumericVector I22vec=tmp1[tmp2];
  std::reverse(I22vec.begin(),I22vec.end());
  NumericMatrix I22(I22vec.length(),1);
  I22(_,0)=I22vec;
  tmp1=D2+r;
  NumericVector I23vec=tmp1[tmp2];
  NumericMatrix I23(I23vec.length(),1);
  I23(_,0)=I23vec;
  
  NumericMatrix I11=J1;
  NumericMatrix tmp3(J1.nrow(),J1.ncol()), tmp6(J1.nrow(),J1.ncol());
  for(int ii=0; ii<(r+1); ii++){
    IntegerVector tmp4=seq(ii,d-r+ii);
    NumericVector tmp5=D1[tmp4];
    NumericVector tmp7=D2[tmp4];
    tmp3(_,ii)=tmp5;
    tmp6(_,ii)=tmp7;
  }
  int loc=r+2; int back=r+1;
  if(r>=1){
    for(int ii=1; ii<(r+1); ii++){
      for(int jj=0; jj<=(r-ii); jj++){
        tmp3(_,(loc-1))=tmp3(_,(loc-back-1))-1;
        tmp6(_,(loc-1))=tmp6(_,(loc-back-1))+1;
        loc++;
      }
      back--;
    }
  }
  NumericMatrix I12=tmp3;
  NumericMatrix I13(J1.nrow(),J1.ncol());
  for(int ii=0; ii<I13.ncol(); ii++){
    NumericVector tmp8=tmp6(_,ii);
    std::reverse(tmp8.begin(),tmp8.end());
    I13(_,ii)=tmp8;
  }
  List I1=List::create(I11,I12,I13);
  List I2=List::create(I21,I22,I23);
  return List::create(Named("I1")=I1,
                      Named("I2")=I2);
}

// [[Rcpp::export]]
List CrArrays(int d, int r) {
  List I1(r+1), I2(r+1);
  for(int ii=0; ii<=r; ii++){
    List Jfit=CrCellArrays(d,ii);
    List J1=as<List>(Jfit["I1"]);
    List J2=as<List>(Jfit["I2"]);
    I1[ii]=J1;
    I2[ii]=J2;
  }
  
  return List::create(Named("I1")=I1,
                      Named("I2")=I2);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# # Test for "CrIndices"
# CrIndices(d=2,r=0)
# CrIndices(d=2,r=1)

# # Test for "CrCellArray"
# CrCellArrays(2,0)
# CrCellArrays(2,1)

# # Test for "CrArray"
# CrArrays(2,0)
# CrArrays(2,1)
*/
