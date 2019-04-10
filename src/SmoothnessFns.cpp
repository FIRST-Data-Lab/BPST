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

// [[Rcpp::export]]
List tdata(NumericMatrix V, NumericMatrix Tr) {
  Function which("which");
  int m=Tr.nrow();
  int n=Tr.ncol();
  int numEdges=0;
  NumericVector Ti(3);
  NumericVector edge1(2), edge2(2), edge3(2);
  NumericVector ind1, ind2, ind3;
  NumericMatrix tmp1(3*m,2);
  NumericMatrix tmp2(m,3*m);
  for(int ii=0; ii<m; ii++){
    Ti=Tr(ii,_);
    edge1(0)=Ti(0); edge1(1)=Ti(1); std::sort(edge1.begin(),edge1.end());
    edge2(0)=Ti(1); edge2(1)=Ti(2); std::sort(edge2.begin(),edge2.end());
    edge3(0)=Ti(0); edge3(1)=Ti(2); std::sort(edge3.begin(),edge3.end());
    ind1=which(edge1(0)==tmp1(_,0) & edge1(1)==tmp1(_,1));
    ind2=which(edge2(0)==tmp1(_,0) & edge2(1)==tmp1(_,1));
    ind3=which(edge3(0)==tmp1(_,0) & edge3(1)==tmp1(_,1));
    if(ind1.length()==0){tmp1(numEdges,_)=edge1; tmp2(ii,numEdges)=1; numEdges++;}
    else{tmp2(ii,(ind1(0)-1))=1;}
    if(ind2.length()==0){tmp1(numEdges,_)=edge2; tmp2(ii,numEdges)=1; numEdges++;}
    else{tmp2(ii,(ind2(0)-1))=1;}
    if(ind3.length()==0){tmp1(numEdges,_)=edge3; tmp2(ii,numEdges)=1; numEdges++;}
    else{tmp2(ii,(ind3(0)-1))=1;}
  }
  
  NumericMatrix E(numEdges,2), TE(m,numEdges);
  // IntegerVector ind=seq(0,(numEdges-1));
  for(int ii=0; ii<2; ii++){
    E(_,ii)=tmp1(_,ii);
  }
  for(int ii=0; ii<m; ii++){
    TE(ii,_)=tmp2(ii,_);
  }
  
  int nV=V.nrow();
  NumericMatrix TV(m,nV);
  for(int ii=0; ii<m; ii++){
    for(int jj=0; jj<3; jj++){
      int ind=Tr(ii,jj)-1;
      TV(ii,ind)=1;
    }
  }
  
  NumericMatrix EV(numEdges,nV);
  for(int ii=0; ii<numEdges; ii++){
    for(int jj=0; jj<2; jj++){
      EV(ii,(E(ii,jj)-1))=1;
    }
  }
  return List::create(Named("E")=E,
                      Named("TE")=TE,
                      Named("TV")=TV,
                      Named("EV")=EV);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
