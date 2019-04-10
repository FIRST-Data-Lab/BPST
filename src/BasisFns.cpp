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
NumericMatrix HQblkdiag(NumericMatrix C, NumericVector cnt) {
  int n=C.nrow();
  int m=C.ncol();
  int k=cnt.length()-1;
  NumericMatrix D(n,k*m);
  for(int ii=0; ii<k; ii++){
    if(cnt(ii)<cnt(ii+1)){
      for(int jj=cnt(ii); jj<cnt(ii+1); jj++){
        for(int kk=0; kk<m; kk++){
          D(jj,(ii*m+kk))=C(jj,kk);
        }
      }
    }
  }
  return D;
}

// [[Rcpp::export]]
NumericMatrix HQbary(NumericMatrix V, NumericMatrix Tr,NumericVector xx, NumericVector yy) {
  Function solve("solve");
  Function kronecker("kronecker");

  int nT=Tr.nrow();
  NumericMatrix A=transpose(Tr);
  NumericMatrix B(nT*3,2);
  for(int ii=0; ii<nT; ii++){
    B((ii*3+0),_)=V((A(0,ii)-1),_);
    B((ii*3+1),_)=V((A(1,ii)-1),_);
    B((ii*3+2),_)=V((A(2,ii)-1),_);
  }
  NumericMatrix C(nT*3,3);
  std::fill(C(_,0).begin(),C(_,0).end(),1);
  C(_,1)=B(_,0); C(_,2)=B(_,1);
  IntegerVector tmp=3*seq(0,nT);
  NumericVector cnt=as<NumericVector>(tmp);
  NumericMatrix D=transpose(HQblkdiag(C,cnt));

  NumericMatrix tmp1(3,xx.length());
  std::fill(tmp1(0,_).begin(),tmp1(0,_).end(),1);
  tmp1(1,_)=xx; tmp1(2,_)=yy;
  NumericMatrix tmp2(nT,1);
  std::fill(tmp2.begin(),tmp2.end(),1);
  NumericMatrix W=kronecker(tmp2,tmp1);
  NumericMatrix Lam=solve(D,W);
  return Lam;
}

// [[Rcpp::export]]
List HQgetInd(NumericMatrix V, NumericMatrix Tr, NumericVector xx, NumericVector yy) {
  Function which("which");

  int nT=Tr.nrow();
  NumericVector cnt(nT+1);
  NumericVector Ind;
  NumericVector tmp1, tmp2, tmp3;
  const double tol=-1e-12;
  NumericMatrix L=HQbary(V,Tr,xx,yy);
  
  for(int ii=0; ii<nT; ii++){
    NumericVector I=which(L((ii*3),_)>=tol & L((ii*3+1),_)>=tol & L((ii*3+2),_)>=tol);
    for(int jj=0; jj<I.length(); jj++){
      tmp1.push_back(L(ii*3,I(jj)-1));
      tmp2.push_back(L(ii*3+1,I(jj)-1));
      tmp3.push_back(L(ii*3+2,I(jj)-1));
      Ind.push_back(I(jj));
    }
    cnt(ii+1)=cnt(ii)+I.length();
  }
  NumericMatrix Lam(Ind.length(),3);
  Lam(_,0)=tmp1; Lam(_,1)=tmp2; Lam(_,2)=tmp3;
  return List::create(Named("Ind")=Ind,
                      Named("cnt")=cnt,
                      Named("Lam")=Lam);
}

// [[Rcpp::export]]
NumericMatrix mtxrbind(NumericMatrix a, NumericMatrix b) {
  int nra=a.nrow();
  int nrb=b.nrow();
  NumericMatrix out(nra+nrb,std::max(a.ncol(),b.ncol()));
  for(int ii=0; ii<nra+nrb; ii++){
    if(ii<nra){
      out(ii,_)=a(ii,_);
    }else{
      out(ii,_)=b(ii-nra,_);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix mtxcbind(NumericMatrix a, NumericMatrix b) {
  int nca=a.ncol();
  int ncb=b.ncol();
  NumericMatrix out(std::max(a.nrow(),b.nrow()),nca+ncb);
  for(int ii=0; ii<nca+ncb; ii++){
    if(ii<nca){
      out(_,ii)=a(_,ii);
    }else{
      out(_,ii)=b(_,ii-nca);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector vbind(NumericVector a, NumericVector b) {
  int na=a.length();
  int nb=b.length();
  NumericVector out(na+nb);
  if(na>0){
    IntegerVector inda=seq(0,(na-1));
    out[inda]=a;
  }
  if(nb>0){
    IntegerVector indb=seq(na,(na+nb-1));
    out[indb]=b;
  }
  return out;
}

// [[Rcpp::export]]
List BSpline2(NumericMatrix V, NumericMatrix Tr, double d, int r, NumericVector xx, NumericVector yy) {
  Function seval("seval");
  Function kronecker("kronecker");

  int nT=Tr.nrow();
  double m=(d+2)*(d+1)/2;
  int nx=xx.length();
  IntegerVector indn=seq(1,nx);
  NumericVector tmp1(nT);
  std::fill(tmp1.begin(),tmp1.end(),1);
  NumericMatrix tmp2=NumericMatrix::diag(m,1);
  double sfold=100.0; int nfold=ceil(nx/sfold);
  NumericVector xxi, yyi;
  IntegerVector indi;
  NumericVector Ind;
  NumericMatrix Bi;
  NumericVector tmp5, tmp6;

  for(int ifold=0; ifold<nfold; ifold++){
    if(ifold<(nfold-1)){
      indi=seq(ifold*sfold,(ifold+1)*sfold-1);
    }
    if(ifold==(nfold-1)){
      indi=seq(ifold*sfold,(nx-1));
    }
    tmp5=as<NumericVector>(indi);

    xxi=xx[indi]; yyi=yy[indi];
    List HQout=HQgetInd(V,Tr,xxi,yyi);
    NumericVector Indifold=as<NumericVector>(HQout["Ind"]);
    NumericVector cnti=as<NumericVector>(HQout["cnt"]);
    NumericMatrix Lami=as<NumericMatrix>(HQout["Lam"]);
    tmp6=tmp5[Indifold-1];
    tmp6=tmp6+1;

    NumericMatrix c=kronecker(tmp1,tmp2);
    int szi=Indifold.length();
    NumericMatrix Ai(szi,m);
    NumericMatrix Bifold;
    IntegerVector IndPi=seq(1,szi);
    for(int ii=0; ii<m; ii++){
      NumericVector tmp4=c(_,ii);
      NumericVector tmp3=seval(V,Tr,tmp4,szi,IndPi,cnti,Lami);
      Ai(_,ii)=tmp3;
    }
    Bifold=HQblkdiag(Ai,cnti);
    Bi=mtxrbind(Bi,Bifold);
    Ind=vbind(Ind,tmp6);
    // for(int jj=0; jj<Indifold.length(); jj++){
    //   Ind.push_back(indi(Indifold(jj)-1)+1);
    // }
  }
  return List::create(Named("nfold")=nfold,
                      Named("Bi")=Bi,
                      Named("Ind")=Ind);
}

// [[Rcpp::export]]
List BSpline(NumericMatrix V, NumericMatrix Tr, double d, int r, NumericVector xx, NumericVector yy) {
  Function seval("seval");
  Function kronecker("kronecker");

  int nT=Tr.nrow();
  double m=(d+2)*(d+1)/2;

  List HQout=HQgetInd(V,Tr,xx,yy);
  NumericVector Ind=as<NumericVector>(HQout["Ind"]);
  NumericVector cnt=as<NumericVector>(HQout["cnt"]);
  NumericMatrix Lam=as<NumericMatrix>(HQout["Lam"]);
  int sz=Ind.length();
  
  NumericVector tmp1(nT);
  std::fill(tmp1.begin(),tmp1.end(),1);
  NumericMatrix tmp2=NumericMatrix::diag(m,1);
  NumericMatrix c=kronecker(tmp1,tmp2);
  
  NumericMatrix A(sz,m);
  NumericMatrix B;
  if(sz>0){
    IntegerVector IndP=seq(1,sz);
    for(int ii=0; ii<m; ii++){
      NumericVector tmp4=c(_,ii);
      NumericVector tmp3=seval(V,Tr,tmp4,sz,IndP,cnt,Lam);
      A(_,ii)=tmp3;
    }
    B=HQblkdiag(A,cnt);
  }
  
  return List::create(Named("Bi")=B,
                      Named("Ind")=Ind);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.

/*** R
# data(V220); data(Tr220); data(Zi);
# bs=BSpline(V,Tr,4,1,Zi[,1],Zi[,2])
*/