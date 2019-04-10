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
double degree(double m) {
  double d;
  d=(-3+sqrt(8*m+1))/2;
  return d;
}

// [[Rcpp::export]]
List indices(int d) {
  int m=(d+1)*(d+2)/2;
  int Mark=1;
  NumericVector I(m);
  NumericVector J(m);
  NumericVector K(m);
  for(int i=d;i>=0;i--){
    for(int j=0;j<=i;j++){
      I(Mark+j-1)=(i-j);
      J(Mark+j-1)=j;
      K(Mark+j-1)=(d-i);
    }
    Mark=Mark+i+1;
  }
  return List::create(Named("I")=I,
                      Named("J")=J,
                      Named("K")=K);
}

// [[Rcpp::export]]
List bary(NumericVector V1, NumericVector V2, NumericVector V3, NumericVector X, NumericVector Y) {
  Function solve("solve");
  
  int px;
  px=X.length();
  NumericVector tmp1(3);
  NumericVector tmp2(3);
  NumericVector tmp3(3);
  NumericVector One(px);
  NumericMatrix A(3,3);
  NumericMatrix c(3,px);
  
  std::fill(tmp1.begin(),tmp1.end(),1);
  tmp2(0)=V1(0); tmp2(1)=V2(0); tmp2(2)=V3(0);
  tmp3(0)=V1(1); tmp3(1)=V2(1); tmp3(2)=V3(1);
  A(0,_)=tmp1;
  A(1,_)=tmp2;
  A(2,_)=tmp3;
  
  std::fill(One.begin(),One.end(),1);
  c(0,_)=One;
  c(1,_)=X;
  c(2,_)=Y;
  
  NumericMatrix lam=solve(A,c);
  NumericVector lam1=lam(0,_);
  NumericVector lam2=lam(1,_);
  NumericVector lam3=lam(2,_);
  return List::create(Named("lam1")=lam1,
                      Named("lam2")=lam2,
                      Named("lam3")=lam3);
}

// [[Rcpp::export]]
NumericVector tcord(NumericVector V1, NumericVector V2, NumericVector V3, NumericVector u) {
  NumericVector u1(1), u2(1);
  std::fill(u1.begin(),u1.end(),u(0));
  std::fill(u2.begin(),u2.end(),u(1));
  List lamv=bary(V1,V2,V3,u1,u2);
  NumericVector w1(1), w2(1);
  std::fill(w1.begin(),w1.end(),0);
  std::fill(w2.begin(),w2.end(),0);
  List a=bary(V1,V2,V3,w1,w2);
  double lam1=as<double>(lamv["lam1"])-as<double>(a["lam1"]);
  double lam2=as<double>(lamv["lam2"])-as<double>(a["lam2"]);
  double lam3=as<double>(lamv["lam3"])-as<double>(a["lam3"]);
  NumericVector lam(3);
  lam(0)=lam1;
  lam(1)=lam2;
  lam(2)=lam3;
  return lam;
}

// [[Rcpp::export]]
double triarea(NumericVector V1, NumericVector V2, NumericVector V3) {
  double x, y, a, b, c, d, A;
  x=V1(0);
  y=V1(1);
  a=V2(0);
  b=V2(1);
  c=V3(0);
  d=V3(1);
  A=((a-x)*(d-y)-(c-x)*(b-y))/2;
  return A;
}

// [[Rcpp::export]]
NumericVector locate(NumericVector I1, NumericVector J1, NumericVector K1, NumericVector I, NumericVector J, NumericVector K) {
  Function match("match");
  NumericMatrix M1=cbind(I1,J1,K1);
  NumericMatrix M=cbind(I,J,K);
  NumericVector Index=match(DataFrame(transpose(M1)),DataFrame(transpose(M)));
  return Index;
}

// [[Rcpp::export]]
NumericMatrix dirder(NumericMatrix Bin, double lam1, double lam2, double lam3) {
  int m=Bin.nrow();
  double d=degree(m);
  List Index=indices(d);
  NumericVector I=as<NumericVector>(Index["I"]);
  NumericVector J=as<NumericVector>(Index["J"]);
  NumericVector K=as<NumericVector>(Index["K"]);
  List Index1=indices(d-1);
  NumericVector I1=as<NumericVector>(Index1["I"]);
  NumericVector J1=as<NumericVector>(Index1["J"]);
  NumericVector K1=as<NumericVector>(Index1["K"]);
  NumericVector indx1=locate(I1+1,J1,K1,I,J,K);
  NumericVector indx2=locate(I1,J1+1,K1,I,J,K);
  NumericVector indx3=locate(I1,J1,K1+1,I,J,K);
  int nr=indx1.length();
  int nc=Bin.ncol();
  NumericMatrix Bout(nr,nc);
  for(int ii=0; ii<nr; ii++){
    Bout(ii,_)=Bin((indx1(ii)-1),_)*lam1;
  }
  for(int ii=0; ii<nr; ii++){
    Bout(ii,_)=Bout(ii,_)+Bin((indx2(ii)-1),_)*lam2;
  }
  for(int ii=0; ii<nr; ii++){
    Bout(ii,_)=Bout(ii,_)+Bin((indx3(ii)-1),_)*lam3;
  }
  NumericMatrix DerBcoeff=Bout*d;
  return DerBcoeff;
}

// [[Rcpp::export]]
NumericMatrix build(double d) {
  List Index=indices(d);
  NumericVector I=as<NumericVector>(Index["I"]);
  NumericVector J=as<NumericVector>(Index["J"]);
  NumericVector K=as<NumericVector>(Index["K"]);
  double m=(d+1)*(d+2)/2;
  NumericMatrix Mat(m,m);
  for(int ii=0; ii<m; ii++){
    for(int jj=0; jj<m; jj++){
      Mat(jj,ii)=Rf_choose((I(ii)+I(jj)),I(ii))*Rf_choose((J(ii)+J(jj)),J(ii))*Rf_choose((K(ii)+K(jj)),K(ii));
    }
  }
  int tmp=Rf_choose(2*d,d)*Rf_choose((2*d+2),2);
  Mat=Mat/tmp;
  return Mat;
}

// [[Rcpp::export]]
NumericMatrix locEng(NumericVector V1, NumericVector V2, NumericVector V3, NumericMatrix Mat, double d) {
  // calling R function crossprod()
  Function crossprod("crossprod"); 
  
  double m=(d+1)*(d+2)/2;
  NumericMatrix Id=NumericMatrix::diag(m,1);
  NumericVector vx(2), vy(2);
  vx(0)=1; vx(1)=0;
  vy(0)=0; vy(1)=1;
  NumericVector lamx=tcord(V1,V2,V3,vx);
  NumericVector lamy=tcord(V1,V2,V3,vy);
  NumericMatrix Dx=dirder(Id,lamx(0),lamx(1),lamx(2));
  NumericMatrix Dxx=dirder(Dx,lamx(0),lamx(1),lamx(2));
  NumericMatrix Dxy=dirder(Dx,lamy(0),lamy(1),lamy(2));
  NumericMatrix Dy=dirder(Id,lamy(0),lamy(1),lamy(2));
  NumericMatrix Dyy=dirder(Dy,lamy(0),lamy(1),lamy(2));
  double tmp=std::abs(triarea(V1,V2,V3));
  NumericMatrix area=NumericMatrix::diag(m,tmp);
  NumericMatrix tmp1=crossprod(Dxx,crossprod(Mat,Dxx));
  NumericMatrix tmp2=crossprod(Dxy,crossprod(Mat,Dxy));
  NumericMatrix tmp3=crossprod(Dyy,crossprod(Mat,Dyy));
  int nr=tmp1.nrow();
  int nc=tmp1.ncol();
  NumericMatrix tmp4(nr,nc);
  for(int jj=0; jj<nc; jj++){
    tmp4(_,jj)=tmp1(_,jj)+tmp2(_,jj)*2+tmp3(_,jj);
  }
  NumericMatrix K=crossprod(area,tmp4);
  return K;
}

// [[Rcpp::export]]
NumericMatrix energy(NumericMatrix V, NumericMatrix Tr, double d) {
  int n=Tr.nrow();
  double m=(d+1)*(d+2)/2;
  NumericMatrix Mat=build(d-2);
  NumericMatrix K(n*m,n*m);
  
  for(int kk=0; kk<n ;kk++){
    NumericVector V1=V((Tr(kk,0)-1),_);
    NumericVector V2=V((Tr(kk,1)-1),_);
    NumericVector V3=V((Tr(kk,2)-1),_);
    NumericMatrix LocK=locEng(V1,V2,V3,Mat,d);
    
    for(int ii=0; ii<LocK.nrow(); ii++){
      for(int jj=0; jj<LocK.ncol(); jj++){
        if(LocK(ii,jj)!=0){
          int ind1=kk*m+jj;
          int ind2=kk*m+ii;
          K(ind1,ind2)=LocK(ii,jj);
        }
      }
    }
  }
  return K;
}

// [[Rcpp::export]]
LogicalVector insideVT(NumericMatrix V, NumericMatrix Tr, NumericVector xx, NumericVector yy) {
  int nT=Tr.nrow();
  double tol=-1e-12;
  LogicalVector ind;
  NumericVector VT,V1,V2,V3;
  NumericVector lam1,lam2,lam3;
  for(int ii=0; ii<nT; ii++){
    VT=Tr(ii,_);
    V1=V((VT(0)-1),_);
    V2=V((VT(1)-1),_);
    V3=V((VT(2)-1),_);
    List a=bary(V1,V2,V3,xx,yy);
    lam1=as<NumericVector>(a["lam1"]);
    lam2=as<NumericVector>(a["lam2"]);
    lam3=as<NumericVector>(a["lam3"]);
    LogicalVector tmp=(lam1>tol & lam2>tol & lam3>tol);
    if(ii==0){
      ind=tmp;
    }else{
      ind=(ind | tmp);
    }
  }
  return ind;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
# xx=c(-0.25,0.75,0.25,1.25)
# yy=c(-0.25,0.25,0.75,1.25)
# V0=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
# Tr0=rbind(c(1,2,3),c(1,3,4))
# ind=as.numeric(insideVT(V0,Tr0,xx,yy))
# ind
# (1:length(xx))[ind==1]
*/
