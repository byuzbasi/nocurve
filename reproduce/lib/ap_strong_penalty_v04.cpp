// Same fixed-penalty normal system, evaluated with two-part binary64 numbers.
// Acceptance is independently checked with exact integer arithmetic in Python.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include <cfenv>
#include <cmath>
#include <limits>
#ifdef __FAST_MATH__
#error This numerical kernel must not use fast-math.
#endif
namespace ap4 {
struct D {
  double h,l;
  D(double a=0,double b=0):h(a),l(b){}
};
D sum(double a,double b) {
  double s=a+b,v=s-a;
  return D(s,(a-(s-v))+(b-v));
}
D add(D a,D b) {
  D s=sum(a.h,b.h),t=sum(a.l,b.l),u=sum(s.l,t.h),v=sum(s.h,u.h);
  return sum(v.h,(v.l+u.l)+t.l);
}
D neg(D a){return D(-a.h,-a.l);}
D sub(D a,D b){return add(a,neg(b));}
D mul(D a,D b) {
  double p=a.h*b.h,e=std::fma(a.h,b.h,-p);
  e+=(a.h*b.l+a.l*b.h);e+=a.l*b.l;
  return sum(p,e);
}
D div(D a,D b) {
  if(b.h==0)Rcpp::stop("Zero pivot in precision solver");
  double q=a.h/b.h;D r=sub(a,mul(b,D(q)));
  double t=r.h/b.h;r=sub(r,mul(b,D(t)));
  return add(add(D(q),D(t)),D(r.h/b.h));
}
D root(D a) {
  if(!(a.h>0))Rcpp::stop("Nonpositive pivot; no penalty clipping or ridge added");
  D q(std::sqrt(a.h));
  for(int j=0;j<2;++j)q=mul(D(.5),add(q,div(a,q)));
  return q;
}
double value(D a){return a.h+a.l;}
}

// [[Rcpp::export]]
Rcpp::List ap4_solve_cpp(const arma::mat& X,const arma::mat& S,
                       const arma::mat& G,const arma::vec& rhs,double c) {
  using namespace ap4;
  static_assert(std::numeric_limits<double>::digits==53,"IEEE binary64 required");
  if(std::fegetround()!=FE_TONEAREST)Rcpp::stop("Round-to-nearest required");
  unsigned p=X.n_cols,n=X.n_rows;
  if(!p||p>96||n>256||G.n_rows!=p||G.n_cols!=p||S.n_rows!=p||S.n_cols!=p||
     rhs.n_elem!=p||!X.is_finite()||!G.is_finite()||!S.is_finite()||
     !rhs.is_finite()||!std::isfinite(c)||c<=0)Rcpp::stop("Invalid precision system");
  if(!arma::approx_equal(G,G.t(),"absdiff",0))Rcpp::stop("Asymmetric Gram matrix");
  // b' S b == b' ((S+S')/2) b exactly. Average the original entries in
  // two-part arithmetic; do not overwrite S, clip eigenvalues, or add a ridge.
  auto penalty=[&](unsigned i,unsigned j){return mul(D(.5),add(D(S(i,j)),D(S(j,i))));};
  std::vector<D> scale(p),L(p*p),z(p),u(p),b(p);
  for(unsigned i=0;i<p;++i)scale[i]=root(add(D(G(i,i)),mul(D(c),penalty(i,i))));
  for(unsigned i=0;i<p;++i)for(unsigned j=0;j<=i;++j) {
    D a=div(add(D(G(i,j)),mul(D(c),penalty(i,j))),mul(scale[i],scale[j]));
    for(unsigned k=0;k<j;++k)a=sub(a,mul(L[i*p+k],L[j*p+k]));
    L[i*p+j]=(i==j)?root(a):div(a,L[j*p+j]);
  }
  for(unsigned i=0;i<p;++i) {
    D a=div(D(rhs[i]),scale[i]);
    for(unsigned j=0;j<i;++j)a=sub(a,mul(L[i*p+j],z[j]));
    z[i]=div(a,L[i*p+i]);
  }
  for(int i=int(p)-1;i>=0;--i) {
    D a=z[i];for(unsigned j=i+1;j<p;++j)a=sub(a,mul(L[j*p+i],u[j]));
    u[i]=div(a,L[i*p+i]);b[i]=div(u[i],scale[i]);
  }
  Rcpp::NumericVector hi(p),lo(p),prediction(n);
  for(unsigned j=0;j<p;++j){hi[j]=b[j].h;lo[j]=b[j].l;}
  for(unsigned i=0;i<n;++i) {
    D a;for(unsigned j=0;j<p;++j)a=add(a,mul(D(X(i,j)),b[j]));
    prediction[i]=value(a);
  }
  if(!Rcpp::is_true(Rcpp::all(Rcpp::is_finite(hi)))||
     !Rcpp::is_true(Rcpp::all(Rcpp::is_finite(lo)))||
     !Rcpp::is_true(Rcpp::all(Rcpp::is_finite(prediction))))Rcpp::stop("Nonfinite precision output");
  return Rcpp::List::create(Rcpp::_["coefficients"]=hi,
    Rcpp::_["coefficient_correction"]=lo,Rcpp::_["prediction"]=prediction,
    Rcpp::_["representation"]="exact sum of two stored binary64 parts",
    Rcpp::_["penalty_representation"]="exact symmetric part of the same quadratic form",
    Rcpp::_["diagonal_equilibration"]=true,Rcpp::_["ridge_added"]=0,
    Rcpp::_["penalty_clipping"]=false);
}
