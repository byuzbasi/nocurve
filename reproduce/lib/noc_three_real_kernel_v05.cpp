// Same linearly binned, tau=4 Gaussian local-linear objective as KernSmooth.
// Only grid nodes consumed by the frozen linear readout are evaluated.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List n3r5_binned_nodes_cpp(Rcpp::IntegerVector occupied,
    Rcpp::NumericVector counts,Rcpp::NumericVector totals,
    Rcpp::IntegerVector wanted,double delta,double h,int radius) {
  int n=occupied.size(),m=wanted.size();
  if(counts.size()!=n||totals.size()!=n||n==0||h<=0||delta<=0)
    Rcpp::stop("Invalid binned kernel design");
  Rcpp::NumericVector value(m),normal(m),slope_out(m);Rcpp::IntegerVector rank(m);
  for(int q=0;q<m;++q) {
    int j=wanted[q];
    int lo=std::lower_bound(occupied.begin(),occupied.end(),j-radius)-occupied.begin();
    int hi=std::upper_bound(occupied.begin(),occupied.end(),j+radius)-occupied.begin();
    long double W=0,U=0,Y=0;
    for(int k=lo;k<hi;++k) {
      if(!R_finite(counts[k])||counts[k]<=0||!R_finite(totals[k]))
        Rcpp::stop("Invalid occupied kernel bin");
      double u=delta*(occupied[k]-j)/h,w=std::exp(-u*u/2);
      long double weight=(long double)counts[k]*w;
      W+=weight;U+=weight*u;Y+=(long double)totals[k]*w;
    }
    if(!(W>0))Rcpp::stop("Consumed grid node has no kernel mass");
    long double mean_u=U/W,mean_y=Y/W,V=0,C=0;
    for(int k=lo;k<hi;++k) {
      double u=delta*(occupied[k]-j)/h,w=std::exp(-u*u/2);
      long double weight=(long double)counts[k]*w,du=(long double)u-mean_u;
      V+=weight*du*du;
      C+=weight*du*((long double)totals[k]/counts[k]-mean_y);
    }
    long double slope=0;
    if(V>0){slope=C/V;rank[q]=2;}
    else if(mean_u==0){rank[q]=1;} // Intercept uniquely estimable; slope need not be.
    else Rcpp::stop("Consumed grid-node intercept is not identifiable");
    double intercept=(double)(mean_y-slope*mean_u),b=(double)slope;
    long double g0=0,g1=0,scale=1;
    for(int k=lo;k<hi;++k) {
      double u=delta*(occupied[k]-j)/h,w=std::exp(-u*u/2);
      long double weight=(long double)counts[k]*w,
        y=(long double)totals[k]/counts[k],p=(long double)intercept+b*u,r=y-p;
      g0+=weight*r;g1+=weight*u*r;scale+=weight*(std::fabs(y)+std::fabs(p));
    }
    double residual=(double)(std::max(std::fabs(g0),std::fabs(g1))/scale);
    if(!R_finite(intercept)||!R_finite(residual)||residual>1e-8)
      Rcpp::stop("Consumed kernel node failed normal-equation check at 1e-8");
    value[q]=intercept;normal[q]=residual;slope_out[q]=b;
  }
  return Rcpp::List::create(Rcpp::Named("prediction")=value,
    Rcpp::Named("normal_residual")=normal,Rcpp::Named("rank")=rank,
    Rcpp::Named("slope")=slope_out);
}
// Replay the normal equations using saved coefficients only; no fit/optimizer.
// [[Rcpp::export]]
Rcpp::NumericVector n3r5_binned_certificate_cpp(Rcpp::IntegerVector occupied,
    Rcpp::NumericVector counts,Rcpp::NumericVector totals,
    Rcpp::IntegerVector wanted,Rcpp::NumericVector intercept,Rcpp::NumericVector slope,
    double delta,double h,int radius) {
  int m=wanted.size();
  if(intercept.size()!=m||slope.size()!=m)Rcpp::stop("Certificate dimensions differ");
  Rcpp::NumericVector result(m);
  for(int q=0;q<m;++q) {
    int j=wanted[q];
    int lo=std::lower_bound(occupied.begin(),occupied.end(),j-radius)-occupied.begin();
    int hi=std::upper_bound(occupied.begin(),occupied.end(),j+radius)-occupied.begin();
    long double g0=0,g1=0,scale=1,W=0;
    for(int k=lo;k<hi;++k) {
      double u=delta*(occupied[k]-j)/h,w=std::exp(-u*u/2);
      long double weight=(long double)counts[k]*w,y=(long double)totals[k]/counts[k],
        p=(long double)intercept[q]+slope[q]*u,r=y-p;
      W+=weight;g0+=weight*r;g1+=weight*u*r;scale+=weight*(std::fabs(y)+std::fabs(p));
    }
    if(!(W>0))Rcpp::stop("Certificate node has no positive kernel mass");
    result[q]=(double)(std::max(std::fabs(g0),std::fabs(g1))/scale);
  }
  return result;
}
